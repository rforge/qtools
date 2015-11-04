c  Define the function to be maximized : this is the maximum score
c  function for any quantile tau between 0 and 1.
      subroutine fcn(y,xx,nobs,k,b,h,tau,sgn)
      integer nobs, k
      double precision b(*), h, tau
      double precision y(nobs), xx(nobs,k), ypred, sgn

      h = 0.D0

      do 100, i = 1, nobs
       ypred = 0.D0
         do j=1, k-1
            ypred = ypred + (b(j)*xx(i,j))
         enddo
         ypred = ypred + sgn * xx(i,k)
        h = (1/dble(nobs)) * (y(i)-(1.D0-tau)) * (ypred/abs(ypred)) + h
100   continue

      return
      end


      subroutine sa(y,xx,nobs,k,b,tau,max,rt,eps,ns,nt,neps,maxevl,lb,
     1              ub,c,iseed1,iseed2,t,vm,bopt,fopt,nacc,nfcnev,
     2              nobds,ier,fstar,xp,nacp,aflag,sgn)

c  type all external variables.
      double precision  xx(nobs,k), y(*), b(k-1), lb(*), ub(*),
     1                  c(*), vm(*), fstar(*), bopt(k-1), xp(k-1), t,
     2                  eps, rt, fopt, tau, sgn
      integer  nacp(*), k, ns, nt, neps, nacc, maxevl,
     1         nobds, ier, nfcnev, iseed1, iseed2, aflag
      logical  max

c  type all internal variables.
      double precision  f, fp, p, pp, ratio
      integer  nup, ndown, nrej, nnew, lnobds, h, i, j, m, k, nobs
      logical  quit

c  type all functions.
      double precision  exprep
      real  ranmar

c  initialize the random number generator ranmar.
      call rmarin(iseed1,iseed2)

c  set initial values.
      nacc = 0
      nobds = 0
      nfcnev = 0
      ier = 99

      do 10, i = 1, k-1
         bopt(i) = b(i)
         nacp(i) = 0
10    continue

      do 20, i = 1, neps
         fstar(i) = 1.0d+20
20    continue 

c  evaluate the function with input b and return value as f.
      call fcn(y,xx,nobs,k,b,f,tau,sgn)

c  if the function is to be minimized, switch the sign of the function.
c  note that all intermediate and final output switches the sign back
c  to eliminate any possible confusion for the user.
      if(.not. max) f = -f
      nfcnev = nfcnev + 1
      fopt = f
      fstar(1) = f

c  start the main loop. note that it terminates if (i) the algorithm
c  succesfully optimizes the function or (ii) there are too many
c  function evaluations (more than maxevl).
100   nup = 0
      nrej = 0
      nnew = 0
      ndown = 0
      lnobds = 0

      do 400, m = 1, nt
         do 300, j = 1, ns
            do 200, h = 1, k-1

c  generate xp, the trial value of b. note use of vm to choose xp.
               do 110, i = 1, k-1
                  if (i .eq. h) then
                     xp(i) = b(i) + (ranmar()*2.- 1.) * vm(i)
                  else
                     xp(i) = b(i)
                  end if

c  if xp is out of bounds, select a point in bounds for the trial.
                  if((xp(i) .lt. lb(i)) .or. (xp(i) .gt. ub(i))) then
                    xp(i) = lb(i) + (ub(i) - lb(i))*ranmar()
                    lnobds = lnobds + 1
                    nobds = nobds + 1
                  end if
110            continue

c  evaluate the function with the trial point xp and return as fp.
               call fcn(y,xx,nobs,k,xp,fp,tau,sgn)
               if(.not. max) fp = -fp
               nfcnev = nfcnev + 1

c  if too many function evaluations occur, terminate the algorithm.
               if(nfcnev .ge. maxevl) then
                  aflag = -1
                  if (.not. max) fopt = -fopt
                  ier = 1
                  return
               end if

c  accept the new point if the function value increases.
               if(fp .ge. f) then
                  do 120, i = 1, k-1
                     b(i) = xp(i)
120               continue
                  f = fp
                  nacc = nacc + 1
                  nacp(h) = nacp(h) + 1
                  nup = nup + 1

c  if greater than any other point, record as new optimum.
                  if (fp .gt. fopt) then
                     do 130, i = 1, k-1
                        bopt(i) = xp(i)
130                  continue
                     fopt = fp
                     nnew = nnew + 1
                  end if

c  if the point is lower, use the metropolis criteria to decide on
c  acceptance or rejection.
               else
                  p = exprep((fp - f)/t)
                  pp = ranmar()
                  if (pp .lt. p) then
                     do 140, i = 1, k-1
                        b(i) = xp(i)
140                  continue
                     f = fp
                     nacc = nacc + 1
                     nacp(h) = nacp(h) + 1
                     ndown = ndown + 1
                  else
                     nrej = nrej + 1
cc                   if(iprint .ge. 3) call prt7(max)
                  end if
               end if

200         continue
300      continue

c  adjust vm so that approximately half of all evaluations are accepted.
         do 310, i = 1, k-1
            ratio = dfloat(nacp(i)) /dfloat(ns)
            if (ratio .gt. .6) then
               vm(i) = vm(i)*(1. + c(i)*(ratio - .6)/.4)
            else if (ratio .lt. .4) then
               vm(i) = vm(i)/(1. + c(i)*((.4 - ratio)/.4))
            end if
            if (vm(i) .gt. (ub(i)-lb(i))) then
               vm(i) = ub(i) - lb(i)
            end if
310      continue

         do 320, i = 1, k-1
            nacp(i) = 0
320      continue

400   continue


c  check termination criteria.
      quit = .false.
      fstar(1) = f
      if ((fopt - fstar(1)) .le. eps) quit = .true.
      do 410, i = 1, neps
         if (abs(f - fstar(i)) .gt. eps) quit = .false.
410   continue

c  terminate sa if appropriate.
      if (quit) then
         do 420, i = 1, k-1
            b(i) = bopt(i)
420      continue
         ier = 0
         if (.not. max) fopt = -fopt
cc       if(iprint .ge. 1) call prt10
         return
      end if

c  if termination criteria is not met, prepare for another loop.
      t = rt*t
      do 430, i = neps, 2, -1
         fstar(i) = fstar(i-1)
430   continue
      f = fopt
      do 440, i = 1, k-1
         b(i) = bopt(i)
440   continue

c  loop again.
      go to 100

      end

      function  exprep(rdum)
c  this function replaces exp to avoid under- and overflows and is
c  designed for ibm 370 type machines. it may be necessary to modify
c  it for other machines. note that the maximum and minimum values of
c  exprep are such that they has no effect on the algorithm.

      double precision  rdum, exprep

      if (rdum .gt. 174.) then
         exprep = 3.69d+75
      else if (rdum .lt. -180.) then
         exprep = 0.0
      else
         exprep = exp(rdum)
      end if

      return
      end

      subroutine rmarin(ij,kl)
c  this subroutine and the next function generate random numbers. see
c  the comments for sa for more information. the only changes from the
c  orginal code is that (1) the test to make sure that rmarin runs first
c  was taken out since sa assures that this is done (this test didn't
c  compile under ibm's vs fortran) and (2) typing ivec as integer was
c  taken out since ivec isn't used. with these exceptions, all following
c  lines are original.

c this is the initialization routine for the random number generator
c     ranmar()
c note: the seed variables can have values between:    0 <= ij <= 31328
c                                                      0 <= kl <= 30081
      real u(97), c, cd, cm
      integer i97, j97
      common /raset1/ u, c, cd, cm, i97, j97
      i = mod(ij/177, 177) + 2
      j = mod(ij    , 177) + 2
      k = mod(kl/169, 178) + 1
      l = mod(kl,     169)
      do 2 ii = 1, 97
         s = 0.0
         t = 0.5
         do 3 jj = 1, 24
            m = mod(mod(i*j, 179)*k, 179)
            i = j
            j = k
            k = m
            l = mod(53*l+1, 169)
            if (mod(l*m, 64) .ge. 32) then
               s = s + t
            endif
            t = 0.5 * t
3        continue
         u(ii) = s
2     continue
      c = 362436.0 / 16777216.0
      cd = 7654321.0 / 16777216.0
      cm = 16777213.0 /16777216.0
      i97 = 97
      j97 = 33
      return
      end

      function ranmar()
      real u(97), c, cd, cm
      integer i97, j97
      common /raset1/ u, c, cd, cm, i97, j97
         uni = u(i97) - u(j97)
         if( uni .lt. 0.0 ) uni = uni + 1.0
         u(i97) = uni
         i97 = i97 - 1
         if(i97 .eq. 0) i97 = 97
         j97 = j97 - 1
         if(j97 .eq. 0) j97 = 97
         c = c - cd
         if( c .lt. 0.0 ) c = c + cm
         uni = uni - c
         if( uni .lt. 0.0 ) uni = uni + 1.0
         ranmar = uni
      return
      end

