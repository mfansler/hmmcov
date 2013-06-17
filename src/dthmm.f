      subroutine loop1(m, n, phi, prob, Pi, logalp, lscale, tmp)
c     first loop (forward eqns) in forwardback.dthmm
      implicit none
      integer i, j, m, n
      double precision phi(m), sumphi, lscale
      double precision prob(n,m*m), Pi(m,m), logalp(n,m)
      double precision tmp(m)
c     the above array occurs in the subroutine call for
c     memory allocation reasons in non gfortran compilers
c     its contents are purely internal to this subroutine
      lscale = 0
      i = 1
      do while(i .le. n)
          if (i .gt. 1) then 
						call multi1(m, phi, Pi, prob, i, n, tmp)
					else 
						phi(1) = phi(1)*prob(1, 1)
						phi(2) = phi(2)*prob(1, 3)
					endif 
				
          j = 1
          sumphi=0.0
          do while(j .le. m)
c              phi(j) = phi(j)*prob(i, j)
              sumphi = sumphi + phi(j)
              j = j+1
          enddo
          j = 1
          do while(j .le. m)
              phi(j) = phi(j)/sumphi
              j = j+1
          enddo
          lscale = lscale + dlog(sumphi)
          j = 1
          do while(j .le. m)
							if (phi(j) .eq. 0.0) phi(j) = 10**(-200)
              logalp(i,j) = dlog(phi(j)) + lscale
              j = j+1
          enddo
          i = i+1
      enddo
      end


      subroutine loop2(m, n, phi, prob, Pi, logbet, lscale, tmp)
c     second loop (backward eqns) in forwardback.dthmm
      implicit none
      integer i, j, m, n
      double precision phi(m), sumphi, lscale
      double precision prob(n,m*m), Pi(m,m), logbet(n,m)
      double precision tmp(m)
c     the above array occurs in the subroutine call for
c     memory allocation reasons in non gfortran compilers
c     its contents are purely internal to this subroutine
      i = n-1
      do while(i .ge. 1)
c          j = 1
c          do while(j .le. m)
c              phi(j) = phi(j)*prob(i+1, j)
c              j = j+1
c          enddo
          call multi2(m, Pi, phi, prob,i, n, tmp)
          j = 1
          sumphi=0.0
          do while(j .le. m)
							if (phi(j) .eq. 0.0) phi(j) = 10**(-200)
              logbet(i,j) = dlog(phi(j)) + lscale
              sumphi = sumphi + phi(j)
              j = j+1
          enddo
          j = 1
          do while(j .le. m)
              phi(j) = phi(j)/sumphi
              j = j+1
          enddo
          lscale = lscale + dlog(sumphi)
          i = i-1
      enddo
      end
