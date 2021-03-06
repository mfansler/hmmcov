      subroutine multi1(m, a, b, p, i, n, c)
c     a is row vector
c     b is (m*m) matrix
c     a is replaced by the matrix product of a*b
      implicit none
      integer m, j, k, i, n
      double precision a(m), b(m, m), c(m), p(n, m**2)
      j = 1
      do while(j .le. m)
          k = 1
          c(j) = 0.0
          do while(k .le. m)
              c(j) = c(j) + a(k)*b(k, j)*p(i, (j-1)*m+k)
              k = k+1
          enddo
          j = j+1
      enddo
      j = 1
      do while(j .le. m)
          a(j) = c(j)
          j = j+1
      enddo
      end


      subroutine multi2(m, a, b, p, i, n, c)
c     a is (m*m) matrix
c     b is column vector
c     b is replaced by the matrix product of a*b
      implicit none
      integer m, j, k, i, n
      double precision a(m, m), b(m), c(m), p(n,m**m)
      j = 1
      do while(j .le. m)
          k = 1
          c(j) = 0.0
          do while(k .le. m)
              c(j) = c(j) + a(j, k)*b(k)*p(i+1, (k-1)*m+j)
              k = k+1
          enddo
          j = j+1
      enddo
      j = 1
      do while(j .le. m)
          b(j) = c(j)
          j = j+1
      enddo
      end


      subroutine multi3(p, q, r, a, b, c)
c     a is (p*q) matrix
c     b is (q*r) matrix
c     c is (p*r) matrix
c     c is the matrix product of a*b
      implicit none
      integer i, j, k, p, q, r
      double precision a(p, q), b(q, r), c(p, r)
      i = 1
      do while(i .le. p)
          j = 1
          do while(j .le. r)
              k = 1
              c(i, j) = 0.0
              do while(k .le. q)
                  c(i, j) = c(i, j) + a(i, k)*b(k, j)
                  k = k+1
              enddo
              j = j+1
          enddo
          i = i+1
      enddo
      end


      subroutine multi4(m, a, b, c, d)
c     a is row vector
c     b is (m*m) matrix
c     c is (m*m) matrix
c     d is a scalar
c     c is the matrix product of diag(exp(d*a))*b
      implicit none
      integer i, j, m
      double precision a(m), b(m, m), c(m, m), d
      i = 1
      do while(i .le. m)
          j = 1
          do while(j .le. m)
              c(i, j) = dexp(d*a(i))*b(i, j)
              j = j+1
          enddo
          i = i+1
      enddo
      end


      subroutine multi5(m, k, a, b, c)
c     a is (m*m) matrix
c     b is (m*m) matrix
c     c is (m*m) matrix
c     c is the matrix product of a*diag(b[,k])
c     where b[,k] is the kth column of b
      implicit none
      integer i, j, m, k
      double precision a(m,m), b(m,m), c(m,m)
      i = 1
      do while(i .le. m)
          j = 1
          do while(j .le. m)
              c(i,j) = a(i,j)*b(j,k)
              j = j+1
          enddo
          i = i+1
      enddo
      end


      subroutine multi6(m, k, a, b, c)
c     a is (m*m) matrix
c     b is (m*m) matrix
c     c is (m*m) matrix
c     c is the matrix product of diag(a[,k])*b
c     where a[k,] is the kth row of a
      implicit none
      integer i, j, m, k
      double precision a(m,m), b(m,m), c(m,m)
      i = 1
      do while(i .le. m)
          j = 1
          do while(j .le. m)
              c(i,j) = a(k,i)*b(i,j)
              j = j+1
          enddo
          i = i+1
      enddo
      end
