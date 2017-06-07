      PROGRAM Main
      IMPLICIT none
      COMMON/DATA/ a, b, step
      COMMON/Grid/ grid
      INTEGER step
      REAL*8 grid(17000), y, Square, Simpson, Gaus_5, a, b
      OPEN (2, FILE = 'square.txt')
      OPEN (3, FILE = 'simpson.txt')
      OPEN (4, FILE = 'gaus5.txt')
      step = 1
      a = -1.d0
      b = 5.d0
7     IF(step.LT.17000)THEN
C          WRITE(2, *) step
C          WRITE(3, *) step
C          WRITE(4, *) step
         CALL MakeGrid(grid)
         y = Square(grid)
         WRITE(2, *)y
         y = Simpson(grid)
         WRITE(3, *)y
         y = Gaus_5(grid)
         WRITE(4, *)y
         step = step*2
         PRINT *, step
         GOTO 7
      ENDIF
      CLOSE(2)
      CLOSE(3)
      CLOSE(4)
      PAUSE
      END
       

      REAL*8 FUNCTION fun(x)
      REAL*8 x
      fun  = 6*x**13+4*x*x+7
      END 


      SUBROUTINE MakeGrid(grid)
      REAL*8 grid(*), part, x, a, b
      INTEGER step, i
      COMMON/DATA/ a, b, step
      part = (b - a) / step
      i = 1
      grid(i) = a
      do i = 1, step + 1 
         x = a + (i - 1)*part
         grid(i) = x
      enddo
      END


      REAL*8 FUNCTION Simpson(grid)
      IMPLICIT NONE
      REAL*8 grid(*), y, part, middle, a, b, fun
      INTEGER step, i
      COMMON/DATA/ a, b, step

      part = ((b - a) / step)
      y = 0.d0
      i = 1
      
      do i = 1, step
         middle = grid(i + 1) + grid(i)
         middle = middle / 2
         y = y + part / 6.d0 * (fun(grid(i)) + 4.d0 * fun(middle) +
     +   fun(grid(i + 1)))
      enddo
      Simpson = y
      END
      
      REAL*8 FUNCTION Square(grid)
      IMPLICIT NONE
      COMMON/DATA/ a, b, step
      REAL*8 grid(*), y, part, middle, a, b, fun
      INTEGER step, i

      part = ((b - a) / step)
      y = 0.d0
      do i = 1, step
         middle = grid(i + 1) + grid(i)
         middle = middle / 2.d0
         y = y + part * fun(middle)
      enddo
      Square = y
      END


      real*8 function Gaus_5(grid)
      implicit none
      COMMON/DATA/ a, b, step
      real*8 grid(*), fun, q(5), x(5), sum1, sum, a, b
      integer j, k, step

      x(1) = -0.906179846
      x(2) = -0.538469310
      x(3) = 0.0
      x(4) = 0.538469310
      x(5) = 0.906179846
      q(1) = 0.236926885
      q(2) = 0.478628670
      q(3) = 0.568888889
      q(4) = 0.478628670
      q(5) = 0.236926885 
      sum = 0.0
      do j = 1, 5
         sum1 = 0.0
         do k=1 , step
            sum1 = sum1 + (grid(k + 1) - grid(k)) * fun((grid(k) + 
     +      grid(k + 1) + x(j) * (grid(k + 1) - grid(k))) / 2.0 )
         enddo
         sum = sum + q(j) * sum1
      enddo
      Gaus_5 = sum / 2.0
      return
      end