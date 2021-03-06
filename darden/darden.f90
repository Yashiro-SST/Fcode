!=======================================================================!
! program  Culculation of  Darden's  F function	parameter			  				!
!	input	 :	input_darden.txt									                        	!
!	output	:	F function parameter.txt								            				!
! Variable List
! A -------------Darden's variable
! a_inf ---------speed of sound @ flight hight
! B -------------Darden's variable
! beta-----------beta = sqrt(abs(M_inf**2.0d0 - 1.0d0))
! C -------------Darden's variable
! C0 ------------initial value of C
! D -------------Darden's variable
! D0 ------------initial value of D
! dn ------------Division number @ integral
! eps -----------allowable error
! errC ----------
! errD ----------
! F -------------Darden's F function (array) 
! F10 -----------
! F1C -----------Partial derivative of F1 with respect to C
! F1D -----------Partial derivative of F1 with respect to D
! F20 -----------
! F2C -----------Partial derivative of F2 with respect to C
! F2D -----------Partial derivative of F2 with respect to D
! FYR_C1 --------
! FYR_C2 --------
! FYR_C3 --------
! FYR_C4 --------
! FYR_int1 ------
! FYR_int2 ------
! FYR_int3 ------
! FYR_int4 ------
! g -------------Gravitational acceleration
! gamma ---------Specific heat ratio
! h_inf ---------Flight hight 
! i -------------loop number
! ite_max--------Maximum number of iterations
! j -------------loop number
! k -------------Variable when calculating S 
! kr ------------Reflection coefficient of the ground
! l -------------Effective length
! M_inf----------Mach number
! mu ------------front/rear shock ratio
! P
! pi ------------Pi 
! Q1 ------------1st interval of integration in F2
! Q2 ------------2nd interval of integration in F2
! Q3 ------------3rd interval of integration in F2
! Q4 ------------4th interval of integration in F2
! QC1 -----------Partial derivative of Q1 with respect to C
! QC2 -----------Partial derivative of Q2 with respect to C
! QC3 -----------Partial derivative of Q3 with respect to C
! QC4 -----------Partial derivative of Q4 with respect to C
! r -------------Gas constant(Air)
! s 
! T_inf ---------Tenparature @ Flight height 
! u_inf ---------Flight speed 
! unit-----------Unit System Change yardpond(0), MKS(1)
! w -------------Weight 
! x 
! y
! yf ------------Darden's variable
! z 
!=======================================================================!

!---module function of solve integral F1, F2---!

module Trapezoidal_integral
  implicit none
  contains

    function FYR_integral1(l, yf, yr, A, dn) Result(FYR_int1)
      
      real(8), intent(in) :: l, yf, yr, A
      integer, intent(in) :: dn
      real(8) dx, x, y, sum, FYR_int1, K
      integer i, n

      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((yf / 2.0d0) / dx)

      do i = 0, n
        x = dx * dble(i)
        K = sqrt(l - x) / (yr - x)
        y = 2.0d0 * x / yf * A * K

        if (i == 0 .or. i ==n) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      FYR_int1 = sum * dx

    end function FYR_integral1

    function FYR_integral2(l, yf, yr, A, C0, dn) Result(FYR_int2)
      
      real(8), intent(in) :: l, yf, yr, A, C0
      integer, intent(in) :: dn
      real(8) dx, x, y, sum, FYR_int2, K
      integer i, n

      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((yf / 2.0d0) / dx)

      do i = 0, n
        x = yf / 2.0d0 + dx * dble(i)
        K = sqrt(l - x) / (yr - x)
        y = K * (C0 * (2.0d0 * x / yf - 1.0d0) - A * (2.0d0 * x / yf -2.0d0))

        if (i == 0 .or. i ==n) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      FYR_int2 = sum * dx

    end function FYR_integral2

    function FYR_integral3(l, yf, yr, lam, B, C0, dn) Result(FYR_int3)
      
      real(8), intent(in) :: l, yf, yr, lam, B, C0
      integer, intent(in) :: dn
      real(8) dx, x, y, sum, FYR_int3, K
      integer i, n

      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((lam - yf) / dx)

      do i = 0, n
        x = yf + dx * dble(i)
        K = sqrt(l - x) / (yr - x)
        y = K * (B * (x - yf) + C0)

        if (i == 0 .or. i ==n) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      FYR_int3 = sum * dx

    end function FYR_integral3

    function FYR_integral4(l, yf, yr, lam, B, D0, dn) Result(FYR_int4)
      
      real(8), intent(in) :: l, yf, yr, lam, B, D0
      integer, intent(in) :: dn
      real(8) dx, x, y, sum, FYR_int4, K ,lx
      integer i, n

      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((l - lam) / dx)

      do i = 0, n
        x = lam + dx * dble(i)
        lx = l - x

        if(lx <= 0.0d0) then
          K = 0.0d0
        else
          K = sqrt(l - x) / (yr - x)
        end if

        y = K * (B * (x - yf) - D0)

        if (i == 0 .or. i ==n) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      FYR_int4 = sum * dx

    end function FYR_integral4

    function F1initial(l, S, yf, yr, B, D0, FYR_int1, FYR_int2, FYR_int3, FYR_int4) Result(F10)

      real(8), intent(in) :: l, S, yf, yr, B, D0, FYR_int1, FYR_int2, FYR_int3, FYR_int4
      real(8) FYR, pi, F10

      pi = 2.0d0 * acos(0.0d0)
      
      FYR = -1.0d0 * (FYR_int1 + FYR_int2 + FYR_int3 + FYR_int4) / (pi * sqrt(abs(yr - l)))
      F10 = S * (yr - l) + B * (l - yf) - D0 - FYR
    
    end function F1initial

    function FYR_partialC1(l, yf, yr, A, C0, S, mu, dn) Result(FYR_C1)

      real(8), intent(in) :: l, yf, yr, A, C0, S, mu
      integer, intent(in) :: dn
      real(8) dx, x, y, sum, FYR_C1, K, KC
      integer i, n

      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((yf / 2.0d0) / dx)

      do i = 0, n
        x = dx * dble(i)
        K = sqrt(l - x) / (yr - x)
        KC = -S * mu * sqrt(l - x) / ((C0 + S * mu * (l - x))**2.0d0)
        y = KC * 2.0d0 * x * A / yf + K * 2.0d0 * x / yf * (2.0d0 * C0 / (S * yf) - 0.5d0)

        if (i == 0 .or. i ==n) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      FYR_C1 = sum * dx

    end function FYR_partialC1

    function FYR_partialC2(l, yf, yr, A, C0, S, mu, dn) Result(FYR_C2)

      real(8), intent(in) :: l, yf, yr, A, C0, S, mu
      integer, intent(in) :: dn
      real(8) dx, x, y, sum, FYR_C2, K, KC
      integer i, n

      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((yf / 2.0d0) / dx)

      do i = 0, n
        x = (yf / 2.0d0) + dx * dble(i)
        K = sqrt(l - x) / (yr - x)
        KC = -S * mu * sqrt(l - x) / ((C0 + S * mu * (l - x))**2.0d0)
        y = KC * (C0 * (2.0d0 * x / yf - 1.0d0) - A * (2.0d0 * x / yf - 2.0d0)) &
            + K * (2.0d0 * x / yf - 1.0d0) - (2.0d0 * C0 / (S * yf) - 0.50d0) &
            * (2.0d0 * x / yf - 2.0d0)

        if (i == 0 .or. i ==n) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      FYR_C2 = sum * dx

    end function FYR_partialC2

    function FYR_partialC3(l, S, yf, yr, lam, B, C0, mu, dn) Result(FYR_C3)

      real(8), intent(in) :: l, S, yf, yr, lam, B, C0, mu
      integer, intent(in) :: dn
      real(8) dx, x, y, sum, FYR_C3, K, KC
      integer i, n

      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((lam - yf) / dx)

      do i = 0, n
        x = yf + dx * dble(i)
        K = sqrt(l - x) / (yr - x)
        KC = -S * mu * sqrt(l - x) / ((C0 + S * mu * (l - x))**2.0d0)
        y = KC * (B * (x - yf) + C0) + K

        if (i == 0 .or. i ==n) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      FYR_C3 = sum * dx

    end function FYR_partialC3

    function FYR_partialC4(l, S, yf, lam, B, C0, D0, mu, dn) Result(FYR_C4)
      
      real(8), intent(in) :: l, S, yf, lam, B, mu, C0, D0
      integer, intent(in) :: dn
      real(8) dx, x, y, sum, FYR_C4, KC, lx
      integer i, n

      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((l - lam) / dx)

      do i = 0, n
        x = lam + dx * dble(i)
        lx = l - x

        if(lx <= 0.0d0) then
          KC = 0.0d0
        else
          KC = -S * mu * sqrt(l - x) / ((C0 + S * mu * (l - x))**2.0d0)
        end if

        y = KC * (B * (x - yf) - D0)

        if (i == 0 .or. i ==n) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      FYR_C4 = sum * dx

    end function FYR_partialC4

    function F1_partialC(l, S, yr, mu, C0, FYR_int1, FYR_int2, FYR_int3, FYR_int4, &
      FYR_C1, FYR_C2, FYR_C3, FYR_C4) Result(F1C)
      
      real(8), intent(in) :: l, S, yr, mu, C0, FYR_int1, FYR_int2, FYR_int3, FYR_int4, &
      FYR_C1, FYR_C2, FYR_C3, FYR_C4
      real(8) pi, FYRC, F1C

      pi = 2.0d0 * acos(0.0d0)
      FYRC = sqrt(S * mu) * (C0**(-1.50d0)) / 2.0d0 * pi &
            * (FYR_int1 + FYR_int2 + FYR_int3 + FYR_int4) &
            - (FYR_C1 + FYR_C2 + FYR_C3 + FYR_C4) / (pi * sqrt(yr - l))
      F1C = 1.0d0 / mu - FYRC 

    end function F1_partialC

    function F1_partialD(l, yr, lam, dn) Result(F1D)

      real(8), intent(in) :: l, yr, lam
      integer, intent(in) :: dn
      real(8) dx, x, y, sum, F1D, lx, K, pi, FYRD
      integer i, n

      pi = 2.0d0 * acos(0.0d0)
      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((l - lam) / dx)

      do i = 0, n
        x = lam + dx * dble(i)
        lx = l - x

        if(lx <= 0.0d0) then
          K = 0.0d0
        else
          K = sqrt(l - x) / (yr - x)
        end if
        y = K

        if (i == 0 .or. i ==n) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      FYRD = - sum * dx / (pi * sqrt(yr - l))
      F1D = - 1.0d0 - FYRD * dx

    end function F1_partialD

    function Qterm1(l, yf, yr, A, dn) result(Q1)
      
      real(8), intent(in) :: l, yf, yr, A
      integer, intent(in) :: dn
      real(8) dx, x, y, t, sum, pi, Q1
      integer i, n

      pi = 2.0d0 * acos(0.0d0)
      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((yf / 2.0d0) / dx)

      do i = 0, n
        x = dx * dble(i)
        t = atan(sqrt((yr - l) / (l - x)))
        y = 2.0d0 * x / yf * A * t

        if (i == 0 .or. i ==n) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      Q1 = - 2.0d0 / pi * sum * dx

    end function Qterm1

    function Qterm2(l, yf, yr, A, C0, dn) result(Q2)
      
      real(8), intent(in) :: l, yf, yr, A, C0
      integer, intent(in) :: dn
      real(8) dx, x, y, t, sum, pi, Q2
      integer i, n

      pi = 2.0d0 * acos(0.0d0)
      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((yf / 2.0d0) / dx)

      do i = 0, n
        x = (yf / 2.0d0) + dx * dble(i)
        t = atan(sqrt((yr - l) / (l - x)))
        y = (C0 * (2.0d0 * x / yf - 1.0d0) - A * (2.0d0 * x / yf -2.0d0) ) * t

        if (i == 0 .or. i ==n) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      Q2 = - 2.0d0 / pi * sum * dx

    end function Qterm2

    function Qterm3(l, yf, yr, lam, B, C0, dn) Result(Q3)
      
      real(8), intent(in) :: l, yf, yr, lam, B, C0
      integer, intent(in) :: dn
      real(8) dx, x, y, t, sum, pi, Q3
      integer i, n

      pi = 2.0d0 * acos(0.0d0)
      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((lam - yf) / dx)

      do i = 0, n
        x = yf + dx * dble(i)
        t = atan(sqrt((yr - l) / (l - x)))
        y = (B * (x - yf) + C0) * t

        if (i == 0 .or. i ==n) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      Q3 = - 2.0d0 / pi * sum * dx

    end function Qterm3

    function Qterm4(l, yf, yr, lam, B, D0, dn) Result(Q4)
      
      real(8), intent(in) :: l, yf, yr, lam, B, D0
      integer, intent(in) :: dn
      real(8) dx, x, y, t, sum, pi, lx, Q4
      integer i, n

      pi = 2.0d0 * acos(0.0d0)
      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((l - lam) / dx)

      do i = 0, n
        x = lam + dx * dble(i)
        lx = l - x
        if(lx <= 0.0d0) then
          t = pi / 2.0d0
        else
          t = atan(sqrt((yr - l) / (l - x)))
        end if
        y = (B * (x - yf) - D0) * t

        if (i == 0 .or. i == n) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      Q4 = - 2.0d0 / pi * sum * dx

    end function Qterm4

    function F2initial(l, S, yf, yr, B, D0, Q1, Q2, Q3, Q4) Result(F20)

      real(8), intent(in) :: l, S, yf, yr, B, D0, Q1, Q2, Q3, Q4
      real(8) Q_sum, F20
      
      Q_sum = Q1 + Q2 + Q3 + Q4
      F20 = (B * (l - yf) - D0 + S / 2.0d0 *(yr - l)) * (yr - l) - Q_sum

    end function F2initial

    function Q_partialC1(l, S, yf, yr, A, C0, mu, dn) result(QC1)
      real(8), intent(in) :: l, S, yf, yr, A, C0, mu
      integer, intent(in) :: dn
      real(8) dx, x, y, t, tc, sum, pi, QC1
      integer i, n

      pi = 2.0d0 * acos(0.0d0)
      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((yf / 2.0d0) / dx)

      do i = 0, n
        x = dx * dble(i)
        t = atan(sqrt((yr - l) / (l - x)))
        tc = sqrt(S * mu * (l - x) / C0) / (2.0d0 * (C0 - S * mu * (l - x)))
        y = 2.0d0 * x / yf * ((2.0d0 * C0 / (S * yf) - 0.5d0) * t + A * tc)

        if (i == 0 .or. i ==n) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      QC1 = - 2.0d0 / pi * sum * dx

    end function Q_partialC1

    function Q_partialC2(l, S, yf, yr, A, C0, mu, dn) result(QC2)
      real(8), intent(in) :: l, S, yf, yr, A, C0, mu
      integer, intent(in) :: dn
      real(8) dx, x, y, t, tc, sum, pi, QC2
      integer i, n

      pi = 2.0d0 * acos(0.0d0)
      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((yf / 2.0d0) / dx)

      do i = 0, n
        x = dble(nint(yf / 2.0d0)) + dx * dble(i)
        t = atan(sqrt((yr - l) / (l - x)))
        tc = sqrt(S * mu * (l - x) / C0) / (2.0d0 * (C0 - S * mu * (l - x)))
        y = ((2.0d0 * x / yf - 1.0d0) - (2.0d0 * C0 / (S * yf) - 0.5d0) &
            * (2.0d0 * x / yf - 2.0d0)) * t &
            + (C0 * (2.0d0 * x / yf - 1.0d0) - A * (2.0d0 * x / yf - 2.0d0)) * tc

        if (i == 0 .or. i ==n) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      QC2 = - 2.0d0 / pi * sum * dx

    end function Q_partialC2

    function Q_partialC3(l, S, yf, yr, lam, B, C0, mu, dn) result(QC3)
      real(8), intent(in) :: l, S, yf, yr, lam, B, C0, mu
      integer, intent(in) :: dn
      real(8) dx, x, y, t, tc, sum, pi, QC3
      integer i, n

      pi = 2.0d0 * acos(0.0d0)
      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((lam - yf) / dx)

      do i = 0, n
        x = yf + dx * dble(i)
        t = atan(sqrt((yr - l) / (l - x)))
        tc = sqrt(S * mu * (l - x) / C0) / (2.0d0 * (C0 - S * mu * (l - x)))
        y = t + (B * (x - yf) + C0) * tc
        
        if (i == 0 .or. i == n) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      QC3 = - 2.0d0 / pi * sum * dx

    end function Q_partialC3

    function Q_partialC4(l, S, yf, lam, B, C0, D0, mu, dn) result(QC4)
      real(8), intent(in) :: l, S, yf, lam, B, C0, D0, mu
      integer, intent(in) :: dn
      real(8) dx, x, y, tc, sum1, sum2, pi, sp, QC4
      integer i, n

      pi = 2.0d0 * acos(0.0d0)
      sp = l - C0 / (S * mu)
      dx =  l / dble(dn)

      if (sp > lam) then
        sum1 = 0.0d0
        sum2 = 0.0d0 
        n = nint((sp - lam) / dx)

        do i = 0, n-1
          x = lam + dx * dble(i)
          tc = sqrt(S * mu * (l - x) / C0) / (2.0d0 * (C0 - S * mu * (l - x)))
          y = (B * (x - yf) - D0) * tc

          if (i == 0 .or. i == n-1) then
          sum1 = sum1 + 0.5d0 * y
          else
          sum1 = sum1 + y
          end if
        enddo

        n = nint((l - sp) / dx)

        do i = 1, n-1
          x = sp + dx * dble(i)
          tc = sqrt(S * mu * (l - x) / C0) / (2.0d0 * (C0 - S * mu * (l - x)))
          y = (B * (x - yf) - D0) * tc

          if (i == 1 .or. i == n-1) then
          sum2 = sum2 + 0.5d0 * y
          else
          sum2 = sum2 + y
          end if
        enddo

        QC4 = - 2.0d0 / pi * (sum1 + sum2) * dx
      
      else if(sp <= lam) then
        sum1 = 0.0d0
        n = nint((l - lam) / dx)

        do i = 0, n-1
          x = lam + dx * dble(i)
          tc = sqrt(S * mu * (l - x) / C0) / (2.0d0 * (C0 - S * mu * (l - x)))
          y = (B * (x - yf) - D0) * tc

          if (i == 0 .or. i == n-1) then
          sum1 = sum1 + 0.5d0 * y
          else
          sum1 = sum1 + y
          end if
        enddo

        QC4 = - 2.0d0 / pi * sum1 * dx

      end if

    end function Q_partialC4

    function F2_partialC(l, S, yf, B, mu, C0, D0, QC1, QC2, QC3, QC4) result(F2C)
      real(8), intent(in) :: l, S, yf, B, mu, C0, D0, QC1, QC2, QC3, QC4
      real(8) QC_sum, smu, F2C

      smu = S * mu
      QC_sum = QC1 + QC2 + QC3 + QC4
      F2C = (B * (l - yf) + C0 / mu - D0) / (S * mu) - QC_sum

    end function F2_partialC

    function F2_partialD(l, yr, lam, dn) Result(F2D)

      real(8), intent(in) :: l, yr, lam
      integer, intent(in) :: dn
      real(8) dx, x, y, sum, pi, QD, F2D
      integer i, n

      pi = 2.0d0 * acos(0.0d0)
      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((l - lam) / dx)

      do i = 0, n-1
        x = lam + dx * dble(i)
        y = atan(sqrt((yr - l) / (l - x))) 

        if (i == 0 .or. i == n-1) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      QD = 2.0d0 / pi * sum * dx
      F2D = -yr + l - QD

    end function F2_partialD

end module Trapezoidal_integral

!---module function of calculate various variables---!

module Darden_variable
  implicit none
  contains

    function cal_beta(M_inf) Result(beta)
      real(8), intent(in) :: M_inf
      real(8) beta

      beta = sqrt(abs(M_inf**2.0d0 - 1.0d0))

    end function cal_beta

    function cal_ainf(gamma, R, T_inf, unit) Result(a_inf)
      real(8), intent(in) :: gamma, R, T_inf
      real(8) a_inf
      integer, intent(in) :: unit
      
      if (unit == 1) then
        a_inf = sqrt(abs(gamma * R * T_inf))
      else if (unit == 0) then
        a_inf = sqrt(abs(gamma * R * T_inf)) * 3.28084
      else
        stop 'cal_a error, invalid unit system number is entered!!!'
      end if

    end function cal_ainf

    function cal_u(M_inf, a_inf) Result(u_inf)
      real(8), intent(in) :: M_inf, a_inf
      real(8) u_inf

      u_inf = M_inf * a_inf

    end function cal_u

    function cal_k(gamma, M_inf, beta) Result(k)
      real(8), intent(in) :: gamma, M_inf, beta
      real(8) k

      k = (gamma + 1) * M_inf**4.0d0 / sqrt(2.0d0 * beta**3.0d0)

    end function cal_k

    function cal_S(k, h_inf, l, unit) Result(S)
      real(8), intent(in) :: k, h_inf, l
      real(8) S
      integer, intent(in) :: unit

      if (unit == 1) then
        S = 1.0d0 / (k * sqrt(abs(h_inf / l)))
      else if (unit == 0) then
        S = 1.0d0 / (k * sqrt(abs(h_inf / l)))
        !S = 0.00026967727
      else
        stop 'cal_S error, invalid unit system number is entered!!!'
      end if

    end function cal_S

    !function cal_S(k, h_inf, l) Result(S)
      !real(8), intent(in) :: k, h_inf, l
     ! real(8) S

     !S = 1.0d0 / (k * sqrt(abs(h_inf / l)))

    !end function cal_S

    function cal_S2(gamma, M_inf, beta, u_inf, h_inf, R, dn, unit) Result(S2)
      real(8), intent(in) :: gamma, M_inf, beta, u_inf, h_inf, R
      real(8) gamma_l, sum, Tc, zc, T_slope, Mh, Th, g, dz, z_sp, x, y, z
      real(8) T_z, Ma_z, beta_z, rt_ratio, p_ratio, rho_ratio, T_ratio
      real(8) ray_int, S2_int, S2
      integer, intent(in) :: dn, unit
      integer i, j

      Tc = 216.65        !---temparature of Stratosphere
      if (unit == 1) then
        T_slope = 0.0065          !---difference of temparature[K] per a hight[m]
      else if (unit == 0) then
        T_slope = 0.0065 / 3.28084   !--- 1m = 3.28084ft
      else
        stop 'cal_T_slope error, invalid unit system number is entered!!!'
      end if

      if (unit == 1) then
        zc = 11000          !---hight of top of Troposphere / bottom of Stratosphere
      else if (unit == 0) then
        zc = 11000 * 3.28084
      else
        stop 'cal_zc error, invalid unit system number is entered!!!'
      end if
      
      Mh = M_inf         !---if hight where calc initial wave form is higher than zc 
      Th = Tc            !---same as above
      gamma_l = (gamma + 1) / 2

      !---integral in the ray tube ratio formura---!

      dz = h_inf / dn
      z_sp = h_inf - zc
      sum = 0.0d0

      do j = 0, dn
        !---z is vertical distance of signal from airplane axis---!
        z = dz * dble(j)
        if (0 <= z .and. z <= z_sp) then               !---temparature value at Stratosphere
          T_z = Tc
        else if (z_sp < z .and. z <= h_inf) then       !---temparature value at Troposphere
          T_z = Tc + T_slope * z
        else 
          stop 'cal_S2 error, variable z is out of range !!!'
        end if
        
        Ma_z = u_inf / sqrt(abs(gamma * R * T_z))   !---cal. Mach number @ z in this loop
        y = 1.0d0 / ((Ma_z**2.0d0 - 1.0d0)**0.50d0) !---cal. integration in ray tube formura

        if (i == 0 .or. i ==dn) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      ray_int = sum * dz

      !---integral in the S2 formura---!

      sum = 0.0d0   !---zero clear---!
      z = 0.0d0
      y = 0.0d0
      !dz = h_inf / dn
      !z_sp = h_inf - zc

      do i = 0, dn
        z = dz * dble(i)
        if (0 <= z .and. z <= z_sp) then             !---temparature value at Stratosphere
          T_z = Tc
        else if (z_sp < z .and. z <= h_inf) then     !---temparature value at Troposphere
          T_z = Tc + T_slope * z
        else 
          stop 'cal_S2 error, variable z is out of range !!!'
        end if
        Ma_z = u_inf / sqrt(abs(gamma * R * T_z))    !---cal. Mach number @ z in this loop
        beta_z = sqrt(abs(Ma_z**2.0d0 - 1.0d0))      !---cal. beta @z in this loop

        rt_ratio = (M_inf * (1.0d0 - 1.0d0 / Ma_z**2.0d0)**0.5d0 * ray_int)**(-1.0d0)                         !---
        p_ratio = exp(g * (z - zc) / (R * Tc))       
        rho_ratio = exp(- g * (z - zc) / (R * Tc))   
        T_ratio = sqrt(Th / T_z)

        x = p_ratio * (rho_ratio * T_ratio)**0.5d0 * (rt_ratio)**0.5d0 * Ma_z / beta_z

        if (i == 0 .or. i ==dn) then
          sum = sum + 0.5d0 * x
        else
          sum = sum + x
        end if
      enddo

      S2_int = sum * dz

      !---calculate Slope of balancing line (limit of parameter B)---!
      S2 = sqrt(2.0d0 * beta) / (gamma_l * Mh**3.0d0 * S2_int)

    end function cal_S2

    function cal_Se_l(beta, W, rho_inf, u_inf, unit) Result(Se_l)
      real(8), intent(in) :: beta, W, rho_inf, u_inf
      real(8) g, Se_l
      integer, intent(in) :: unit

      if (unit == 1) then
        g = 9.80665
      else if (unit == 0) then
        g = 9.80665 * 3.28084
      else
        stop 'cal_Ae error, invalid unit system number is entered!!!'
      end if

      Se_l = beta * W * g /(rho_inf * u_inf**2.0d0)

    end function cal_Se_l

    function set_B(Sn, percen, sig) Result(B)
      real(8), intent(in) ::Sn, percen
      real(8) B
      integer, intent(in) :: sig
      
      if (sig == 0) then
        B = 0.0d0
      else if (sig == 1) then
        B = percen * Sn
      else
        stop 'set_B error! please check sig is integer. percen less than 1. '
      end if
    
    end function set_B


    function cal_A(C0, S, yf) Result(A)
      real(8), intent(in) :: C0, S, yf
      real(8) A

      A = ( C0**2.0d0 / ( S * yf ) ) - C0 / 2.0d0

    end function cal_A

    function cal_yr(l, C0, S, mu) Result(yr)
      real(8), intent(in) :: l, C0, S, mu
      real(8) yr

      yr = l + C0 / ( S * mu )

    end function cal_yr

    function cal_lyf(l, yf) Result(lyf)
      real(8), intent(in) :: l, yf
      real(8) lyf

      lyf = l - yf

    end function cal_lyf

    function cal_lyfh(l, yf) Result(lyfh)
      real(8), intent(in) :: l, yf
      real(8) lyfh

      lyfh = l - (yf / 2.0d0)

    end function cal_lyfh

    function cal_lam1(A, yf, l) Result(lam1)
      real(8), intent(in) :: A, yf, l
      real(8) lam1

      lam1 = 32d0 * A * (l**2.5) / (15d0 * yf)

    end function cal_lam1

    function cal_lam2(A, C0, yf, lyfh) Result(lam2)
      real(8), intent(in) :: A, C0, yf, lyfh
      real(8) lam2

      lam2 = 32d0 * (C0 - 2.0d0 * A) * (lyfh**2.5) / (15d0 * yf)

    end function cal_lam2

    function cal_lam3(A, B, C0, yf, lyf) Result(lam3)
      real(8), intent(in) :: A, B, C0, yf, lyf
      real(8) lam3

      lam3 = 16d0 * (2.0d0 * A + B * yf - 2.0d0 * C0) * (lyf**2.5) / (15d0 * yf)

    end function cal_lam3

    function cal_lam4(lam1, lam2, lam3, Se_l) Result(lam4)
      real(8), intent(in) :: lam1, lam2, lam3, Se_l
      real(8) lam4

      lam4 = lam1 + lam2 + lam3 - Se_l

    end function cal_lam4

    function cal_lam(l, lam4,C0, D0) Result(lam)
      real(8), intent(in) :: l, lam4,C0, D0
      real(8) lam

      lam = l - (3.0d0 * (lam4) / (8.0d0 * (C0 + D0))) ** (2.0d0 / 3.0d0)

    end function cal_lam
  
end module Darden_variable

!---main program of calculate Ffunc parameter---!

program darden
  use Trapezoidal_integral
  use Darden_variable
  implicit none

  !---Variable definition---!

  real(8) A, C, D, lam, yr
  real(8) B, mu, yf, W, gamma, M_inf, h_inf, l, R, T_inf, rho_inf
  real(8) S, S2, Sn, k, beta, u_inf, a_inf, Se_l, percen
  real(8) C0, dC, D0, dD
  real(8) C0_n, D0_n, CaddD, A_n, lam1_n, lam2_n, lam3_n, lam4_n, lam_n
  real(8) eps, errC, errD
  real(8) lyf, lyfh, lam1, lam2, lam3, lam4
  real(8) F10, F20, F1C, F2C, F1D, F2D
  real(8) FYR_int1, FYR_int2, FYR_int3, FYR_int4
  real(8) FYR_C1, FYR_C2, FYR_C3, FYR_C4
  real(8) Q1, Q2, Q3, Q4
  real(8) QC1, QC2, QC3, QC4
  integer dn, i, j, ite_max, sig, unit, Scal

  !---Reading paramator---!
  open(10,file='input_darden.txt')

  read(10,*) mu                  !---front/rear shock ratio
  read(10,*) sig                 !---0:minimum-overpressure signature, 1:minimum-shock
  read(10,*) percen              !---if M.S.sig, B = percen * S 
  read(10,*) Scal                !---0:S calc on Isothermal atmosphere, 1: Real atmosphere
  read(10,*) M_inf               !---Mach number
  read(10,*) h_inf               !---Flight height
  read(10,*) W                   !---Weight [kg]
  read(10,*) l                   !---Effective length
  read(10,*) yf                  !---yf position
  read(10,*) rho_inf             !---density @Flight height
  read(10,*) T_inf               !---Tenparature @Flight height
  read(10,*) gamma               !---Specific heat ratio
  read(10,*) R                   !---Gas constant(Air)
  read(10,*) C0                  !---initial C
  read(10,*) D0                  !---initial D
  read(10,*) dn                  !---Division number @ integral
  read(10,*) ite_max             !---Maximum number of iterations
  read(10,*) unit                !---Unit System Change yardpond(0), MKS(1)
  read(10,*) eps                 !---allowable error "epsilon"

  close(10)

  !---Output input data---!

  write(*,*) '-----------------------input data-------------------------'

  write(*,*) 'mu = ', mu
  write(*,*) 'sig = ', sig
  write(*,*) 'percen = ', percen
  write(*,*) 'Scal = ', Scal
  write(*,*) 'M_inf = ', M_inf
  write(*,*) 'h_inf = ', h_inf
  write(*,*) 'W = ', W
  write(*,*) 'l = ', l
  write(*,*) 'yf = ', yf
  write(*,*) 'rho_inf = ', rho_inf
  write(*,*) 'T_inf = ', T_inf
  write(*,*) 'gamma = ', gamma
  write(*,*) 'R = ', R
  write(*,*) 'C0 = ', C0
  write(*,*) 'D0 = ', D0
  write(*,*) 'dn = ', dn
  write(*,*) 'ite_max = ', ite_max
  write(*,*) 'unit = ', unit
  write(*,*) 'eps = ', eps

  !---Calculate Known function---!

  write(*,*) '---------------------Known function---------------------------'

  beta  = cal_beta(M_inf)
  a_inf = cal_ainf(gamma, R, T_inf, unit)
  u_inf = cal_u(M_inf, a_inf)
  k     = cal_k(gamma, M_inf, beta)
  S     = cal_S(k, h_inf, l, unit)
  S2    = cal_S2(gamma, M_inf, beta, u_inf, h_inf, R, dn, unit)
  Se_l  = cal_Se_l(beta, W, rho_inf, u_inf, unit)
  
  if(Scal == 0) then
    Sn = S
  else if(Scal ==1) then
    Sn = S2
  end if

  B = set_B(Sn, percen, sig)

  !---Output known function---!

  write(*,*) 'beta = ', beta
  write(*,*) 'a_inf = ', a_inf
  write(*,*) 'u_inf = ', u_inf
  write(*,*) 'k = ', k
  write(*,*) 'S = ', S
  write(*,*) 'S2 = ', S2
  write(*,*) 'Se_l = ', Se_l
  write(*,*) 'B = ' , B

  !---Newton method---!

  write(*,*) '-------------------Newton loop started-----------------------------'

  open(30, file='result.txt')
  write(30,*) 'i,   dC,   dD,   C0,   D0,   lam4,   lam'

  do i = 0, ite_max

    write(*,*) '__________________________________________________________'

    write(*,*) 'iteration' , i+1
    write(*,*) 'C0 = ', C0
    write(*,*) 'D0 = ', D0

    !---Calculate unknown funcion---!
    A    = cal_A(C0, Sn, yf)
    yr   = cal_yr(l, C0, Sn, mu)
    lyf  = cal_lyf(l, yf)
    lyfh = cal_lyfh(l, yf)
    lam1 = cal_lam1(A, yf, l)
    lam2 = cal_lam2(A, C0, yf, lyfh)
    lam3 = cal_lam3(A, B, C0, yf, lyf)
    lam4 = cal_lam4(lam1, lam2, lam3, Se_l)
    lam  = cal_lam(l, lam4,C0, D0)

    write(*,*) 'A = ', A
    write(*,*) 'yr = ', yr
    !write(*,*) 'lam1 = ', lam1
    !write(*,*) 'lam2 = ', lam2
    !write(*,*) 'lam3 = ', lam3
    write(*,*) 'lam4 = ', lam4
    write(*,*) 'lam = ', lam

    if(lam4 <= 0) then
      write(*,*) 'lam = NaN!!!'
      stop
    end if

    !---Calculation F10, F20, F1 partial CorD, and F2 partial CorD---!

    FYR_int1 = FYR_integral1(l, yf, yr, A, dn)
    !write(*,*) 'FYR_int1 = ', FYR_int1
    FYR_int2 = FYR_integral2(l, yf, yr, A, C0, dn)
    !write(*,*) 'FYR_int2 = ', FYR_int2
    FYR_int3 = FYR_integral3(l, yf, yr, lam, B, C0, dn)
    !write(*,*) 'FYR_int3 = ', FYR_int3
    FYR_int4 = FYR_integral4(l, yf, yr, lam, B, D0, dn)
    !write(*,*) 'FYR_int4 = ', FYR_int4
    F10 = F1initial(l, Sn, yf, yr, B, D0, FYR_int1, FYR_int2, FYR_int3, FYR_int4)
    write(*,*) 'F10 = ', F10

    FYR_C1 = FYR_partialC1(l, yf, yr, A, C0, Sn, mu, dn)
    !write(*,*) 'FYR_partialC1 = ', FYR_C1
    FYR_C2 = FYR_partialC2(l, yf, yr, A, C0, Sn, mu, dn)
    !write(*,*) 'FYR_partialC2 = ', FYR_C2
    FYR_C3 = FYR_partialC3(l, Sn, yf, yr, lam, B, C0, mu, dn)
    !write(*,*) 'FYR_partialC3 = ', FYR_C3
    FYR_C4 = FYR_partialC4(l, Sn, yf, lam, B, C0, D0, mu, dn)
    !write(*,*) 'FYR_partialC4 = ', FYR_C4
    F1C = F1_partialC(l, Sn, yr, mu, C0, FYR_int1, FYR_int2, FYR_int3, FYR_int4, &
      FYR_C1, FYR_C2, FYR_C3, FYR_C4)
    write(*,*) 'F1C = ', F1C
    
    F1D = F1_partialD(l, yr, lam, dn)
    write(*,*) 'F1D = ', F1D
    
    Q1 = Qterm1(l, yf, yr, A, dn)
    !write(*,*) 'Qterm1 = ', Q1
    Q2 = Qterm2(l, yf, yr, A, C0, dn)
    !write(*,*) 'Qterm2 = ', Q2
    Q3 = Qterm3(l, yf, yr, lam, B, C0, dn)
    !write(*,*) 'Qterm3 = ', Q3
    Q4 = Qterm4(l, yf, yr, lam, B, D0, dn)
    !write(*,*) 'Qterm4 = ', Q4
    F20 = F2initial(l, Sn, yf, yr, B, D0, Q1, Q2, Q3, Q4)
    write(*,*) 'F20 = ', F20

    QC1 = Q_partialC1(l, Sn, yf, yr, A, C0, mu, dn)
    !write(*,*) 'QC1 = ', QC1
    QC2 = Q_partialC2(l, Sn, yf, yr, A, C0, mu, dn)
    !write(*,*) 'QC2 = ', QC2
    QC3 = Q_partialC3(l, Sn, yf, yr, lam, B, C0, mu, dn)
    !write(*,*) 'QC3 = ', QC3
    QC4 = Q_partialC4(l, Sn, yf, lam, B, C0, D0, mu, dn)
    !write(*,*) 'QC4 = ', QC4
    F2C = F2_partialC(l, Sn, yf, B, mu, C0, D0, QC1, QC2, QC3, QC4)
    write(*,*) 'F2C = ', F2C

    F2D = F2_partialD(l, yr, lam, dn)
    write(*,*) 'F2D = ', F2D
    

    !---calculation delC and delD @ this iteration---!

    dC = - (F2D * F10 - F1D * F20) / (F1C * F2D - F1D * F2C)
    dD = - (-F2C * F10 + F1C * F20) / (F1C * F2D - F1D * F2C)

    write(*,*) 'delC =', dC
    write(*,*) 'delD =', dD
    write(*,*) ' '

    !Next lamda judgement. If next lam will be Complex number, dC dD must be changed.

    do j = 0, 100
      write(*,*) 'lamda-judge-iteration =', j+1
      C0_n = C0 + dC
      D0_n = D0 + dD
      CaddD = C0_n + D0_n
      A_n    = cal_A(C0_n, Sn, yf)
      lam1_n = cal_lam1(A_n, yf, l)
      lam2_n = cal_lam2(A_n, C0_n, yf, lyfh)
      lam3_n = cal_lam3(A_n, B, C0_n, yf, lyf)
      lam4_n = cal_lam4(lam1_n, lam2_n, lam3_n, Se_l)
      lam_n  = cal_lam(l, lam4, C0_n , D0_n)
      write(*,*) 'C0_n =', C0_n
      write(*,*) 'D0_n =', D0_n
      write(*,*) 'lam4_n =', lam4_n
      write(*,*) 'lam_n =', lam_n
      write(*,*) 'Cn + Dn = ', CaddD

      if (C0_n < 0.0d0) then
        write(*,*) '!!!---------------next_C will be minus---------------!!!'
        dC = - 0.20d0 * dC 
        write(*,*) 'dC is changed'
        write(*,*) 'delC =', dC
      
      else if (D0_n < 0.0d0) then
          write(*,*) '!!!---------------next_D will be minus---------------!!!'
          dD = - 0.20d0 * dD 
          write(*,*) 'dD is changed'
          write(*,*) 'delD =', dD
      
      else if (C0_n >= 0.0d0 .and. lam4_n <= 0.0d0) then
        write(*,*) '!!!--------------next_lam4 will be minus-------------!!!'
        dC = - 0.010d0 * dC 
        write(*,*) 'dC is changed'
        write(*,*) 'delC =', dC

      !else if (lam4_n > 0.0d0 .and. D0_n < - C0_n) then
        !write(*,*) '!!!------lam4 is OK, but next_C+D will be minus------!!!'
        !dD =  -  dD
        !write(*,*) 'dD is changed'
        !write(*,*) 'delD =', dD
      
      else
        exit
      
      end if

    end do

    !---Output variables for result file @ current loop---!

    write(30,*) i+1, dC, dD, C0, D0, lam4, lam

    !---Convergence judgment---!

    !abs(dC) < eps .and. abs(dD) < eps
    !abs(dC) / abs(C0) < eps  .and. abs(dD) / abs(D0) < eps
    
    errC = abs(dC) / abs(C0)
    errD = abs(dD) / abs(D0)

    write(*,*) 'errC = ', errC
    write(*,*) 'errD = ', errD

    if (errC < eps .and. errD < eps) then
      C0 = C0 + dC
      D0 = D0 + dD
      write(*,*) ' '
      write(*,*) '__________________________________________________________'
      write(*,*) 'delD and delC is Converged !!!'
      write(*,*) 'Iteration number =', i+1
      exit

    else if (errC >= eps .and. errD >= eps) then
      C0 = C0 + dC
      D0 = D0 + dD
    
    else if (errC < eps .and. errD >= eps) then
      D0 = D0 + dD
    
    else if (errC >= eps .and.errD < eps) then
      C0 = C0 + dC

    else
      write(*,*) 'dC = ', dC
      write(*,*) 'dD = ', dD
      stop 'something is wrong !!'
    endif

    dC = dC
    dD = dD

    if(i == ite_max) then
      write(*,*) ' '
      write(*,*) '!!!!!----------Iteration time out...----------!!!!!'
    end if

  enddo

  close(30)

  C = C0
  D = D0

  !---output parameter for cal_Ffunction---!

  write(*,*) 'mu = ', mu
  write(*,*) 'B = ', B
  write(*,*) 'M = ', M_inf
  write(*,*) 'h = ', h_inf
  write(*,*) 'W = ', W
  write(*,*) 'l = ', l
  write(*,*) 'yf = ', yf
  write(*,*) 'rho = ', rho_inf
  write(*,*) 'T = ', T_inf
  write(*,*) 'gamma = ', gamma
  write(*,*) 'R = ', R
  write(*,*) 'beta = ', beta
  write(*,*) 'a_inf = ', a_inf
  write(*,*) 'u_inf = ', u_inf
  write(*,*) 'k = ', k
  write(*,*) 's = ', s
  write(*,*) 's2 = ', s2
  write(*,*) 'sn = ', sn
  write(*,*) 'Se_l = ', Se_l

  write(*,*) 'A = ', A
  write(*,*) 'yr = ', yr
  write(*,*) 'lam = ', lam
  write(*,*) 'C = ', C
  write(*,*) 'D = ', D

  open(20, file='F function parameter.txt')

  write(20,'(f16.10)', advance = 'no') l
  write(20,*) '  !---length---!'
  write(20,'(f16.10)', advance = 'no') A
  write(20,*) '  !---paramater A---!'
  write(20,'(f16.10)', advance = 'no') B
  write(20,*) '  !---paramater B---!'
  write(20,'(f16.10)', advance = 'no') C
  write(20,*) '  !---paramater C---!'
  write(20,'(f16.10)', advance = 'no') D
  write(20,*) '  !---paramater D---!'
  write(20,'(f16.10)', advance = 'no') yf
  write(20,*) '  !---paramater yf---!'
  write(20,'(f16.10)', advance = 'no') lam
  write(20,*) '  !---paramater lamda---!'
  write(20,'(i5)', advance = 'no') dn
  write(20,*) '             !---Division number---!'
  write(20,'(f16.10)', advance = 'no') yr
  write(20,*) '  !---paramater yr(%)---!'
  write(20,'(f16.10)', advance = 'no') Se_l
  write(20,*) '  !---calculation result of Se @x=l---!'
  write(20,'(f16.10)', advance = 'no') S
  write(20,*) '  !---slope of balancing point---!'
  write(20,'(f16.10)', advance = 'no') S2
  write(20,*) '  !---slope of balancing point 2---!'
  write(20,'(f16.10)', advance = 'no') Sn
  write(20,*) '  !---slope of balancing point n---!'
  write(20,'(i5)', advance = 'no') i
  write(20,*) '  !---Iteration number of calculation---!'


  close(20)


end program darden