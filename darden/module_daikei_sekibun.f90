module daikei_sekibun
  implicit none
  contains

    function FYR_integral1(l, yf, yr, A, dn) Result(FYR_int1)
      
      real(8), intent(in) :: l, yf, yr, A
      integer, intent(in) :: dn
      real(8) dx, x, y, sum, FYR_int1
      integer i, n

      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((yf / 2.0d0) / dx)

      do i = 0, n
        x = dx * dble(i)
        y = 2.0d0 * x * sqrt(abs(l - x)) * A / (yf * (yr - x))

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
      real(8) dx, x, y, sum, FYR_int2
      integer i, n

      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((yf / 2.0d0) / dx)

      do i = 0, n
        x = yf / 2.0d0 + dx * dble(i)
        y = sqrt(abs(l - x)) / (yr - x) &
          * (C0 * (2.0d0 * x / yf - 1.0d0) - A * (2.0d0 * x / yf -2.0d0) )

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
      real(8) dx, x, y, sum, FYR_int3
      integer i, n

      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((lam - yf) / dx)

      do i = 0, n
        x = yf + dx * dble(i)
        y = sqrt(abs(l - x)) * (B * (x - yf) + C0) / (yr - x)

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
      real(8) dx, x, y, sum, FYR_int4
      integer i, n

      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((l - lam) / dx)

      do i = 0, n
        x = lam + dx * dble(i)
        y = sqrt(abs(l - x)) * (B * (x - yf) - D0) / (yr - x) 

        if (i == 0 .or. i ==n) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      FYR_int4 = sum * dx

    end function FYR_integral4

    function F1initial(l, S, yr, B, D0, FYR_int1, FYR_int2, FYR_int3, FYR_int4) Result(F10)

      real(8), intent(in) :: l, S, yr, B, D0, FYR_int1, FYR_int2, FYR_int3, FYR_int4
      real(8) FYR, pi, F10

      pi = 2.0d0 * acos(0.0d0)
      
      FYR = -(FYR_int1 + FYR_int2 + FYR_int3 + FYR_int4) / (pi * sqrt(abs(yr - l)))
      F10 = S * (yr - l) + B * (l - yr) - D0 - FYR
    
    end function F1initial

    function FYR_partialC1(l, yf, yr, A, C0, S, mu, dn) Result(FYR_C1)

      real(8), intent(in) :: l, yf, yr, A, C0, S, mu
      integer, intent(in) :: dn
      real(8) dx, x, y, sum, FYR_C1
      integer i, n

      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((yf / 2.0d0) / dx)

      do i = 0, n
        x = dx * dble(i)
        y = - 2.0d0 * x * sqrt(abs(l - x)) / (yf * (yr - x)) &
            * (A / (S * mu *(yr - x)) -(2.0d0 * C0 / (S * yf) - 0.50d0))

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
      real(8) dx, x, y, sum, FYR_C2
      integer i, n

      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((yf / 2.0d0) / dx)

      do i = 0, n
        x = (yf / 2.0d0) + dx * dble(i)
        y = - sqrt(abs(l - x)) / (yr - x) * (1.0d0 / (S * mu * (yr - x)) &
            * (C0 * (2.0d0 * x / yf - 1.0d0) - A * (2.0d0 * x / yf - 2.0d0)) &
            + (2.0d0 * x / yf - 1.0d0) - (2.0d0 * C0 / (S * yf) - 0.50d0) &
            * (2.0d0 * x / yf - 2.0d0))

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
      real(8) dx, x, y, sum, FYR_C3
      integer i, n

      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((lam - yf) / dx)

      do i = 0, n
        x = yf + dx * dble(i)
        y = - sqrt(abs(l - x)) / (yr - x) * ((B * (x - yf) + C0) &
            / (S * mu * (yr - x)) + 1.0d0)

        if (i == 0 .or. i ==n) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      FYR_C3 = sum * dx

    end function FYR_partialC3

    function FYR_partialC4(l, S, yf, yr, lam, B, D0, mu, dn) Result(FYR_C4)
      
      real(8), intent(in) :: l, S, yf, yr, lam, B, mu, D0
      integer, intent(in) :: dn
      real(8) dx, x, y, sum, FYR_C4
      integer i, n

      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((l - lam) / dx)

      do i = 0, n
        x = lam + dx * dble(i)
        y = sqrt(abs(l - x)) / (S * mu * (yr - x)**2.0d0) * (B * (x - yf) - D0)  

        if (i == 0 .or. i ==n) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      FYR_C4 = sum * dx

    end function FYR_partialC4

    function F1_partialC(l, S, yr, B, mu, C0, FYR_int1, FYR_int2, FYR_int3, FYR_int4, &
      FYR_C1, FYR_C2, FYR_C3, FYR_C4) Result(F1C)
      
      real(8), intent(in) :: l, S, yr, B, mu, C0, FYR_int1, FYR_int2, FYR_int3, FYR_int4, &
      FYR_C1, FYR_C2, FYR_C3, FYR_C4
      real(8) pi, FYRC, F1C

      pi = 2.0d0 * acos(0.0d0)
      FYRC = - sqrt(S * mu) * C0**(-1.50d0) &
            * (FYR_int1 + FYR_int2 + FYR_int3 + FYR_int4) / 2.0d0 * pi &
            + (FYR_C1 + FYR_C2 + FYR_C3 + FYR_C4) / (pi * sqrt(yr - l))
      F1C = 1.0d0 / mu + B / (S * mu) - FYRC 

    end function F1_partialC

    function F1_partialD(l, yr, lam, dn) Result(F1D)

      real(8), intent(in) :: l, yr, lam
      integer, intent(in) :: dn
      real(8) dx, x, y, sum, F1D
      integer i, n

      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((l - lam) / dx)

      do i = 0, n
        x = lam + dx * dble(i)
        y = - sqrt(abs(l - x)) / (yr - x) 

        if (i == 0 .or. i ==n) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      F1D = - 1.0d0 - sum * dx

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
      real(8) dx, x, y, t, sum, pi, Q4
      integer i, n

      pi = 2.0d0 * acos(0.0d0)
      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((l - lam) / dx)

      do i = 0, n
        x = lam + dx * dble(i)
        t = atan(sqrt((yr - l) / (l - x)))
        y = (B * (x - yf) - D0) * t

        if (i == 0 .or. i == n) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      Q4 = - 2.0d0 / pi * sum * dx

    end function Qterm4

    function F2initial(l, S, yr, B, D0, Q1, Q2, Q3, Q4) Result(F20)

      real(8), intent(in) :: l, S, yr, B, D0, Q1, Q2, Q3, Q4
      real(8) Q_sum, F20
      
      Q_sum = Q1 + Q2 + Q3 + Q4
      F20 = (B * (l - yr) - D0 + S / 2.0d0 *(yr - D0)) * (yr - l) - Q_sum

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
        tc = 1.0d0 / (2.0d0 * sqrt(S * mu * C0 * (l - x)) + C0)
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
        tc = 1.0d0 / (2.0d0 * sqrt(S * mu * C0 * (l - x)) + C0)
        y = (2.0d0 * x / yf - 1.0d0) - (2.0d0 * C0 / (S * yf) - 0.5d0) &
            * (2.0d0 * x / yf - 2.0d0) * t &
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
        tc = 1.0d0 / (2.0d0 * sqrt(S * mu * C0 * (l - x)) + C0)
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
      real(8) dx, x, y, tc, sum, pi, QC4
      integer i, n

      pi = 2.0d0 * acos(0.0d0)
      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((l - lam) / dx)

      do i = 0, n
        x = lam + dx * dble(i)
        tc = 1.0d0 / (2.0d0 * sqrt(S * mu * C0 * (l - x)) + C0)
        y = (B * (x - yf) - D0) * tc

        if (i == 0 .or. i == n) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      QC4 = - 2.0d0 / pi * sum * dx

    end function Q_partialC4

    function F2_partialC(l, S, B, mu, C0, D0, QC1, QC2, QC3, QC4) result(F2C)
      real(8), intent(in) :: l, S, B, mu, C0, D0, QC1, QC2, QC3, QC4
      real(8) QC_sum, smu, F2C

      smu = S * mu
      QC_sum = QC1 + QC2 + QC3 + QC4
      F2C = C0 * (B / smu + (S / (2.0d0 * smu)) * (1.0d0 / smu + 1.0d0)) &
            - D0 + S * (l - D0) / 2.0d0 - QC_sum

    end function F2_partialC

    function F2_partialD(l, S, yr, lam, dn) Result(F2D)

      real(8), intent(in) :: l, S, yr, lam
      integer, intent(in) :: dn
      real(8) dx, x, y, sum, pi, QD, F2D
      integer i, n

      pi = 2.0d0 * acos(0.0d0)
      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((l - lam) / dx)

      do i = 0, n
        x = lam + dx * dble(i)
        y = atan(sqrt((yr - l) / (l - x))) 

        if (i == 0 .or. i == n) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      QD = 2.0d0 / pi * sum * dx
      F2D = (- 1.0d0 - S / 2.0d0) * (yr - l) - QD

    end function F2_partialD

    !function Ae_cal(l, A, B, C, D, yf, dn) Result(Ae)
      !integer, intent(in) :: dn
      !real(8), intent(in) :: l, A, B, C, D, yf
      !real(8) Ae(n)

      !dy = 0.0d0

    !end function Ae_cal


end module daikei_sekibun