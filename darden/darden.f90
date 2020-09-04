!=======================================================================!
! program  Culculation of  Darden's  F function								!
!	input	 :	input_darden.txt										!
!	output	:	F.txt, AE.txt												!
!=======================================================================!
program darden
  use daikei_sekibun
  use Ffunc
  use Ae_calculation
  implicit none

  !---Variable definition---!

  real(8) A, C, D, lam, yr
  real(8) B, mu, yf, W, gamma, M_inf, h_inf, l, R, T_inf, rho_inf, AMW
  real(8) S, k, beta, u_inf, a_inf, Ae_l
  real(8) C0,dC,D0,dD
  real(8) err
  real(8) lam1, lam2, lam3
  real(8) dx
  real(8) F10, F20, F1C, F2C, F1D, F2D
  real(8) FYR_int1, FYR_int2, FYR_int3, FYR_int4
  real(8) FYR_C1, FYR_C2, FYR_C3, FYR_C4
  real(8) Q1, Q2, Q3, Q4
  real(8) QC1, QC2, QC3, QC4
  real(8), allocatable :: x(:), F(:), Ae(:)
  integer dn, i, ite_max

  !---Reading paramator---!
  open(10,file='input_darden.txt')

  read(10,*) mu, B               !---front/rear shock ratio, sig paramata B
  read(10,*) M_inf, h_inf        !---Mach number, Flight height
  read(10,*) W, l, yf            !---Weight, Effective length, nose length
  read(10,*) rho_inf, T_inf      !---density, Tenparature @Flight height
  read(10,*) gamma, R, AMW       !---Specific heat ratio，Gas constant, Average molecular weight
  read(10,*) C0, D0              !---initial C, D
  read(10,*) dn                  !---Division number @ integral
  read(10,*) ite_max             !---Maximum number of iterations

  close(10)

  !---Calculation Known function---!

  beta = sqrt(abs(M_inf**2.0d0 - 1.0d0))
  a_inf = sqrt(abs(gamma * R * T_inf / (AMW * 0.001)))
  u_inf = M_inf * a_inf
  k = (gamma + 1) * M_inf**4.0d0 / sqrt(2.0d0 * beta**3.0d0)
  S = l / (k * sqrt(abs(h_inf / l)))
  Ae_l = beta * W * 9.80665 /(rho_inf * u_inf**2.0d0)

  !---Output---!

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
  write(*,*) 'AMW = ', AMW
  write(*,*) 'beta = ', beta
  write(*,*) 'a_inf = ', a_inf
  write(*,*) 'u_inf = ', u_inf
  write(*,*) 'k = ', k
  write(*,*) 's = ', s
  write(*,*) 'Ae_l = ', Ae_l
  write(*,*) 'C0 = ', C0
  write(*,*) 'D0 = ', D0


  !---Newton method---!
  err = 1.0d0 * 10d0**(-4d0)

  do i = 0, ite_max

    write(*,*) '_____________________________'

    write(*,*) 'iteration' , i+1
    write(*,*) 'C0 = ', C0
    write(*,*) 'D0 = ', D0

    !---Calculation unknown funcion---!
    A = ( C0**2.0d0 / ( S * yf ) ) - C0 / 2.0d0
    yr = l + C0 / ( S * mu )
    lam1 = 32d0 * A * (l**2.5) / (15d0 * yf)
    lam2 = 32d0 * (C - 2.0d0 * A) * ((l - yf / 2.0d0) ** 2.5d0) / (15d0 * yf)
    lam3 = 16d0 * (2.0d0 * A + B * yf - 2.0d0 * C) * ((l - yf) ** 2.5d0) / (15d0 * yf)
    lam = l - (3.0d0 * (lam1 + lam2 + lam3 - Ae_l) / (8.0d0 * (C0 + D0))) ** (2.0d0 / 3.0d0)

    write(*,*) 'A = ', A
    write(*,*) 'yr = ', yr
    write(*,*) 'lam1 = ', lam1
    write(*,*) 'lam2 = ', lam2
    write(*,*) 'lam3 = ', lam3
    write(*,*) 'lam = ', lam

    !---Calculation F10, F20, F1 partial CorD, and F2 partial CorD---!

    FYR_int1 = FYR_integral1(l, yf, yr, A, dn)
    !write(*,*) 'FYR_int1 = ', FYR_int1
    FYR_int2 = FYR_integral2(l, yf, yr, A, C0, dn)
    !write(*,*) 'FYR_int2 = ', FYR_int2
    FYR_int3 = FYR_integral3(l, yf, yr, lam, B, C0, dn)
    !write(*,*) 'FYR_int3 = ', FYR_int3
    FYR_int4 = FYR_integral4(l, yf, yr, lam, B, D0, dn)
    !write(*,*) 'FYR_int4 = ', FYR_int4
    F10 = F1initial(l, S, yr, B, D0, FYR_int1, FYR_int2, FYR_int3, FYR_int4)
    !write(*,*) 'F10 = ', F10

    FYR_C1 = FYR_partialC1(l, yf, yr, A, C0, S, mu, dn)
    !write(*,*) 'FYR_partialC1 = ', FYR_C1
    FYR_C2 = FYR_partialC2(l, yf, yr, A, C0, S, mu, dn)
    !write(*,*) 'FYR_partialC2 = ', FYR_C2
    FYR_C3 = FYR_partialC3(l, S, yf, yr, lam, B, C0, mu, dn)
    !write(*,*) 'FYR_partialC3 = ', FYR_C3
    FYR_C4 = FYR_partialC4(l, S, yf, yr, lam, B, D0, mu, dn)
    !write(*,*) 'FYR_partialC4 = ', FYR_C4
    F1C = F1_partialC(l, S, yr, B, mu, C0, FYR_int1, FYR_int2, FYR_int3, FYR_int4, &
                      FYR_C1, FYR_C2, FYR_C3, FYR_C4)
    !write(*,*) 'F1C = ', F1C
    
    F1D = F1_partialD(l, yr, lam, dn)
    !write(*,*) 'F1D = ', F1D
    
    Q1 = Qterm1(l, yf, yr, A, dn)
    !write(*,*) 'Qterm1 = ', Q1
    Q2 = Qterm2(l, yf, yr, A, C0, dn)
    !write(*,*) 'Qterm2 = ', Q2
    Q3 = Qterm3(l, yf, yr, lam, B, C0, dn)
    !write(*,*) 'Qterm3 = ', Q3
    Q4 = Qterm4(l, yf, yr, lam, B, D0, dn)
    !write(*,*) 'Qterm4 = ', Q4
    F20 = F2initial(l, S, yr, B, D0, Q1, Q2, Q3, Q4)
    !write(*,*) 'F20 = ', F20

    QC1 = Q_partialC1(l, S, yf, yr, A, C0, mu, dn)
    !write(*,*) 'QC1 = ', QC1
    QC2 = Q_partialC2(l, S, yf, yr, A, C0, mu, dn)
    !write(*,*) 'QC2 = ', QC2
    QC3 = Q_partialC3(l, S, yf, yr, lam, B, C0, mu, dn)
    !write(*,*) 'QC3 = ', QC3
    QC4 = Q_partialC4(l, S, yf, lam, B, C0, D0, mu, dn)
    !write(*,*) 'QC4 = ', QC4
    F2C = F2_partialC(l, S, B, mu, C0, D0, QC1, QC2, QC3, QC4)
    !write(*,*) 'F2C = ', F2C

    F2D = F2_partialD(l, S, yr, lam, dn)
    !write(*,*) 'F2D = ', F2D
    

    !---calculation delC and delD @ this iteration---!

    dC = - (F2D * F10 - F1D * F20) / (F1C * F2D - F1D * F2C)
    dD = - (-F2C * F10 + F1C * F20) / (F1C * F2D - F1D * F2C)

    write(*,*) 'delC =', dC
    write(*,*) 'delD =', dD

    !---Convergence judgment---!

    if (abs(dC) < err .and. abs(dD) < err) then
      write(*,*) ' '
      write(*,*) 'delD and delC is Converged!!!'
      !write(*,*) 'delC =', dC
      !write(*,*) 'delD =', dD
      write(*,*) 'Iteration number =', i
      exit

    else if (abs(dC) >= err .and. abs(dD) >= err) then
      C0 = C0 + dC
      D0 = D0 + dD
    
    else if (abs(dC) < err .and. abs(dD) >= err) then
      D0 = D0 + dD
    
    else if (abs(dC) >= err .and. abs(dD) < err) then
      C0 = C0 + dC

    else
      write(*,*) 'dC = ', dC
      write(*,*) 'dD = ', dD
      stop 'something is wrong !!'
    endif

  enddo

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
  write(*,*) 'AMW = ', AMW
  write(*,*) 'beta = ', beta
  write(*,*) 'a_inf = ', a_inf
  write(*,*) 'u_inf = ', u_inf
  write(*,*) 'k = ', k
  write(*,*) 's = ', s
  write(*,*) 'Ae_l = ', Ae_l

  write(*,*) 'A = ', A
  write(*,*) 'yr = ', yr
  write(*,*) 'lam = ', lam
  write(*,*) 'C = ', C
  write(*,*) 'D = ', D

  open(20, file='F function parameter.txt')
  write(20,*) l
  write(20,*) A
  write(20,*) B
  write(20,*) C
  write(20,*) D
  write(20,*) yf
  write(20,*) lam
  write(20,*) dn
  close(20)
  
  !---Calculate F function---!
  allocate (x(dn))
  allocate (F(dn))

  x(:) = 0.0d0
  dx = l / dble(dn)
  do i = 0, dn
    x(i) = dx * dble(i)

  enddo

  F(:) = cal_Ffunc(x, l, A, B, C, D, yf, lam, dn)

  open(30, file='Darden Ffunc.txt')

      write (30,*) 'i,   x(i),   F(i)'
      do i = 0, dn
        write (30,*)i, x(i), F(i)
      end do

  close(30)
  write(*,*) 'output Ffunc is completed !'

  !---Calculate Equivalent Area---!
  allocate (Ae(dn))

  Ae(:) = Ae_cal(x, F, l, dn)

  open(40, file='Equivalent Area Destribution.txt')

      write (40,*) 'i,   x(i),   Ae(i)'
      do i = 0, dn
        write (40,*)i, x(i), Ae(i)
      end do

  close(40)
  write(*,*) 'output Ae is completed !'
  write(*,*) ' '

end program darden



!---module function of F1, F2---!

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

      do i = 0, n-1
        x = lam + dx * dble(i)
        t = atan(sqrt((yr - l) / (l - x)))
        y = (B * (x - yf) - D0) * t

        if (i == 0 .or. i == n-1) then
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

      do i = 0, n-1
        x = lam + dx * dble(i)
        tc = 1.0d0 / (2.0d0 * sqrt(S * mu * C0 * (l - x)) + C0)
        y = (B * (x - yf) - D0) * tc

        if (i == 0 .or. i == n-1) then
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
      F2D = (- 1.0d0 - S / 2.0d0) * (yr - l) - QD

    end function F2_partialD

end module daikei_sekibun


module Ffunc
  implicit none
  contains

    function cal_Ffunc(x, l, A, B, C, D, yf, lam, dn) result(F)
      integer, intent(in) :: dn
      real(8), intent(in) :: x(dn)
      real(8), intent(in) :: l, A, B, C, D, yf, lam 
      real(8) F(dn), dx
      integer i

      F(:) = 0.0d0
      dx = l / dble(dn)
    
      do i = 0, dn+1
        if(x(i) <= (yf / 2.0d0)) then
          F(i) = 2.0d0 * x(i) * A / yf
        
        else if((yf / 2.0d0) < x(i) .and. x(i) < yf ) then
          F(i) = C * (2.0d0 * x(i) / yf - 1.0d0) - A * (2.0d0 * x(i) / yf - 2.0d0)
    
        else if(yf <= x(i) .and. x(i) < lam) then 
          F(i) = B * (x(i) - yf) + C
    
        else if(lam <= x(i) .and. x(i) < l) then
          F(i) = B * (x(i) - yf) - D
        end if
      enddo

    end function cal_Ffunc

end module Ffunc

module Ae_calculation
  implicit none
  contains
  
    function Ae_cal(x, F, l, dn) Result(Ae)
      integer, intent(in) :: dn
      real(8), intent(in) :: x(dn), F(dn)
      real(8), intent(in) :: l
      real(8) Ae(dn), y, z, dy, sum
      integer i, j

      Ae(:) = 0.0
      dy = l / dble(dn)
      sum = 0.0d0

      do i = 0, dn+1
        do j = 0, i
          y = dy * dble(j)
          z = F(j) * sqrt(x(i)-y)

          if (j == 0 .or. j ==i) then
            sum = sum + 0.5d0 * z
          else
            sum = sum + z
          end if
        end do

        Ae(i) = 4 * sum * dy
      enddo

    end function Ae_cal

end module Ae_calculation