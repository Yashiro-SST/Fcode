!=======================================================================!
! program  Culculation of  Darden's  F function								!
!	input	 :	input_darden.txt										!
!	output	:	F.txt, AE.txt												!
!=======================================================================!
program darden
  implicit none

  real(8) A, C, D, lam, yr, m
  real(8) B, mu, yf, W, gamma, M_inf, h_inf, l, R, T_inf, rho_inf, AMW
  real(8) S, k, beta, u_inf, a_inf, Ae_l
  real(8) C0,dC,D0,dD
  real(8) err
  real(8) F1, F2, delta
  real(8) ita_max
  real(8) mach_ang
  real(8) lam1, lam2, lam3
  integer  dn


  !---Read paramator---!
  open(10,file='input_darden.txt')

  read(10,*) mu, B               !---front/rear shock ratio, sig paramata B
  read(10,*) M_inf, h_inf        !---Mach number, Flight height
  read(10,*) W, l, yf            !---Weight, Effective length
  read(10,*) rho_inf, T_inf      !---density, Tenparature @Flight height
  read(10,*) gamma, R, AMW       !---Specific heat ratio，Gas constant, Average molecular weight
  read(10,*) C0, D0              !---initial C, D
  read(10,*) dn                   !---Division number @ integral

  close(10)

  !---Calculate Known function---!

  beta = sqrt(abs(M_inf**2.0 - 1.0))
  a_inf = sqrt(abs(gamma * R * T_inf / (AMW * 0.001)))
  u_inf = M_inf * a_inf
  k = (gamma + 1) * M_inf**2.0 / sqrt(abs(2.0 * beta**3.0))
  S = l / (k * sqrt(abs(h_inf / l)))
  Ae_l = beta * W /(rho_inf * u_inf**2.0)

  !---Cal unknown funcion---!
  
  A = ( C0**2.0 / ( S * yf ) ) - C0 / 2
  
  yr = l + C0 / ( S * mu )

  lam1 = 32.0 * A * (l**2.5) / 15.0 * yf
  lam2 = 32.0 * (C0 - 2.0d0 * A) * (l - yf / 2.0) ** 2.50 / (15.0 * yf)
  lam3 = 16.0 * (B - 2.0 * (C0 - A) / yf) * (l - yf) ** 2.50
  lam = l - (3.0 * ( lam1 + lam2 + lam3 - Ae_l  ) / ((8 * (C0 + D0) ) ** (2.0d0 /3.0d0)))



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

  !---Newton method---!
  
  !---Culculate Equivalent Area---!

end program darden


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

      FYR_int1 = sum

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
        x = dble(nint(yf / 2.0d0)) + dx * dble(i)
        y = sqrt(abs(l - x)) / (yr - x) &
          * (C0 * (2.0d0 * x / yf - 1.0d0) - A * (2.0d0 * x / yf -2.0d0) )

        if (i == 0 .or. i ==n) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      FYR_int2 = sum

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
        x = dble(nint((lam - yf))) + dx * dble(i)
        y = sqrt(abs(l - x)) * (B * (x - yf) + C0) / (yr - x)

        if (i == 0 .or. i ==n) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      FYR_int3 = sum

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
        x = dble(nint((l - lam))) + dx * dble(i)
        y = sqrt(abs(l - x)) * (B * (x - yf) + D0) / (yr - x) 

        if (i == 0 .or. i ==n) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      FYR_int4 = sum

    end function FYR_integral4


    function F10(l, S, yr, B, D0, FYR_int1, FYR_int2, FYR_int3, FYR_int4) Result(F1ini)

      real(8), intent(in) :: l, S, yr, B, D0, FYR_int1, FYR_int2, FYR_int3, FYR_int4
      real(8) FYR, pi, F1ini

      pi = 2.0d0 * acos(0.0d0)
      
      FYR = -(FYR_int1 + FYR_int2 + FYR_int3 + FYR_int4) / (pi * sqrt(abs(yr - l)))
      F1ini = S * (yr - l) + B * (l - yr) - D0 - FYR
    
    end function F10

    function FYR_percialC1(l, yf, yr, A, C0, S, mu, dn) Result(FYR_C1)

      real(8), intent(in) :: l, yf, yr, A, C0, S, mu
      integer, intent(in) :: dn
      real(8) dx, x, y, sum, FYR_C1
      integer i, n

      dx =  l / dble(dn) 
      sum = 0.0d0
      n = nint((yf / 2.0d0) / dx)

      do i = 0, n
        x = dx * dble(i)
        y = 2.0d0 * sqrt(abs(l - x)) / (yf * (yr - x)) &
            * (A / (S * mu *(yr - x)) -(2.0d0 * C0 / (S * yf) - 0.50d0))

        if (i == 0 .or. i ==n) then
          sum = sum + 0.5d0 * y
        else
          sum = sum + y
        end if
      enddo

      FYR_C1 = - sum

    end function FYR_percialC1

    !function FYR_C2
    !end function FYR_C2

    !function FYR_C3
    !end function FYR_C3

    !function FYR_C4
    !end function FYR_C4

    !function FYR_Csum
    !end function FYR_Csum

    !function FYR_D
    !end function FYR_D

    !function Q_C1
    !end function Q_C1

    !function Q_C2
    !end function Q_C2

    !function Q_C3
    !end function Q_C3

    !function Q_C4
    !end function Q_C4

    !function Q_Csum
    !end function Q_Csum

    !unction F2_C
    !end function F2_C

    !function F2_D
    !end function F2_D

!!!0630!!
!!!06/30 17:09!!!
!!!06/30 17:14!!!


end module daikei_sekibun