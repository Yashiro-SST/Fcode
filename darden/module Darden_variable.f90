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
        T_slope = 0.0065 / 3.28084
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
      
      Mh = M_inf         !---if hight where cal initial wave form is higher from zc 
      Th = Tc            !---same as above
      gamma_l = (gamma + 1) / 2

      !---integral in the ray tube ratio formura---!

      dz = h_inf / dn
      z_sp = h_inf - zc
      sum = 0.0d0

      do j = 0, dn
        !---z is vertical distance of signal from airplane axis---!
        z = dz * dble(j)
        if (0 <= z .and. z <= z_sp) then               !---temparature value at Stratospere
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
        if (0 <= z .and. z <= z_sp) then             !---temparature value at Stratospere
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