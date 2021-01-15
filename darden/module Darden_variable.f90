module Darden_variable
  implicit none
  contains

    function cal_beta(M_inf) Result(beta)
      real(8), intent(in) :: M_inf
      real(8) beta

      beta = sqrt(abs(M_inf**2.0d0 - 1.0d0))

    end function cal_beta

    function cal_ainf(gamma, R, T_inf) Result(a_inf)
      real(8), intent(in) :: gamma, R, T_inf
      real(8) a_inf

      a_inf = sqrt(abs(gamma * R * T_inf))

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

    function cal_S(k, h_inf, l) Result(S)
      real(8), intent(in) :: k, h_inf, l
      real(8) S

      S = 1.0d0 / (k * sqrt(abs(h_inf / l)))

    end function cal_S

    function cal_S2(gamma, M_inf, beta, dn) Result(S2)
      real(8), intent(in) :: gamma, M_inf, beta
      real(8) S2, gamma_l, sum, S2_int
      integer, intent(in) :: dn
      integer i, n

     ! gamma_l = (gamma + 1) / 2

      do i = 0, n
      !  x = dx * dble(i)
      !  K = sqrt(l - x) / (yr - x)
      !  y = 2.0d0 * x / yf * A * K

      !  if (i == 0 .or. i ==n) then
      !    sum = sum + 0.5d0 * y
       ! else
      !    sum = sum + y
      !  end if
      enddo

      S2_int = sum * dx

      S2 = sqrt(2.0d0 * beta) / (gamma_l * M_inf**3.0d0 * S2_int)
    
    end function cal_S2

    function cal_Se_l(beta, W, rho_inf, u_inf, unit) Result(Ae_l)
      real(8), intent(in) :: beta, W, rho_inf, u_inf
      real(8) g, Ae_l
      integer, intent(in) :: unit

      if (unit == 1) then
        g = 9.80665
      else if (unit == 0) then
        g = g * 3.28084
      else
        stop 'cal_Ae error, invalid unit system number is entered!!!'
      end if

      Ae_l = beta * W * g /(rho_inf * u_inf**2.0d0)

    end function cal_Se_l

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

    function cal_lam4(lam1, lam2, lam3, Ae_l) Result(lam4)
      real(8), intent(in) :: lam1, lam2, lam3, Ae_l
      real(8) lam4

      lam4 = lam1 + lam2 + lam3 - Ae_l

    end function cal_lam4

    function cal_lam(l, lam4,C0, D0) Result(lam)
      real(8), intent(in) :: l, lam4,C0, D0
      real(8) lam

      lam = l - (3.0d0 * (lam4) / (8.0d0 * (C0 + D0))) ** (2.0d0 / 3.0d0)

    end function cal_lam
  
  end module Darden_variable