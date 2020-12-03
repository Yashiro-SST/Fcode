program tanx
    implicit none
    integer i, n
    real(8) f, df, x0, x, xn, eps, err

    eps = 1.0 * 10d0**(-8)

    write(*,*) 'input n : '
    read(*,*) n
    write(*,*) 'input x0 : '
    read(*,*) x0

    x = x0

    do i=0, n

        f = tan(x) - 1.0d0 / x
        df = (1 / (cos(x)**2)) + (1.0d0 / (x**2))
        xn = x - f / df
        err = (x - xn) / xn

        if (err <= eps) then
            x = xn
            exit
        end if
        x = xn
    end do

    write(*,*) 'x = ', x

end program tanx