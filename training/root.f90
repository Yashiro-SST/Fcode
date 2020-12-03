program root
    implicit none
    integer i, n, C
    real(8) x0, x, xn, eps, err

    eps = 1.0 * 10d0**(-8)

    write(*,*) 'input C : '
    read(*,*) C
    write(*,*) 'input n : '
    read(*,*) n
    write(*,*) 'input x0 : '
    read(*,*) x0

    x = x0

    do i=0, n
        xn = (x + C / x) / 2.0d0
        err = (x - xn) / xn

        if (err <= eps) then
            x = xn
            exit
        end if
        x = xn
    end do

    write(*,*) 'square root of C = ', x

end program root