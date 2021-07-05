program daikei_sekibun
 implicit none
    real(8) dx,x,y,s
    integer i,n
    write(*,'(a)')'input n :'
    read(*,*) n
    if (n < 1) stop 'stop,n<1'
    dx=1.0d0 / dble(n)
    s = 0.0d0
    do i=0,n
        x = dx * dble(i)
        y = x * (1.0d0 -x)
        if (i == 0 .or. i ==n) then
            s = s + 0.5d0 * y
        else
            s = s + y
        end if
    enddo
    s = 6.0d0 * s * dx
    write(*,*)'s = ', s
end program daikei_sekibun