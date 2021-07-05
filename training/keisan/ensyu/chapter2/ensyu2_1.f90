program ensyu2_1
    implicit none
    integer i
    real(8) u(3), v(3), dotp
    open(10, file = 'input2_1.txt')
    read(10, *) (u(i), i = 1, 3)
    read(10, *) v(1:3)
    close(10)
    write(*, '(3e10.2)') (u(i), i=1,3)
    write(*, '(3e10.2)') v(1:3)
    dotp = 0.0d0
    do i = 1, 3
        dotp = dotp + u(i) * v(i)
    enddo
    write(*,*) 'dot product = ', dotp
end program ensyu2_1