program dotp2
  implicit none
  integer i
  integer, parameter :: n = 2
  real(8) u(n), v(n), dotp
  write(*,'(a)') 'input vector u'
  read(*,*) u(:)
  write(*,'(a)') 'input vector v'
  read(*,*) v(:)
  dotp = 0.0d0
  do i = 1, n
    dotp = dotp + u(i) * v(i)
  enddo
  write(*,*) 'dot product = ', dotp
end program dotp2