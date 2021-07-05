program dotp3
  implicit none
  real(8), allocatable :: u(:), v(:)
  integer :: i, n = 2
  real(8) dotp
  allocate (u(n), v(n))
  !
  write(*,'(a)') 'input vector u'
  read(*,*) u(:)
  write(*,'(a)') 'input vector v'
  read(*,*) v(:)
  dotp = 0.0d0
  do i = 1, n
    dotp = dotp + u(i) * v(i)
  enddo
  write(*,*) 'dot product = ', dotp
  !
  deallocate (u,v)
end program dotp3