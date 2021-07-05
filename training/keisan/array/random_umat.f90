program random_umat
  implicit none
  real(8), allocatable :: a(:,:)
  integer n, i, j
  write(*,'(a)', advance = 'no') 'input n (1<=n<=100)'
  read(*,*) n
  if (n<1 .or. 100<n) stop 'stop, n is invalid...'
  allocate (a(n,n))
  call random_seed
  do j = 1, n
    call random_number(a(1:j, j))
    a(j+1:n, j) = 0.0d0
  enddo
  do i = 1, n
    write(*,'(100e12.4)') a(i, 1:n)
  enddo
end program random_umat