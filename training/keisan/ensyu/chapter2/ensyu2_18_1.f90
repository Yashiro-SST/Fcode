program ensyu2_18_1
  implicit none
  integer n, i, j, is
  real(8), allocatable :: a(:,:), x(:), y(:)
  real(8) t1,t2
  write(*, '(a)') 'input n : '
  read(*,*) n
  allocate (a(n,n),stat = is)
  if(is /= 0) stop 'cannot allocate'
  allocate (x(n))
  call random_seed
  call random_number (a(1:n,1:n))
  call random_number (x(1:n))
  do i=1,n
    write(*,*) 'a(n,n) = ', a(i,1:n)
  enddo
  write(*,*) 'x(n) = ', x(:)
  call cpu_time(t1)
  do i = 1, n
    y(i) = 0.0d0
    do j = 1, n
      y(i) = y(i) + a(i, j) * x(j)
    enddo
  enddo
  call cpu_time(t2)
  write(*,*) 'matmul =', y(:)
  write(*,*) 'cpu time = ', t2 - t1
end program ensyu2_18_1