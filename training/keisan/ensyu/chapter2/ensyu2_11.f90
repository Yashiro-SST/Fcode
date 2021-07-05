program ensyu2_11
  implicit none
  real(8), allocatable :: a(:)
  integer, allocatable :: x(:)
  integer n
  write(*,'(a)',advance = 'no') 'input n (>= 1) : '
  read(*,*) n
  if (n<1) stop 'stop n < 1'
  allocate (a(n), x(n))
  call random_seed
  call random_number(a(1:n))
  x(1:n) = int(10 * a(1:n))
  write(*,*) 'a = ', a(1:n)
  write(*,*) 'x = ', x(1:n)
  deallocate (a,x)
end program ensyu2_11