program ensyu2_9
  implicit none
  real(8), allocatable :: u(:), v(:)
  integer :: n
  open(10,file = 'input2_9.txt')
  read(10,*) n
  allocate (u(n), v(n))
  read(10,*) u(1:n)
  read(10,*) v(1:n)
  close(10)
  open(20,file = 'output2_9.txt')
  write(20,*) 'dp = ', dot_product(u, v)
  close(20)
  deallocate (u, v)
end program ensyu2_9