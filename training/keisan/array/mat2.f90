program mat2
  implicit none
  integer :: n, i, fi = 10
  real(8), allocatable :: a(:,:)
  open(fi, file = 'mat.txt')
  read(fi, *) n
  allocate (a(n,n))
  do i = 1, n
    read(fi, *) a(i, 1:n)
  enddo
  close(fi)
  do i = 1, n
    write(*, '(100e12.4)') a(i, 1:n)
  enddo
  deallocate(a)
end program mat2