program simple_sort
  implicit none
  integer, parameter :: n = 10
  real(8) a(n), am
  integer i, m
  call random_seed
  call random_number(a(:))
  do i = 1, n-1
    am = minval(a(i+1:n))
    m  = minloc(a(i+1:n), 1) + i
    if (a(i) > am) then
      a(m) = a(i)
      a(i) = am
    endif
  enddo
  write(*,*) a(:)
end program simple_sort