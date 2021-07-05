program ensyu2_3alt
  implicit none
  integer a(3), b(3), c(3), i
  a(1:3) = (/ 4,5,6 /)
  b(1:3) = (/ 1,2,3 /)
  c(:) = 0
  do i = 1, 3
    c(i) = ( a(i)-b(i) ) ** 2
  enddo
  write(*,*) c(:)
end program ensyu2_3alt