program ensyu2_3
  implicit none
  integer a(3), b(3), c(3)
  a(1:3) = (/ 4,5,6 /)
  b(1:3) = (/ 1,2,3 /)
  c(:) = 0
  c(:) = ( a(:)-b(:) ) ** 2
  write(*,*) c(:)
end program ensyu2_3