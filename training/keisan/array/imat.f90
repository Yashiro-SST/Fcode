program imat
  implicit none
  integer a(2, 2)
  a(1, 1:2) = (/ 11, 12 /)
  a(2, 1:2) = (/ 21, 22 /)
  write(*,*) a(:,:)
end program imat