program ensyu2_2_2
  implicit none
  integer ia(4)
  ia(1:4) = (/ 1, 2, 3, 4 /)
  ia(2:4) = ia(1:3)
  write(*,*) ia(1:4)
end program ensyu2_2_2