program ensyu2_2_1
  implicit none
  integer ia(4), i
  ia(1:4) = (/ 1, 2, 3, 4 /)
  do i =2,4
    ia(i) = ia(i-1)
  enddo
  write(*,*) ia(1:4)
end program ensyu2_2_1