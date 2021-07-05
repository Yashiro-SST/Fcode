program ensyu2_4
  implicit none
  real(8) u(3), l
  u(1:3) = (/ 1.0d0, 2.0d0, 3.0d0 /)
  l = sqrt( u(1)**2 + u(2)**2 + u(3)**2 )
  write(*,*) 'l = ', l
end program ensyu2_4