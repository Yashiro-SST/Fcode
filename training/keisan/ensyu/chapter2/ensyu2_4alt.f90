program ensyu2_4alt
  implicit none
  real(8) u(3)
  u(1:3) = (/ 1.0d0, 2.0d0, 3.0d0 /)
  write(*,*) 'l = ', sqrt(dot_product(u,u))
end program ensyu2_4alt