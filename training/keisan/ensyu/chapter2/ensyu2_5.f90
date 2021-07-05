program ensyu2_5
  implicit none
  real(8) u(3), l, un(3)
  u(:) = (/ 1.0d0,2.0d0,3.0d0 /)
  l = sqrt(dot_product(dble(u),dble(u)))
  if (l == 0) stop 'vector length is 0...'
  un(1:3) = u(1:3) / dble(l)
  write(*,*) 'normalize u is ', un(1), un(2), un(3)
  write(*,*) 'l = ', l
end program ensyu2_5