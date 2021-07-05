program ensyu2_6
  implicit none
  real(8) u(3), l, yogen(3), sqar
  integer :: x(3)=(/ 1,0,0 /), y(3)=(/ 0,1,0 /), z(3)=(/ 0,0,1 /)
  write(*,'(a)') 'input vector u'
  read (*,*) u(:)
  l = sqrt(dot_product(dble(u),dble(u)))
  if (l == 0) stop 'vector length is 0...'
  yogen(1) = dot_product(dble(u),x) / dble(l)
  yogen(2) = dot_product(dble(u),y) / dble(l)
  yogen(3) = dot_product(dble(u),z) / dble(l)
  write(*,*) 'l = ',l
  write(*,*) 'cos(a,b,g) = ',yogen(1), yogen(2), yogen(3)
  sqar = yogen(1)**2 + yogen(2)**2 + yogen(3)**2
  write(*,*) 'cosa**2 + cosb**2 + cosg**2 = ',sqar
end program ensyu2_6