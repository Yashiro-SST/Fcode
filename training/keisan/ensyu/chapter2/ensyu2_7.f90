program ensyu2_7
  implicit none
  real(8) u(3), v(3),  i, j, k
  write(*,'(a)') 'input vector u'
  read(*,*) u(:)
  write(*,'(a)') 'input vector v'
  read(*,*) v(:)
  i = u(2)*v(3) - u(3)*v(2)
  j = u(3)*v(1) - u(1)*v(3)
  k = u(1)*v(2) - u(2)*v(1)
  write(*,*) 'u cross v = ',i, j, k
end program ensyu2_7