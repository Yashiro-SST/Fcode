program ensyu2_8
  implicit none
  real(8) p(3), q(3), c(3), A, l
  write(*,'(a)') 'input coordinate of point P'
  read(*,*) p(:)
  write(*,'(a)') 'input coordinate of point Q'
  read(*,*) q(:)
  c(1) = dble(p(2))*dble(q(3)) - dble(p(3))*dble(q(2))
  c(2) = dble(p(3))*dble(q(1)) - dble(p(1))*dble(q(3))
  c(3) = dble(p(1))*dble(q(2)) - dble(p(2))*dble(q(1))
  l = sqrt(abs(dot_product(dble(c),dble(c))))
  A = dble(l) / 2.0d0
  write(*,*) 'A is ', A
end program ensyu2_8