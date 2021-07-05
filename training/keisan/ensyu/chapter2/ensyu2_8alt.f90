program ensyu2_8alt
  implicit none
  real(8) p(3), q(3), c(3), la, lb, lc, s, A
  write(*,'(a)') 'input coordinate of point P'
  read(*,*) p(:)
  write(*,'(a)') 'input coordinate of point Q'
  read(*,*) q(:)
  c(:) = dble(q(:)) - dble(p(:))
  la = sqrt(abs(dot_product(dble(p),dble(p))))
  lb = sqrt(abs(dot_product(dble(q),dble(q))))
  lc = sqrt(abs(dot_product(dble(c),dble(c))))
  s = ( dble(la) + dble(lb) + dble(lc) ) / 2.0d0
  A = sqrt(dble(s) * (dble(s)-dble(la)) * (dble(s)-dble(lb)) * (dble(s)-dble(lc)) )
  write(*,*) 'A is ', A
end program ensyu2_8alt