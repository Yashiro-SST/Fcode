program ensyu2_10
  implicit none
  real(8), allocatable :: r(:)
  integer n,i
  real(8) ave, bs, sd
  write(*,'(a)',advance = 'no') 'input n (>= 1) : '
  read(*,*) n
  if (n<1) stop 'stop n < 1'
  allocate (r(n))
  call random_seed
  call random_number(r(1:n))
  r(1:n) = 2.0d0 * r(1:n) - 1.0d0
  write(*,*) 'ramdom number = ', r(1:n)
  !
  ave = 0.0d0
  do i = 1, n
    ave = ave + r(i)
  enddo
  ave = ave / n
  bs = 0.0d0
  do i = 1, n
    bs = bs + ( (r(i)-ave) ** 2) /n
  enddo
  sd = sqrt(dble(bs))
  write(*,*) 'ave, variance, sd = ', ave, bs, sd
end program ensyu2_10