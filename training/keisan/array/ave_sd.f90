program ensyu2_10
  implicit none
  real(8), allocatable :: r(:), sd(:), bs(:)
  integer n,i
  real(8) ave
  write(*,'(a)',advance = 'no') 'input n (>= 1) : '
  read(*,*) n
  if (n<1) stop 'stop n < 1'
  allocate (r(n), sd(n), bs(n))
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
  bs(1:n) = 0.0d0
  do i = 1, n
    bs(i) = bs(i) + ( (r(i)-ave) ** 2) /n
  enddo
  sd(1:n) = sqrt(bs(1:n))
  write(*,*) 'average, standard deviation = ', ave, sd(1:n)
end program ensyu2_10