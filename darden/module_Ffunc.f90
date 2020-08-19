module cal_Ffunc
  implicit none

  contains

  real(8) A, C, D, lam, yr
  real(8), allocatable :: x(:), F(:)
  integer dn, i

  open(20, file = 'F_function.txt')

  x(:) = 0.0d0
  F(:) = 0.0d0
  
  dx = l / dble(dn)
  allocate (x(dn))

  do i = 0, dn
    x(i) = dx * dble(i)

    if(x(i) <= (yf / 2.0d0)) then
      F(i) = 2.0d0 * x(i) * A / yf
    
    else if((yf / 2.0d0) < x(i) .and. x(i) < yf ) then
      F(i) = C * (2.0d0 * x(i) / yf - 1.0d0) - A * (2.0d0 * x(i) / yf - 2.0d0)

    else if(yf <= x(i) .and. x(i) < lam) then 
      F(i) = B * (x(i) - yf) + C

    else if(lam <= x(i) .and. x(i) < l) then
      F(i) = B * (x(i) - yf) - D
    end if

  enddo

  write (20,*) 'i,   x(i),   F(i)'
  do i = 0, dn
    write (20,*) i, x(i), F(i)
  end do

  close(20)

end module cal_Ffunc