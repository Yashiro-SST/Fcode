module Ffunc
  implicit none
  contains

    function cal_Ffunc(x, l, A, B, C, D, yf, lam, dn) result(F)
      integer, intent(in) :: dn
      real(8), intent(in) :: x(dn)
      real(8), intent(in) :: l, A, B, C, D, yf, lam 
      real(8) F(dn), dx
      integer i

      F(:) = 0.0d0
      dx = l / dble(dn)
    
      do i = 0, dn
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

    end function cal_Ffunc
    
end module Ffunc