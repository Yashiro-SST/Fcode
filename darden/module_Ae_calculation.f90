module Ae_calculation
  implicit none
  contains
  
    function Ae_cal(x, F, l, dn) Result(Ae)
      integer, intent(in) :: dn
      real(8), intent(in) :: x(dn), F(dn)
      real(8), intent(in) :: l
      real(8) Ae(dn), y, z, dy, sum
      integer i, j

      Ae(:) = 0.0
      dy = l / dble(dn)
      sum = 0.0d0

      do i = 0, dn+1
        do j = 0, i
          y = dy * dble(j)
          z = F(j) * sqrt(x(i)-y)

          if (j == 0 .or. j ==i) then
            sum = sum + 0.5d0 * z
          else
            sum = sum + z
          end if
        end do

        Ae(i) = 4 * sum * dy
      enddo

    end function Ae_cal

end module Ae_calculation