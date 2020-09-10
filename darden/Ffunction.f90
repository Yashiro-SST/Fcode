!==============================================================================!
! program  Culculation   Darden's  F function	and Equivalent Area Destribution !
!	input   :	F function parameter.txt !
!	output　:	Darden Ffunc.txt, Equivalent Area Destribution.txt !
!==============================================================================!
module Ffunc
  implicit none
  contains

    function cal_Ffunc(x, l, A, B, C, D, yf, lam, dn) result(F)
      integer, intent(in) :: dn
      real(8), intent(in) :: x(0:dn)
      real(8), intent(in) :: l, A, B, C, D, yf, lam 
      real(8) F(0:dn), dx
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

module Ae_calculation
  implicit none
  contains
  
    function Ae_cal(x, F, l, dn) Result(Ae)
      integer, intent(in) :: dn
      real(8), intent(in) :: x(0:dn), F(0:dn)
      real(8), intent(in) :: l
      real(8) Ae(0:dn), y, z, dy, sum
      integer i, j

      Ae(:) = 0.0
      dy = l / dble(dn)

      do i = 0, dn
        sum = 0.0d0
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


program Ffunction
  use Ffunc
  use Ae_calculation
  implicit none

  !---Variable definition---!

  real(8) l, A, B, C, D, yf, lam
  real(8) dx
  real(8), allocatable :: x(:), F(:), Ae(:)
  integer dn, i

  !---Reading paramator---!
  open(10,file='F function parameter.txt')

  read(10,*) l       !---Length---!
  read(10,*) A       !---paramater A---!
  read(10,*) B       !---paramater B---!
  read(10,*) C       !---paramater C---!
  read(10,*) D       !---paramater D---!
  read(10,*) yf      !---paramater yf---!
  read(10,*) lam     !---paramater lamda---!
  read(10,*) dn      !---Division number---!

  close(10)

  write(*,*) 'l = ', l
  write(*,*) 'A = ', A
  write(*,*) 'B = ', B
  write(*,*) 'C = ', C
  write(*,*) 'D = ', D
  write(*,*) 'yf = ', yf
  write(*,*) 'lam = ', lam
  write(*,*) 'dn = ', dn
  
  !---Calculate F function---!
  allocate (x(0:dn))
  allocate (F(0:dn))

  x(:) = 0.0d0
  dx = l / dble(dn)
  do i = 0, dn
    x(i) = dx * dble(i)

  enddo

  F(:) = cal_Ffunc(x, l, A, B, C, D, yf, lam, dn)
  write(*,*) 'calculation Ffunc... '
  do i = 0, dn
    write (*,*)i, x(i), F(i)
  end do

  !---Output Ffunction---!

  open(20, file='Darden Ffunc.txt')

  write(20,*) 'i,   x(i),   F(i)'
  do i = 0, dn
    write(20,*) i, x(i), F(i)
  end do

  close(20)
  write(*,*) 'output Ffunc is completed !'

  !---Calculate Equivalent Area---!
  allocate (Ae(0:dn))

  Ae(:) = Ae_cal(x, F, l, dn)
  write(*,*) 'calculation Ae... '
  do i = 0, dn
    write (*,*)i, x(i), Ae(i)
  end do

  !---Output Equivalent Area Destribution---!

  open(30, file='Equivalent Area Destribution.txt')

  write (30,*) 'i,   x(i),   Ae(i)'
  do i = 0, dn
    write(30,*) i, x(i), Ae(i)
  end do

  close(30)
  write(*,*) 'output Ae is completed !'
  write(*,*) ' '

end program Ffunction