!=======================================================================!
! program  Interpolation of Volume destribution			              			!
!	input	 :	sl_input.txt, sv_input.txt								              		!
!	output	:	sv_interpolated.txt										    	              	!
!=======================================================================!

Module basic_module

  Implicit none
  Contains
  !=================================================
  Function IP(fx, fy, lx, ly, x) Result(y)
  Implicit none
  Double precision fx, fy, lx, ly, x, y
  
    y = fy + (x - fx) * (ly - fy) / (lx - fx)
  
  End Function IP
  !=================================================
  
  !=================================================
  Function OPl(fx, fy, lx, ly, x) Result(y)
  Implicit none
  Double precision fx, fy, lx, ly, x, y
    
    y = ly + (x - lx) * (ly - fy) / (lx - fx)
    
  End Function OPl
  !=================================================
  
  !=================================================
  Function OPf(fx, fy, lx, ly, x) Result(y)
  Implicit none
  Double precision fx, fy, lx, ly, x, y
    
    y = fy + (x - fx) * (ly - fy) / (lx - fx)
    
  End Function OPf
  !=================================================
  
End Module basic_module
  
Program Linear_Interpolation
  Use basic_module
  Implicit none
  Integer i, j
  Integer n1, n2
  real(8) n
  real(8) dum
  real(8),allocatable :: x1(:), y1(:)
  real(8),allocatable :: x2(:), y2(:)
  real(8),allocatable :: iy2(:)
  real(8) gr, l, alpha, Ma, md, deg, pi, lv, lr
  gr = 1000
  l = 30.77
  pi = acos(-1.0d0)
  
  write(*,*) 'alpha of level flight ='
  read(*,*) alpha
  write(*,*) 'Flight Mach number ='
  read(*,*) Ma
  md = asin(1/Ma)*180/pi
  write(*,*) 'Mach degree = ',md
  deg = alpha + md
  write(*,*) 'deg = ', deg
  write(*,*) 'Vertical length (grid dimension) ='
  read(*,*) lv

  !一個目の2次元プロットデータ読み込み　!揚力等価断面積読み込み
    Open(100,file="sl_input.txt")
    n1 = 0
    Do
      Read (100, *, end=100) dum
      n1 = n1 + 1
    End do
  
  100 Continue
    Rewind(100)
    Allocate(x1(n1), y1(n1), iy2(n1))
    Do i = 1, n1
      Read(100,*) x1(i), y1(i)
    End do
    Close(100)

    write(*,*) 'sl_input.txt read was done !'
  
    !二個目のプロットデータ読み込み !断面積分布読み込み
    Open(110,file="sv_input.txt")
    n2 = 0
    n = 0
    read (110, *) n
    n2 = int(n)
    allocate(x2(n2), y2(n2))

    do i = 1, n2
      read (110, *) x2(i), y2(i)
    end do
    Close(110)

    write(*,*) 'sv_input.txt read was done !'

!!!---rotate volume data---!!!    
    
    lr = x2(n2) - x2(1)
    write(*,*) 'lr = ', lr
    do i = 1, n2
      x2(i) = x2(i) * lv / lr - (x2(1) - x2(n2))
      y2(i) = y2(i) * 
    end do
  
  !一個目のプロットデータのx座標に合わせて
  !二個目のプロットデータのy座標を線形内挿・外挿
    Do i = 1, n1
      Do j = 1, n2
        If(x1(i) < x2(j)) Then
          If(j == 1) Then
            iy2(i) = OPf(x2(1), y2(1), x2(2), y2(2), x1(i))
            Exit
          Else
            iy2(i) = IP(x2(j-1), y2(j-1), x2(j), y2(j), x1(i))
            Exit
          End if
        End if
        If(j == n2) Then
          iy2(i) = OPl(x2(n2-1), y2(n2-1), x2(n2), y2(n2), x1(i))
        End if
      End do
    End do
    
    write(*,*) 'interpolation was done !'
  
  !二個目のプロットデータを一個目のプロットデータのx座標に合わせたプロットデータ
    Open(120,file="sv_interpolated.txt")
      Do i = 1, n1
        Write(120,*) x1(i), iy2(i)
      End do
    Close(120)

    write(*,*) 'Program Finish !'
  
End program Linear_Interpolation
  