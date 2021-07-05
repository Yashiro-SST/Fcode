program mat3d
  implicit none
  integer i, j, k, is, n
  real(8) t1, t2
  integer, allocatable :: a(:, :, :)
  write(*, '(a)', advance = 'no') 'input n : '
  read (*,*) n
  allocate (a(n,n,n), stat = is)
  if (is /= 0) stop 'cannot allocate (n is too large)' !割付け不可なら停止
  ! --- メモリ上で近接する要素に順に0を代入
  call cpu_time(t1)
  do k = 1, n
    do j = 1, n
      do i = 1, n
        a(i, j, k) = 0
      enddo
    enddo
  enddo
  call cpu_time(t2)
  write(*,*) 'cpu time = ', t2 -t1
  ! --- メモリ上で離れた位置にある要素に順に0を代入
  call cpu_time(t1)
  do i = 1, n
    do j = 1, n
      do k = 1, n
        a(i, j, k) = 0
      enddo
    enddo
  enddo
  call cpu_time(t2)
  write(*,*) 'cpu time = ', t2 - t1
  deallocate(a)
end program mat3d