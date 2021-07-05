program loop_err
 implicit none !宣言をしない
 integer i,wa
 wa=0
 do i=1,100
  wa=va+i
 enddo
 write(*,*)'wa=',wa
end program loop_err