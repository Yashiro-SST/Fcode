program ensyu1_2_1
 implicit none
 integer i,k,k2,k3
 k=0
 k2=0
 k3=0
 do i=1,20
  k = k+i
  k2=k2+i**2
  k3=k3+i**3
 enddo
 write(*,*)'k, k2, k3= ', k, k2, k3
end program ensyu1_2_1