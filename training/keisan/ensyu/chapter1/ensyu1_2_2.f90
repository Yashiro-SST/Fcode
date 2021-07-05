program ensyu1_2_2
 implicit none
 integer i,k,k2,k3
 k=0
 k2=0
 k3=0
 do i=1,20
  k = i*(i+1)/2
  k2= i*(i+1)*(2*i+1)/6
  k3= ((i**2)*((i+1)**2))/4
 enddo
 write(*,*)'k, k2, k3= ', k, k2, k3
end program ensyu1_2_2