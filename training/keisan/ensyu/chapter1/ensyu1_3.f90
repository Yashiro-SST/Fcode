program ensyu1_3
 implicit none
 integer a(10),i
 a(1)=1
 a(2)=2
 do i=3,10
  a(i)=a(i-1)+a(i-2)
 enddo
 do i=1,10
  write(*,*) a(i)
 end do
end program ensyu1_3