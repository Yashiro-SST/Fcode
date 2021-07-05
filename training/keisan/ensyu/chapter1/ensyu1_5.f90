program ensyu1_5
 implicit none
 integer wa,m,n,i
 do 
    write(*,*) 'input m (if m>n,stop!) :'
    read (*,*) m
    write(*,*) 'input n (if m>n,stop!) :'
    read (*,*) n
    if (m > n) stop 'good bye ...'
    wa=0
    do i=m,n
        wa = wa + i
    enddo
    write(*,*) 'wa = ', wa
 enddo
end program ensyu1_5