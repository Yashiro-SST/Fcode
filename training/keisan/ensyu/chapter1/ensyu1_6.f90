program ensyu1_6
 implicit none
 integer wa,m,n,i
 do 
    write(*,*) 'input m (if input negative m,stop!) :'
    read (*,*) m
    if (n < 0) stop 'good bye ...'
    write(*,*) 'input n :'
    read (*,*) n
    wa=0
    if (m <= n) then
        do i=m,n
            wa = wa + i
        enddo
    else if (m > n) then
        do i=n,m
            wa = wa + i
        enddo
    else
        stop 'something is wrong !'
    endif
    write(*,*) 'wa = ', wa
 enddo
end program ensyu1_6