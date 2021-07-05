program loop_exit
 implicit none
 integer wa,n,i
 do 
    write(*,'(a)') 'input n (input 0 to stop) :'
    read (*,*) n
    if (n == 0) then
        exit
    else if (n < 0) then
        write(*,*) 'sorry, input positive n ...'
        cycle
    end if
    wa=0
    do i=1,n
        wa = wa + i
    enddo
    write(*,*) 'wa = ', wa
 enddo
 write(*,*) 'exit from do-loop ...'
end program loop_exit