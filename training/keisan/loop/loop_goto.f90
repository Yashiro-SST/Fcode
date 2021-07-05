program loop_goto
 implicit none
 integer wa,n,i
 do 
    write(*,'(a\)') 'input n (input 0 to stop) :'
    read (*,*) n
    if (n == 0) then
        goto 1
    else if (n < 0) then
        write(*,*) 'sorry, input positive n ...'
    else
        wa=0
        do i=1,n
            wa = wa + i
        enddo
        write(*,*) 'wa = ', wa
    end if
 enddo
1 continue
 write(*,*) 'exit from do-loop ...'
end program loop_goto