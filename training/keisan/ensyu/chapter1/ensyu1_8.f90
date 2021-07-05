program ensyu1_8
    implicit none
    integer m,n,i,p
do
    write(*,*)'input n(natural number,<1,000,000):'
    read(*,*) n
    if(n<2) then
        write(*,*)'Its 1,1 is not prime number...'
    else if(n==2) then
        write(*,*)'n is prime number!'
    else
        do i=2,int(sqrt(dble(n)))
            if(mod(n,i)/=0) then
                continue
            else
                stop 'n is not prime number...'
            end if
        enddo
        write(*,*)'n is prime number!!'
    end if
enddo
end program ensyu1_8