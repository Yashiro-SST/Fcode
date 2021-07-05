program ensyu1_7
 implicit none
 integer P,C,r,n,i,sa,nf,rf,saf
 do 
    write(*,*) 'input n (please input n (0<n<=10) ):'
    read (*,*) n
    write(*,*) 'input r (please input r (0<r<n)   ):'
    read (*,*) r
    !
    if (n<=0.or.10<n) then
        write (*,*)'sorry,input n which meet 0<n<=10'
        exit
    else if(r<0.or.n<r) then
        write (*,*)'sorry,input r which meet 0<r<n'
        exit
    end if
    !
    sa=n-r
    nf=1
    rf=1
    saf=1
    do i=1,n
        nf=nf*i
    enddo
    do i=1,r
        rf=rf*i
    enddo
    do i=1,sa
        saf=saf*i
    enddo
    P=nf/saf
    C=nf/(rf*saf)
    write(*,*) 'P, C = ', P, C
 enddo
 write(*,*) 'exit from do-loop...'
end program ensyu1_7