program newton_sqrt
    implicit none
    real(8) ::x1,x2,a,er,er0 = 1.0d-6        !er0は誤差の許容値(しきい値)
    integer ::k,km=100                       !kmは最大反復回数
    write(*,'(a)') 'input a : '
    read(*,*) a
    if(a <= 0.0d0) stop 'a <= 0.0d0'
    x1 = a                                   !近似解の初期値をaとする
    do k=1,km
        x2 = x1 - 0.5d0 * (x1 ** 2 -a) / x1  !新しい近似解x2を計算する
        er = abs(x2-x1)                      !x1とx2の差の絶対値をerとする
        if (er < er0) exit
        x1 = x2
    enddo
    write(*,*) 'kai, k, er = ', x2, k, er
end program newton_sqrt