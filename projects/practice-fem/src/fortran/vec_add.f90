program vec_add
    implicit none
    integer, parameter :: N = 10000000  ! 1000万要素
    real(8) :: a(N), b(N), c(N)
    integer :: i
    real :: t_start, t_end

    ! 初期化
    do i = 1, N
        a(i) = 1.0d0
        b(i) = 2.0d0
    end do

    ! CPU計算
    call cpu_time(t_start)
    do i = 1, N
        c(i) = a(i) + b(i)
    end do
    call cpu_time(t_end)
    print *, "CPU時間：", t_end - t_start, "秒"

    ! GPU計算
    call cpu_time(t_start)
    !$acc parallel loop
    do i = 1, N
        c(i) = a(i) + b(i)
    end do
    !$acc end parallel loop
    call cpu_time(t_end)
    print *, "GPU時間：", t_end - t_start, "秒"

    print *, "c(1) =", c(1)  ! 答えの確認（3.0になるはず）

end program vec_add