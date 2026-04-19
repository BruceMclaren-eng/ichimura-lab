module fem_types
    implicit none
    type :: pair
        double precision :: x,y
    end type pair
end module fem_types

program fem_gauss
    use fem_types
    implicit none

    !==================================
    !【Step0:配列の確保】
    !==================================
    integer :: nx, ny
    integer :: n_nodes
    integer :: n_elems
    double precision :: x_min, x_max, y_min, y_max
    integer :: n_edges

    type(pair), allocatable :: nodes(:)
    integer, allocatable :: elems(:,:)
    double precision, allocatable :: F_vec(:)
    double precision, allocatable :: u_vec(:)
    double precision, allocatable :: K_local(:,:)
    double precision, allocatable :: K_semilocal(:,:,:)
    double precision, allocatable :: K_dense(:,:)   ! 全体剛性行列（密行列）
    logical, allocatable :: is_fixed(:)

    integer :: e, i, j
    integer :: global_nodes(3)
    real :: t_start, t_end
    real :: t_gauss_start, t_gauss_end

    call cpu_time(t_start)

    !===================================
    !【Step1:メッシュの作成】
    !===================================
    x_min = 0.0d0
    x_max = 1.0d0
    y_min = 0.0d0
    y_max = 2.0d0
    nx = 50
    ny = 50
    n_nodes = (nx + 1) * (ny + 1)
    n_elems = nx * ny * 2
    n_edges = 3

    !==================================
    !【Step2:配列の初期化】
    !==================================
    allocate(nodes(n_nodes))
    allocate(elems(n_elems, n_edges))
    allocate(F_vec(n_nodes))
    allocate(u_vec(n_nodes))
    allocate(K_local(n_edges, n_edges))
    allocate(K_dense(n_nodes, n_nodes))
    allocate(is_fixed(n_nodes))

    !==================================
    !【Step3:メッシュ生成】
    !==================================
    call generate_nodes(nx, ny, x_max, x_min, y_max, y_min, nodes)
    call generate_elems(nx, ny, elems)
    print *, "Mesh Generated"
    print *, "Generated Nodes : ", n_nodes
    print *, "Generated Elements: ", n_elems

    !==================================
    !【Step4:全体剛性行列の組み立て（密行列）】
    !==================================
    allocate(K_semilocal(n_elems, n_edges, n_edges))
    K_dense = 0.0d0

    ! 各要素の剛性行列を計算してK_semilocalに保存
    do e = 1, n_elems
        global_nodes = elems(e, :)
        call calc_elem_stiffness(nodes(global_nodes(1)), nodes(global_nodes(2)), nodes(global_nodes(3)), K_local)
        K_semilocal(e, :, :) = K_local
    end do

    ! K_semilocalからグローバル剛性行列（密行列）にアセンブリ
    call assemble_global_dense(K_semilocal, elems, n_elems, n_edges, n_nodes, K_dense)
    print *, "Global Stiffness Matrix Assembled"

    !==================================
    !【Step5:ベクトルFの作成・境界条件の適用】
    !==================================
    F_vec = 0.0d0
    u_vec = 0.0d0
    is_fixed = .false.

    do i = 1, n_nodes
        if (abs(nodes(i)%y - 0.0d0) < 1.0d-6) then
            is_fixed(i) = .true.
            u_vec(i) = 0.0d0
            call dirichlet_dense(i, n_nodes, K_dense, F_vec, 0.0d0)
        else if (abs(nodes(i)%y - 2.0d0) < 1.0d-6) then
            is_fixed(i) = .true.
            u_vec(i) = 1.0d0
            call dirichlet_dense(i, n_nodes, K_dense, F_vec, 1.0d0)
        end if
    end do

    !==================================
    !【Step6:ガウスの消去法で解く】
    !==================================

    call cpu_time(t_gauss_start)
    call gauss_elimination(n_nodes, K_dense, F_vec, u_vec)
    call cpu_time(t_gauss_end)
    !print *, "Gauss Elimination Done"
    print *, "Gauss Elimination Time: ", t_gauss_end - t_gauss_start, " seconds"

    !==================================
    !【Step7:結果の出力】
    !==================================
    !print *, "========================================"
    !print *, " FINAL RESULTS (Nodal Values)"
    !print *, "========================================"
    !print *, " Node ID |   X      Y    |   U(x,y)  | Exact(y/2)"
    !print *, "---------+---------------+-----------+-----------"
    !do i = 1, n_nodes
     !   write(*, '(I6, "   |", F6.2, " ", F6.2, "  |", F10.5, " |", F10.5)') &
     !       i, nodes(i)%x, nodes(i)%y, u_vec(i), nodes(i)%y / 2.0d0
    !end do

    deallocate(nodes, elems, F_vec, u_vec, K_local, K_dense, K_semilocal, is_fixed)

    call cpu_time(t_end)
    !print *, "Total Computation Time: ", t_end - t_start, " seconds"

contains

    !----------------------------------
    ! 節点座標の生成
    !----------------------------------
    subroutine generate_nodes(nx, ny, x_max, x_min, y_max, y_min, nodes)
        implicit none
        integer, intent(in) :: nx, ny
        double precision, intent(in) :: x_max, x_min, y_max, y_min
        type(pair), intent(out) :: nodes(:)
        integer :: i, j, node_id
        double precision :: dx, dy

        dx = (x_max - x_min) / dble(nx)
        dy = (y_max - y_min) / dble(ny)
        node_id = 0
        do j = 0, ny
            do i = 0, nx
                node_id = node_id + 1
                nodes(node_id)%x = x_min + dble(i) * dx
                nodes(node_id)%y = y_min + dble(j) * dy
            end do
        end do
    end subroutine

    !----------------------------------
    ! 要素コネクティビティの生成
    !----------------------------------
    subroutine generate_elems(nx, ny, elems)
        implicit none
        integer, intent(in) :: nx, ny
        integer, intent(out) :: elems(:,:)
        integer :: i, j, elem_id
        integer :: n1, n2, n3, n4

        elem_id = 0
        do j = 0, ny - 1
            do i = 0, nx - 1
                n1 = j * (nx + 1) + i + 1
                n2 = n1 + 1
                n3 = n1 + nx + 1
                n4 = n3 + 1

                elem_id = elem_id + 1
                elems(elem_id, 1) = n1
                elems(elem_id, 2) = n2
                elems(elem_id, 3) = n4

                elem_id = elem_id + 1
                elems(elem_id, 1) = n1
                elems(elem_id, 2) = n4
                elems(elem_id, 3) = n3
            end do
        end do
    end subroutine

    !----------------------------------
    ! 要素剛性行列の計算
    !----------------------------------
    subroutine calc_elem_stiffness(p1, p2, p3, k_mat)
        implicit none
        type(pair), intent(in) :: p1, p2, p3
        double precision, intent(out) :: k_mat(3,3)
        double precision :: x(3), y(3)
        double precision :: mat_J(2,2), mat_L(2,2)
        double precision :: detJ, val_m, area
        double precision :: v(3,2), res(3,2)
        integer :: i, j

        x(1) = p1%x; y(1) = p1%y
        x(2) = p2%x; y(2) = p2%y
        x(3) = p3%x; y(3) = p3%y

        mat_J(1,1) = -x(1) + x(2)
        mat_J(1,2) = -y(1) + y(2)
        mat_J(2,1) = -x(1) + x(3)
        mat_J(2,2) = -y(1) + y(3)

        detJ  = mat_J(1,1) * mat_J(2,2) - mat_J(1,2) * mat_J(2,1)
        val_m = 1.0d0 / detJ
        area  = abs(detJ) * 0.5d0

        mat_L(1,1) =  val_m * mat_J(2,2)
        mat_L(1,2) = -val_m * mat_J(1,2)
        mat_L(2,1) = -val_m * mat_J(2,1)
        mat_L(2,2) =  val_m * mat_J(1,1)

        v(1,1) = -1.0d0; v(1,2) = -1.0d0
        v(2,1) =  1.0d0; v(2,2) =  0.0d0
        v(3,1) =  0.0d0; v(3,2) =  1.0d0

        do i = 1, 3
            res(i,1) = mat_J(1,1) * v(i,1) + mat_L(2,1) * v(i,2)
            res(i,2) = mat_J(1,2) * v(i,1) + mat_L(2,2) * v(i,2)
        end do

        do j = 1, 3
            do i = 1, 3
                k_mat(i,j) = (res(i,1)*res(j,1) + res(i,2)*res(j,2)) * area
            end do
        end do
    end subroutine

    !----------------------------------
    ! 全体剛性行列へのアセンブリ（密行列版）
    ! K_local(i,j) → K_dense(global_nodes(i), global_nodes(j)) に足し込む
    !----------------------------------
    subroutine assemble_global_dense(K_semilocal, elems, n_elems, n_edges, n_nodes, K_dense)
        implicit none
        integer, intent(in) :: n_elems, n_edges, n_nodes
        double precision, intent(in) :: K_semilocal(:,:,:)
        integer, intent(in) :: elems(:,:)
        double precision, intent(inout) :: K_dense(:,:)
        integer :: e, i, j, gi, gj

        do e = 1, n_elems
            do j = 1, n_edges
                gj = elems(e, j)       ! 局所j → グローバル節点番号
                do i = 1, n_edges
                    gi = elems(e, i)   ! 局所i → グローバル節点番号
                    K_dense(gi, gj) = K_dense(gi, gj) + K_semilocal(e, i, j)
                end do
            end do
        end do
    end subroutine

    !----------------------------------
    ! ディリクレ境界条件の適用（密行列版）
    ! 該当行を単位行に変換し、F_vec(node) = bc_val に設定
    !----------------------------------
    subroutine dirichlet_dense(node, n, K_dense, F_vec, bc_val)
        implicit none
        integer, intent(in) :: node, n
        double precision, intent(inout) :: K_dense(:,:)
        double precision, intent(inout) :: F_vec(:)
        double precision, intent(in) :: bc_val
        integer :: j

        do j = 1, n
            if (j == node) then
                K_dense(node, j) = 1.0d0
            else
                K_dense(node, j) = 0.0d0
            end if
        end do
        F_vec(node) = bc_val
    end subroutine

    !----------------------------------
    ! ガウスの消去法（部分ピボット選択付き）
    ! Ax = b を解く
    !----------------------------------
    subroutine gauss_elimination(n, A, b, x)
        implicit none
        integer, intent(in) :: n
        double precision, intent(inout) :: A(n,n), b(n)
        double precision, intent(out) :: x(n)
        double precision :: factor, tmp
        integer :: i, j, k, pivot_row
        double precision :: max_val

        ! 前進消去（部分ピボット選択）
        do k = 1, n - 1

            ! ピボット選択：k列でk行以下の最大要素を探す
            max_val = abs(A(k,k))
            pivot_row = k
            do i = k + 1, n
                if (abs(A(i,k)) > max_val) then
                    max_val = abs(A(i,k))
                    pivot_row = i
                end if
            end do

            ! ピボット行と現在行を交換
            if (pivot_row /= k) then
                do j = 1, n
                    tmp = A(k,j)
                    A(k,j) = A(pivot_row,j)
                    A(pivot_row,j) = tmp
                end do
                tmp = b(k)
                b(k) = b(pivot_row)
                b(pivot_row) = tmp
            end if

            if (abs(A(k,k)) < 1.0d-15) then
                print *, "ERROR: Zero pivot at row", k
                stop
            end if

            ! 消去
            do i = k + 1, n
                factor = A(i,k) / A(k,k)
                do j = k, n
                    A(i,j) = A(i,j) - factor * A(k,j)
                end do
                b(i) = b(i) - factor * b(k)
            end do
        end do

        ! 後退代入
        x(n) = b(n) / A(n,n)
        do i = n - 1, 1, -1
            x(i) = b(i)
            do j = i + 1, n
                x(i) = x(i) - A(i,j) * x(j)
            end do
            x(i) = x(i) / A(i,i)
        end do
    end subroutine

end program