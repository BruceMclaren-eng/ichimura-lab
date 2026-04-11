module fem_types
    implicit none
    type :: pair
        double precision :: x,y
    end type pair

end module fem_types

program fem_practice
    use fem_types
    implicit none


    !==================================
    !【Step0:配列の確保】
    !==================================
    !【メッシュ情報】
    integer :: nx, ny ! x,y方向の分割数
    integer :: n_nodes ! 全ノード数
    integer :: n_elems ! 全要素数
    double precision :: x_min, x_max, y_min, y_max !領域の定義
    integer :: n_edges !メッシュの辺の数
    
    ! 【動的配列】
    type(pair), allocatable :: nodes(:)!ノードの座標
    integer, allocatable :: elems(:,:) !要素のノード（要素の番号,要素上の接点番号）
    double precision, allocatable ::F_vec(:) !力ベクトル
    double precision, allocatable :: u_vec(:)!未知数ベクトル
    double precision, allocatable :: K_local(:,:)!要素剛性行列
    integer,          allocatable :: coo_row(:), coo_col(:)
    double precision, allocatable :: coo_val(:)
    integer,          allocatable :: row_ptr(:), col_idx(:)
    double precision, allocatable :: values(:)

    !　【静的配列】
    integer :: coo_nnz
    integer :: coo_ptr
    integer :: e, i, j
    integer :: global_nodes(3)
    logical, allocatable :: is_fixed(:)
    real :: t_start, t_end

    call cpu_time(t_start)
    !===================================
    !【Step1:メッシュの作成】
    !===================================    
    ! 領域の定義
    x_min = 0.0d0
    x_max = 1.0d0
    y_min = 0.0d0
    y_max = 2.0d0
    !分割数の定義
    nx = 500
    ny = 1000
    !全ノード数
    n_nodes = (nx + 1) * (ny + 1)
    !全要素数
    n_elems = nx * ny * 2
    !メッシュの辺の数
    n_edges = 3


    !==================================
    !【Step2:配列の初期化】
    !==================================
    allocate(nodes(n_nodes))
    allocate(elems(n_elems,n_edges))
    allocate(F_vec(n_nodes))
    allocate(u_vec(n_nodes))
    allocate(K_local(n_edges,n_edges))


    !==================================
    !【Step3:メッシュ生成】
    !==================================
    call generate_nodes(nx,ny,x_max,x_min,y_max,y_min,nodes)
    call generate_elems(nx,ny,elems)
    print *, "Mesh Generated"
    print *, "Generated Nodes : " , n_nodes
    print *, "Generated Elements" , n_elems


    !==================================
    !【Step4:要素剛性行列作成】
    !==================================
    coo_nnz = n_elems * 9 !非ゼロ要素の最大値で定義
    allocate(coo_row(coo_nnz))
    allocate(coo_col(coo_nnz))
    allocate(coo_val(coo_nnz))
    coo_ptr = 0

    do e = 1, n_elems
        global_nodes = elems(e,:)
        call calc_elem_stiffness(nodes(global_nodes(1)), nodes(global_nodes(2)), nodes(global_nodes(3)), K_local)
        do i = 1, 3
            do j = 1,3
                coo_ptr = coo_ptr + 1
                coo_row(coo_ptr) = global_nodes(i)
                coo_col(coo_ptr) = global_nodes(j)
                coo_val(coo_ptr) = K_local(i,j)
            end do
        end do
    end do

    print *, "COO Assembly Done"

    call coo_to_crs(coo_row, coo_col, coo_val, coo_ptr, n_nodes, row_ptr, col_idx, values)
    print *, "CRS Conversion Done: nnz =", row_ptr(n_nodes + 1) - 1
    deallocate(coo_row, coo_col, coo_val)

    !==================================
    !【Step5:ベクトルFの作成】
    !==================================
    allocate(is_fixed(n_nodes))
    is_fixed = .false.
    do i = 1, n_nodes
        if(abs(nodes(i)%y - 0.0d0) < 1.0d-6) then
            is_fixed(i) = .true.
            U_vec(i) = 0.0d0
            call dirichlet_row(i, n_nodes, values, col_idx, row_ptr, F_vec, 0.0d0)
        else if(abs(nodes(i)%y - 2.0d0) < 1.0d-6) then
            is_fixed(i) = .true.
            U_vec(i) = 1.0d0
            call dirichlet_row(i, n_nodes, values, col_idx, row_ptr, F_vec, 1.0d0)
        end if
    end do

    !==================================
    !【Step6:方程式を解く】
    !==================================
    
    call cg_solver(n_nodes, values, col_idx, row_ptr, F_vec, U_vec)

    print *, "========================================"
    print *, " FINAL RESULTS (Nodal Values)"
    print *, "========================================"
    print *, " Node ID |   X      Y    |   U(x,y)  | Exact(y/2)"
    print *, "---------+---------------+-----------+-----------"
    do i = 1, n_nodes
        write(*, '(I6, "   |", F6.2, " ", F6.2, "  |", F10.5, " |", F10.5)') &
            i, nodes(i)%x, nodes(i)%y, U_vec(i), nodes(i)%y / 2.0d0
    end do

    deallocate(nodes, elems, F_vec, U_vec, is_fixed)
    deallocate(row_ptr, col_idx, values)


    !==================================
    !【Step7:計算時間の出力】
    !==================================
    call cpu_time(t_end)
    print *, "Total Computation Time: ", t_end - t_start, " seconds"


contains
    !節点の座標と要素のノード番号の定義
    subroutine generate_nodes(nx,ny,x_max,x_min,y_max,y_min,nodes)
        implicit none
        integer, intent(in):: nx, ny !the number of node split
        double precision, intent(in) :: x_max, x_min, y_max, y_min ! maximam and minimam x, y
        type(pair), intent(out) :: nodes(:) !nodes' coordinate
        integer :: i, j, node_id
        double precision :: dx, dy

        
        dx = (x_max - x_min) / dble(nx)
        dy = (y_max - y_min) / dble(ny)

        node_id = 0
        do j = 0, ny
            do i = 0, nx
                node_id = node_id +1
                nodes(node_id)%x = x_min + dble(i) * dx
                nodes(node_id)%y = y_min + dble(j) * dy
            end do
        end do
    end subroutine

    !要素コネクティビティの生成
    subroutine generate_elems(nx,ny,elems)
        implicit none
        integer, intent(in) :: nx, ny
        integer, intent(out) :: elems(:,:)
        integer :: i,j, elem_id
        integer :: n1, n2, n3, n4

        elem_id = 0
        do j = 0, ny-1
            do i = 0, nx-1
                n1 = j * (nx + 1) + i + 1  !left lower node
                n2 = n1 + 1             !right lower node
                n3 = n1 + nx + 1        !left upper node
                n4 = n3 + 1             !right upper node

                !triangle 1: (n1, n2, n4)
                elem_id = elem_id + 1
                elems(elem_id,1) = n1
                elems(elem_id,2) = n2
                elems(elem_id,3) = n4

                !triangle 2: (n1,n4,n3)
                elem_id = elem_id + 1
                elems(elem_id,1) = n1
                elems(elem_id,2) = n4
                elems(elem_id,3) = n3
            end do
        end do
    end subroutine

    !要素剛性行列の計算
    subroutine calc_elem_stiffness(p1,p2,p3, k_mat)
        implicit none
        type(pair), intent(in) :: p1,p2,p3
        double precision, intent(out) :: k_mat(3,3)
        double precision :: x(3), y(3)
        double precision :: mat_J(2,2), mat_L(2,2)
        double precision :: detJ, val_m, area
        double precision :: v(3,2), res(3,2)
        integer :: i,j_idx

        !節点座標をセット
        x(1) = p1%x; y(1) = p1%y
        x(2) = p2%x; y(2) = p2%y
        x(3) = p3%x; y(3) = p3%y

        !ヤコビアンを作成
        mat_J(1,1) = - x(1) + x(2)
        mat_J(1,2) = - y(1) + y(2)
        mat_J(2,1) = - x(1) + x(3)
        mat_J(2,2) = - y(1) + y(3)


        detJ = mat_J(1,1) * mat_J(2,2) - mat_J(1,2) * mat_J(2,1) !行列式の計算
        val_m = 1.0d0 / detJ !逆行列の係数
        area = abs(detJ) * 0.5d0 !三角形の実面積
        
        !ヤコビアンの逆行列計算
        mat_L(1,1) = val_m * mat_J(2,2)
        mat_L(1,2) = -val_m * mat_J(1,2)
        mat_L(2,1) = -val_m * mat_J(2,1)
        mat_L(2,2) = val_m * mat_J(1,1)

        !局所微分の計算
        v(1,1) = - 1.0d0; v(1,2) = - 1.0d0
        v(2,1) =   1.0d0; v(2,2) =   0.0d0
        v(3,1) =   0.0d0; v(3,2) =   1.0d0

        !全体座標微分への変換
        do i = 1, 3
            res(i,1) = mat_J(1,1) * v(i,1) + mat_L(2,1) * v(i,2)
            res(i,2) = mat_J(2,1) * v(i,1) + mat_L(2,2) * v(i,2)
        end do

        !要素剛性行列 k_matの計算
        do i = 1, 3
            do j_idx = 1,3
                k_mat(i, j_idx) = (res(i,1)*res(j_idx,1) + res(i,2)*res(j_idx,2)) * area
            end do
        end do
    end subroutine

    !要素剛性行列COOから全体剛性行列CRSへの変換
    subroutine coo_to_crs(coo_row, coo_col, coo_val, nnz_in, n, row_ptr, col_idx, values)
        integer,                        intent(in) :: coo_row(:), coo_col(:), nnz_in, n 
        double precision,               intent(in) :: coo_val(:)
        integer, allocatable,          intent(out) :: row_ptr(:), col_idx(:)
        double precision, allocatable, intent(out) :: values(:)

        integer          :: i, k, col, nnz
        integer          :: used_cols(20)
        integer          :: used_count
        integer,          allocatable :: pos(:)
        integer,          allocatable :: sorted_col(:)
        integer,          allocatable :: row_ptr_orig(:)
        double precision, allocatable :: sorted_val(:)
        double precision, allocatable :: dense_buffer(:)


        ! row_ptrの作成
        allocate(row_ptr(n+1))
        row_ptr = 0

        do k = 1, nnz_in
            row_ptr(coo_row(k) + 1) = row_ptr(coo_row(k) + 1) + 1
        end do

        row_ptr(1) = 1
        do i = 2, n+1
            row_ptr(i) = row_ptr(i) + row_ptr(i-1)
        end do

        !cooの暫定準備
        allocate(pos(n))
        allocate(sorted_col(nnz_in))
        allocate(sorted_val(nnz_in))

        pos = row_ptr(1:n)

        do k = 1, nnz_in
            i = coo_row(k)
            sorted_col(pos(i)) = coo_col(k)
            sorted_val(pos(i)) = coo_val(k)
            pos(i) = pos(i) + 1
        end do
        deallocate(pos)

        !重複の処理
        allocate(row_ptr_orig(n+1))
        row_ptr_orig = row_ptr

        allocate(dense_buffer(n))
        allocate(col_idx(nnz_in))
        allocate(values(nnz_in))

        dense_buffer = 0.0d0

        nnz = 0
        row_ptr(1) = 1

        do i = 1, n
            used_count = 0

            do k = row_ptr_orig(i), row_ptr_orig(i+1) - 1
                col = sorted_col(k)

                if (dense_buffer(col) == 0.0d0) then
                    used_count = used_count + 1
                    used_cols(used_count) = col
                end if

                dense_buffer(col) = dense_buffer(col) + sorted_val(k)
            end do

            do k = 1, used_count
                col = used_cols(k)
                nnz = nnz + 1
                col_idx(nnz) = col
                values(nnz) = dense_buffer(col)
                dense_buffer(col) = 0.0d0
            end do 

            row_ptr(i+1) = nnz + 1
        end do

        deallocate(row_ptr_orig,dense_buffer,sorted_col,sorted_val)

    end subroutine

    subroutine matvec(n, values, col_idx, row_ptr, x, y)
        integer, intent(in) :: n
        double precision, intent(in) :: values(:), x(:)
        integer, intent(in) :: col_idx(:), row_ptr(:)
        double precision, intent(inout) :: y(:)

        integer :: i, k
        !$acc kernels
        do i = 1, n
            y(i) = 0.0d0
            do k = row_ptr(i), row_ptr(i+1) - 1
                y(i) = y(i) + values(k) * x(col_idx(k))
            end do
        end do
        !$acc end kernels
    end subroutine

    
    subroutine dirichlet_row(node, n, valuesm, col_idx, row_ptr, F_vec, bc_val)
        implicit none
        integer, intent(in) :: node, n
        double precision, intent(inout) :: valuesm(:)   
        integer, intent(in) :: col_idx(:), row_ptr(:)
        double precision, intent(inout) :: F_vec(:)     
        double precision, intent(in) :: bc_val
        integer :: k

        do k = row_ptr(node), row_ptr(node + 1) - 1
            if (col_idx(k) == node) then
                valuesm(k) = 1.0d0   
            else
                valuesm(k) = 0.0d0  
            end if
        end do

        F_vec(node) = bc_val
    end subroutine

    subroutine cg_solver(n, values, col_idx, row_ptr, b, x)
        integer, intent(in) :: n
        double precision, intent(in) :: values(:), b(:)
        integer, intent(in) :: col_idx(:), row_ptr(:)
        double precision, intent(inout) :: x(:)

        double precision :: r(n), p(n), Ap(n)
        double precision :: alpha, beta, rr, rr_new, pAp
        integer :: iter, ii, kk

        !====================================
        !【データ転送（CPU→GPU）】
        !====================================
        !$acc data copyin(values, col_idx, row_ptr, b) copy(x) create(r, p, Ap, alpha, beta, rr, rr_new, pAp)


        call matvec(n, values, col_idx, row_ptr, x, Ap)

        !【初期残差と初期探索方向の計算】
        !$acc parallel loop
        do i = 1, n
            r(i) = b(i) - Ap(i)
            p(i) = r(i)
        end do
        !$acc end parallel loop


        !【初期残差の大きさ計算】
        rr = 0.0d0
        !$acc parallel loop reduction(+:rr)
        do i = 1, n
            rr = rr + r(i) * r(i)
        end do
        !$acc end parallel loop


        do iter = 1,2 *n
            if (rr < 1.0d0-20) exit
                
            call matvec(n, values, col_idx, row_ptr, p, Ap)
            !【alphaの計算】
            pAp = 0.0d0
            !$acc parallel loop reduction(+:pAp)
            do i = 1, n
                pAp = pAp + p(i) * Ap(i)
            end do
            !$acc end parallel loop
            alpha = rr / pAp

            !【xとrの更新】
            !$acc parallel loop
            do i = 1, n
                x(i) = x(i) + alpha * p(i)
                r(i) = r(i) - alpha * Ap(i)
            end do
            !$acc end parallel loop

            !【新しい残差の大きさ計算】
            rr_new = 0.0d0
            !$acc parallel loop reduction(+:rr_new)
            do i = 1, n
                rr_new = rr_new + r(i) * r(i)
            end do
            !$acc end parallel loop

            !【betaの計算】
            beta = rr_new / rr

            !【pの更新】
            !$acc parallel loop
            do i = 1, n
                p(i) = r(i) + beta * p(i)
            end do
            !$acc end parallel loop

            !【残差の大きさの更新】
            rr = rr_new
        end do
        !$acc end data
    end subroutine

end program 