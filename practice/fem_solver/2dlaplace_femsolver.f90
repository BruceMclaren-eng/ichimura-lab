module fem_types
    implicit none
    type :: pair
        real(4) :: x,y
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
    real(4) :: x_min, x_max, y_min, y_max !領域の定義
    integer :: n_edges !メッシュの辺の数
    
    ! 【動的配列】
    type(pair), allocatable :: nodes(:)!ノードの座標
    integer, allocatable :: elems(:,:) !要素のノード（要素の番号,要素上の接点番号）
    real(4), allocatable ::F_vec(:) !力ベクトル
    real(4), allocatable :: u_vec(:)!未知数ベクトル
    real(4), allocatable :: K_local(:,:)!要素剛性行列
    integer,          allocatable :: row_ptr(:), col_idx(:)
    real(4), allocatable :: values(:)
    real(4), allocatable :: diag(:)
    real(4), allocatable :: K_semilocal(:, :, :)
    logical, allocatable :: is_fixed(:)

    !　【静的配列】
    integer :: e, i, j, k
    integer :: global_nodes(3)
    real :: t_start, t_end
    integer(8) :: t1, t2, rate, a1, a2, b1, b2, c1, c2, d1, d2, e1, e2, f1, f2
    integer :: nnz_global
    real :: t_pcg_start, t_pcg_end
    
    call system_clock(t1, rate)
    !===================================
    !【Step1:メッシュの作成】
    !===================================    
    ! 領域の定義
    x_min = 0.0e0
    x_max = 1.0e0
    y_min = 0.0e0
    y_max = 2.0e0
    !分割数の定義
    nx = 4000
    ny = 4000
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
    call system_clock(a1, rate)
    call generate_nodes(nx,ny,x_max,x_min,y_max,y_min,nodes)
    call generate_elems(nx,ny,elems)
    call system_clock(a2, rate)
    print *, "Mesh Generated"
    print *, "Generated Nodes : " , n_nodes
    print *, "Generated Elements" , n_elems


    !==================================
    !【Step4:要素剛性行列作成】
    !==================================
    allocate(K_semilocal(n_elems, n_edges, n_edges))

    call system_clock(c1, rate)
    do e = 1, n_elems
        global_nodes = elems(e,:)
        call calc_elem_stiffness(nodes(global_nodes(1)), nodes(global_nodes(2)), nodes(global_nodes(3)), K_local)
        K_semilocal(e,:,:) = K_local
    end do
    call system_clock(c2, rate)

    call system_clock(d1, rate)
    call calc_global_stiffness(K_semilocal, n_elems, n_edges, row_ptr, col_idx, values)
    call system_clock(d2, rate)
    print *, "CRS Conversion Done: nnz =", nnz_global

    !==================================
    !【Step5:ベクトルFの作成】
    !==================================
    allocate(is_fixed(n_nodes))
    is_fixed = .false.
    call system_clock(e1, rate)
    do i = 1, n_nodes
        if(abs(nodes(i)%y - 0.0e0) < 1.0e-6) then
            is_fixed(i) = .true.
            U_vec(i) = 0.0e0
            call dirichlet_row(i, n_nodes, values, col_idx, row_ptr, F_vec, 0.0e0)
        else if(abs(nodes(i)%y - 2.0e0) < 1.0d-6) then
            is_fixed(i) = .true.
            U_vec(i) = 1.0e0
            call dirichlet_row(i, n_nodes, values, col_idx, row_ptr, F_vec, 1.0e0)
        end if
    end do
    call system_clock(e2, rate)

    allocate(diag(n_nodes))
    do i = 1, n_nodes
        diag(i) = 0.0e0
        do j = row_ptr(i), row_ptr(i+1) - 1
            if (col_idx(j) == i) then
                diag(i) = values(j)
                exit
            end if
        end do
    end do
    !==================================
    !【Step6:方程式を解く】
    !==================================
    call system_clock(f1, rate)
    call pcg_solver(n_nodes, values, col_idx, row_ptr, F_vec, U_vec, diag)
    call system_clock(f2, rate)
    deallocate(diag)
    !print *, "Calculation Finished"
    print *, "CG Method Wall Time", real(f2-f1, 8) / real(rate, 8), "seconds"
    print *, "generation of nodes and mesh time: ", real(a2-a1, 8) / real(rate, 8), "seconds"
    print *, "Calculation of Element Stiffness Matrices time: ", real(c2-c1, 8) / real(rate, 8), "seconds"
    print *, "Calculation of Global Stiffness Matrix time: ", real(d2-d1, 8) / real(rate, 8), "seconds"
    print *, "Application of Boundary Conditions time: ", real(e2-e1, 8) / real(rate, 8), "seconds"
    print *, "Solution of Matrix Equation time: ", real(f2-f1, 8) / real(rate, 8), "seconds"
    !print *, "========================================"
    !print *, " FINAL RESULTS (Nodal Values)"
    !print *, "========================================"
    !print *, " Node ID |   X      Y    |   U(x,y)  | Exact(y/2)"
    !print *, "---------+---------------+-----------+-----------"
    !do i = 1, n_nodes
    !    write(*, '(I6, "   |", F6.2, " ", F6.2, "  |", F10.5, " |", F10.5)') &
    !        i, nodes(i)%x, nodes(i)%y, U_vec(i), nodes(i)%y / 2.0d0
    !end do

    !deallocate(nodes, elems, F_vec, U_vec, is_fixed)
    !deallocate(row_ptr, col_idx, values)


    !==================================
    !【Step7:計算時間の出力】
    !==================================
    call system_clock(t2, rate)
    print *, "Total Time: ", real(t2 - t1, 8) / real(rate, 8), " seconds"


contains
    !節点の座標と要素のノード番号の定義
    subroutine generate_nodes(nx,ny,x_max,x_min,y_max,y_min,nodes)
        implicit none
        integer, intent(in):: nx, ny !the number of node split
        real(4), intent(in) :: x_max, x_min, y_max, y_min ! maximam and minimam x, y
        type(pair), intent(out) :: nodes(:) !nodes' coordinate
        integer :: i, j, node_id
        real(4) :: dx, dy

        
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
        real(4), intent(out) :: k_mat(3,3)
        real(4) :: x(3), y(3)
        real(4) :: mat_J(2,2), mat_L(2,2)
        real(4) :: detJ, val_m, area
        real(4) :: v(3,2), res(3,2)
        integer :: i,j

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
        val_m = 1.0e0 / detJ !逆行列の係数
        area = abs(detJ) * 0.5e0 !三角形の実面積
        
        !ヤコビアンの逆行列計算
        mat_L(1,1) = val_m * mat_J(2,2)
        mat_L(1,2) = -val_m * mat_J(1,2)
        mat_L(2,1) = -val_m * mat_J(2,1)
        mat_L(2,2) = val_m * mat_J(1,1)

        !局所微分の計算
        v(1,1) = - 1.0e0; v(1,2) = - 1.0e0
        v(2,1) =   1.0e0; v(2,2) =   0.0e0
        v(3,1) =   0.0e0; v(3,2) =   1.0e0

        !全体座標微分への変換
        do i = 1, 3
            res(i,1) = mat_J(1,1) * v(i,1) + mat_L(2,1) * v(i,2)
            res(i,2) = mat_J(1,2) * v(i,1) + mat_L(2,2) * v(i,2)
        end do

        !要素剛性行列 k_matの計算
        do i = 1, 3
            do j = 1,3
                k_mat(i, j) = (res(i,1)*res(j,1) + res(i,2)*res(j,2)) * area
            end do
        end do
    end subroutine

    subroutine calc_global_stiffness(K_semilocal, n_elems, n_edges, row_ptr, col_idx, values)
        integer,                        intent(in) :: n_edges, n_elems
        real(4),               intent(in) :: K_semilocal(:, :, :)
        integer, allocatable,          intent(out) :: row_ptr(:), col_idx(:)
        real(4), allocatable, intent(out) :: values(:)

        integer  :: i, j, k, gi, gj, ptr, crs_pos, max_nnz
        integer, allocatable ::col_idx_temp(:,:), nnz_row(:)
        logical :: already_exists

        max_nnz = 10
        allocate(nnz_row(n_nodes))
        allocate(col_idx_temp(n_nodes, max_nnz))
        nnz_row = 0
        col_idx_temp = 0

        do e = 1, n_elems
            do i = 1, n_edges
                gi = elems(e, i)
                do j = 1, n_edges
                    gj = elems(e, j)
                    already_exists = .false.
                    do k = 1, nnz_row(gi)
                        if (col_idx_temp(gi, k) == gj) then 
                            already_exists = .true.
                            exit
                        end if
                    end do
                    if (.not. already_exists) then
                        nnz_row(gi) = nnz_row(gi) + 1
                        if (nnz_row(gi) > max_nnz) then
                            print *, "Error: Exceeded maximum non-zeros per row"
                            stop
                        end if
                    col_idx_temp(gi, nnz_row(gi)) = gj
                    end if
                end do
            end do
        end do

        !row_ptrの作成
        allocate(row_ptr(n_nodes + 1))
        row_ptr(1) = 1
        do gi = 1, n_nodes
            row_ptr(gi + 1) = row_ptr(gi) + nnz_row(gi)
        end do
        nnz_global = row_ptr(n_nodes + 1) - 1


        !col_idxの構築
        allocate(col_idx(nnz_global))
        crs_pos = 0
        do gi = 1, n_nodes
            do ptr = 1, nnz_row(gi)
                if (col_idx_temp(gi, ptr) /= 0) then
                    crs_pos = crs_pos + 1
                    col_idx(crs_pos) = col_idx_temp(gi, ptr)
                end if
            end do
        end do
        deallocate(nnz_row, col_idx_temp)

        !valuesの構築
        allocate(values(nnz_global))
        values = 0.0e0
        do e = 1, n_elems
            do i = 1, n_edges
                gi = elems(e, i)
                do j = 1, n_edges
                    gj = elems(e, j)
                    do k = row_ptr(gi), row_ptr(gi + 1) - 1
                        if (col_idx(k) == gj) then
                            values(k) = values(k) + K_semilocal(e, i, j)
                            exit
                        end if
                    end do
                end do
            end do
        end do

    end subroutine


    subroutine dirichlet_row(node, n, valuesm, col_idx, row_ptr, F_vec, bc_val)
        implicit none
        integer, intent(in) :: node, n
        real(4), intent(inout) :: valuesm(:)   
        integer, intent(in) :: col_idx(:), row_ptr(:)
        real(4), intent(inout) :: F_vec(:)     
        real(4), intent(in) :: bc_val
        integer :: k

        do k = row_ptr(node), row_ptr(node + 1) - 1
            if (col_idx(k) == node) then
                valuesm(k) = 1.0e0   
            else
                valuesm(k) = 0.0e0  
            end if
        end do

        F_vec(node) = bc_val
    end subroutine

    subroutine matvec(n, values, col_idx, row_ptr, x, Ap)
        integer, intent(in) :: n
        real(4), intent(in) :: values(:), x(:)
        integer, intent(in) :: col_idx(:), row_ptr(:)
        real(4), intent(inout) :: Ap(:)
        integer :: i, k

        !$acc parallel loop present(values, col_idx, row_ptr, x, Ap)
        do i = 1, n
            Ap(i) = 0.0e0
            do k = row_ptr(i), row_ptr(i+1) - 1
                Ap(i) = Ap(i) + values(k) * x(col_idx(k))
            end do
        end do
        !$acc end parallel loop
        !$acc wait
    end subroutine

   
    subroutine pcg_solver(n, values, col_idx, row_ptr, b, x, diag)
    integer, intent(in) :: n
    real(4), intent(in) :: values(:), b(:), diag(:)
    integer, intent(in) :: col_idx(:), row_ptr(:)
    real(4), intent(inout) :: x(:)

    real(4) :: r(n), p(n), Ap(n), z(n),r0_norm, r_norm_sq, r0_sq
    real(4) :: alpha, beta, rz, rz_new, pAp
    integer :: iter, i
    integer(8) :: a1,a2,b1,b2,c1,c2, d1, d2, rate
    real(8) :: t_spmv, t_pap, t_update, t_pvec

    t_spmv = 0.0d0
    t_pap = 0.0d0
    t_update = 0.0d0
    t_pvec = 0.0d0

    !$acc data copyin(values, col_idx, row_ptr, b, diag) &
    !$acc      copy(x) &
    !$acc      create(r, p, Ap, z)

    ! 初期残差
    call matvec(n, values, col_idx, row_ptr, x, Ap)
    
    r0_sq = 0.0e0
    !$acc parallel loop reduction(+:r0_sq)
    do i = 1, n
        r(i) = b(i) - Ap(i)
        z(i) = r(i) / diag(i)   ! Jacobi前処理
        p(i) = z(i)
        r0_sq = r0_sq + r(i) * r(i)
    end do
    !$acc end parallel loop

    rz = 0.0e0
    !$acc parallel loop reduction(+:rz)
    do i = 1, n
        rz = rz + r(i) * z(i)
    end do
    !$acc end parallel loop

    do iter = 1, 10

        call system_clock(a1, rate)
        call matvec(n, values, col_idx, row_ptr, p, Ap)
        call system_clock(a2, rate)

        t_spmv = t_spmv + real(a2 - a1, 8) / real(rate, 8)

        call system_clock(b1, rate)
        pAp = 0.0e0
        !$acc parallel loop reduction(+:pAp)
        do i = 1, n
            pAp = pAp + p(i) * Ap(i)
        end do
        !$acc end parallel loop
        !$acc wait
        call system_clock(b2, rate)

        t_pap = t_pap + real(b2 - b1, 8) / real(rate, 8)

        alpha = rz / pAp

        call system_clock(c1, rate)
        rz_new = 0.0e0
        r_norm_sq = 0.0e0
        !$acc parallel loop reduction(+:rz_new)
        do i = 1, n
            x(i) = x(i) + alpha * p(i)
            r(i) = r(i) - alpha * Ap(i)
            z(i) = r(i) / diag(i)   ! Jacobi前処理
            rz_new = rz_new + r(i) * z(i)
            r_norm_sq = r_norm_sq + r(i) * r(i)
        end do
        !$acc end parallel loop
        !$acc wait
        call system_clock(c2, rate)

        t_update = t_update + real(c2-c1, 8) / real(rate, 8)

        beta = rz_new / rz

        call system_clock(d1, rate)
        !$acc parallel loop
        do i = 1, n
            p(i) = z(i) + beta * p(i)
        end do
        !$acc end parallel loop
        !$acc wait
        call system_clock(d2, rate)

        

        t_pvec = t_pvec + real(d2-d1, 8) / real(rate, 8)

        rz = rz_new

        if (sqrt(r_norm_sq) / sqrt(r0_sq) < 1.0e-8) exit

    end do

    !$acc end data

    print *, "Number of iterations:", iter - 1
    print *, "t_spmv : ", t_spmv, "seconds"
    print *, "t_pap : ", t_pap, "seconds"
    print *, "t_update : ", t_update, "seconds"
    print *, "t_pvec : ", t_pvec, "seconds"
end subroutine

end program 