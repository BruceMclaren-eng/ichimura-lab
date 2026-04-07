!=============================================================================
! ラプラス方程式の有限要素法 (FEM) 解法  ── CRS版
!=============================================================================

module fem_types
    implicit none
    type :: pair
        double precision :: x, y
    end type pair
end module fem_types


program solve_laplace_fem
    use fem_types
    implicit none

    integer :: nx, ny
    integer :: n_nodes, n_elems
    integer :: i, j, k, e
    integer :: global_nodes(3)

    double precision :: x_min, x_max, y_min, y_max

    type(pair), allocatable :: nodes(:)
    integer,    allocatable :: elems(:,:)
    double precision, allocatable :: F_vec(:), U_vec(:)
    logical,          allocatable :: is_fixed(:)
    double precision :: K_local(3,3)
    double precision :: bc_value

    integer,          allocatable :: coo_row(:), coo_col(:)
    double precision, allocatable :: coo_val(:)
    integer :: coo_nnz
    integer :: coo_ptr

    integer,          allocatable :: row_ptr(:), col_idx(:)
    double precision, allocatable :: values(:)

    x_min = 0.0d0;  x_max = 1.0d0
    y_min = 0.0d0;  y_max = 2.0d0
    nx = 100
    ny = 200

    n_nodes = (nx + 1) * (ny + 1)
    n_elems = nx * ny * 2

    allocate(nodes(n_nodes))
    allocate(elems(n_elems, 3))
    allocate(F_vec(n_nodes))
    allocate(U_vec(n_nodes))
    allocate(is_fixed(n_nodes))

    F_vec    = 0.0d0
    U_vec    = 0.0d0
    is_fixed = .false.

    coo_nnz = n_elems * 9
    allocate(coo_row(coo_nnz))
    allocate(coo_col(coo_nnz))
    allocate(coo_val(coo_nnz))
    coo_ptr = 0

    call generate_mesh(nx, ny, x_min, x_max, y_min, y_max, nodes, elems)
    print *, "--- Mesh Generated: Nodes =", n_nodes, " Elements =", n_elems

    do e = 1, n_elems
        global_nodes = elems(e, :)
        call calc_elem_stiffness( &
            nodes(global_nodes(1)), nodes(global_nodes(2)), nodes(global_nodes(3)), K_local)
        do i = 1, 3
            do j = 1, 3
                coo_ptr = coo_ptr + 1
                coo_row(coo_ptr) = global_nodes(i)
                coo_col(coo_ptr) = global_nodes(j)
                coo_val(coo_ptr) = K_local(i, j)
            end do
        end do
    end do
    print *, "--- COO Assembly Done: entries (with duplicates) =", coo_ptr

    call coo_to_crs(coo_row, coo_col, coo_val, coo_ptr, n_nodes, row_ptr, col_idx, values)
    print *, "--- CRS Conversion Done: nnz (unique) =", row_ptr(n_nodes+1)-1

    deallocate(coo_row, coo_col, coo_val)
    print *, "--- COO deallocated"

    print *, "--- Applying Boundary Conditions ---"

    do i = 1, n_nodes
        bc_value = 0.0d0
        if (abs(nodes(i)%y - 0.0d0) < 1.0d-6) then
            is_fixed(i) = .true.
            bc_value    = 0.0d0
        else if (abs(nodes(i)%y - 2.0d0) < 1.0d-6) then
            is_fixed(i) = .true.
            bc_value    = 1.0d0
        end if
        if (is_fixed(i)) then
            U_vec(i) = bc_value
            do k = 1, n_nodes
                call crs_subtract_col(k, i, bc_value, &
                                      values, col_idx, row_ptr, F_vec)
            end do
        end if
    end do

    do i = 1, n_nodes
        if (is_fixed(i)) then
            call apply_dirichlet_row(i, n_nodes, values, col_idx, row_ptr, F_vec, U_vec(i))
        end if
    end do

    print *, "--- Solving (CG method) ---"
    U_vec = 0.0d0
    do i = 1, n_nodes
        if (is_fixed(i)) then
            if (abs(nodes(i)%y - 2.0d0) < 1.0d-6) U_vec(i) = 1.0d0
        end if
    end do

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


contains

    !==========================================================================
    ! subroutine: generate_mesh
    !
    ! 【目的】
    !   矩形領域を nx×ny の格子に分割し、各矩形を対角線（↗方向）で
    !   2つの三角形要素に分けて nodes と elems を生成する。
    !
    ! 【節点番号の並び方（nx=1, ny=2の場合）】
    !   5─6   ← j=2 (y=y_max)
    !   |╲|
    !   3─4   ← j=1
    !   |╲|
    !   1─2   ← j=0 (y=y_min)
    !   左下から右へ、下段から上段へ連番を振る。
    !
    ! 【引数】
    !   nx, ny         : x・y方向の分割数            [in]
    !   x_min..y_max   : 領域の範囲                  [in]
    !   nodes          : 節点座標配列 (n_nodes)       [out]
    !   elems          : コネクティビティ (n_elems×3) [out]
    !==========================================================================
    subroutine generate_mesh(nx, ny, x_min, x_max, y_min, y_max, nodes, elems)
        integer,          intent(in)  :: nx, ny
        double precision, intent(in)  :: x_min, x_max, y_min, y_max
        type(pair),       intent(out) :: nodes(:)
        integer,          intent(out) :: elems(:,:)
        double precision :: dx, dy
        integer :: i, j, node_id, elem_id, n1, n2, n3, n4

        ! ── 格子間隔の計算 ──
        ! 領域をnx, ny等分したときの1格子あたりの幅・高さ
        dx = (x_max - x_min) / dble(nx)
        dy = (y_max - y_min) / dble(ny)

        ! ── 節点座標の生成 ──
        ! y方向（外ループ）→ x方向（内ループ）の順で番号を振る。
        ! これにより節点番号が左下から右上へ連番になる。
        node_id = 0
        do j = 0, ny
            do i = 0, nx
                node_id = node_id + 1
                ! x座標：x_minからdx刻みでi番目
                nodes(node_id)%x = x_min + dble(i) * dx
                ! y座標：y_minからdy刻みでj番目
                nodes(node_id)%y = y_min + dble(j) * dy
            end do
        end do

        ! ── 要素コネクティビティの生成 ──
        ! 各矩形セル(i,j)の4隅の節点番号を計算し、対角線で2三角形に分割する。
        !
        !   n3─n4      対角線はn1─n4（↗方向）
        !   |╲ |
        !   n1─n2
        !
        elem_id = 0
        do j = 1, ny
            do i = 1, nx
                ! 矩形の4隅の節点番号を計算
                n1 = (j - 1) * (nx + 1) + i   ! 左下
                n2 = n1 + 1                     ! 右下
                n3 = n1 + (nx + 1)              ! 左上
                n4 = n3 + 1                     ! 右上

                ! 三角形1：右下三角形（n1, n2, n4）
                ! n1→n2→n4の順（反時計回り）
                elem_id = elem_id + 1
                elems(elem_id, 1) = n1
                elems(elem_id, 2) = n2
                elems(elem_id, 3) = n4

                ! 三角形2：左上三角形（n1, n4, n3）
                ! n1→n4→n3の順
                elem_id = elem_id + 1
                elems(elem_id, 1) = n1
                elems(elem_id, 2) = n4
                elems(elem_id, 3) = n3
            end do
        end do

    end subroutine generate_mesh


    !==========================================================================
    ! subroutine: calc_elem_stiffness
    !
    ! 【目的】
    !   1次三角形要素の剛性行列 k_mat(3×3) を計算する。
    !   弱形式から導かれる定義：
    !     k_mat(i,j) = ∫_Ωe (∇Ni · ∇Nj) dA
    !
    ! 【1次三角形要素の特徴】
    !   形状関数 N が1次式なので ∇N が要素内で定数になる。
    !   そのため積分が「定数 × 面積」に簡略化できる。
    !   （高次要素ではGauss積分が必要になる）
    !
    ! 【引数】
    !   p1, p2, p3 : 三角形の3節点座標  [in]
    !   k_mat      : 要素剛性行列 3×3  [out]
    !==========================================================================
    subroutine calc_elem_stiffness(p1, p2, p3, k_mat)
        type(pair),       intent(in)  :: p1, p2, p3
        double precision, intent(out) :: k_mat(3,3)
        double precision :: x(3), y(3)
        double precision :: mat_J(2,2), mat_L(2,2)
        double precision :: detJ, val_m, area
        double precision :: v(3,2), res(3,2)
        integer :: i, j_global_nodes

        ! 節点座標を配列にセット（以降の計算でループを使いやすくするため）
        x(1) = p1%x;  y(1) = p1%y
        x(2) = p2%x;  y(2) = p2%y
        x(3) = p3%x;  y(3) = p3%y

        ! ── STEP1: ヤコビアン行列 J の構成 ──
        ! 局所座標(ξ,η)から全体座標(x,y)への変換行列。
        ! x = x1*N1 + x2*N2 + x3*N3 をξ,ηで偏微分すると：
        !   ∂x/∂ξ = x2-x1,  ∂x/∂η = x3-x1
        !   ∂y/∂ξ = y2-y1,  ∂y/∂η = y3-y1
        mat_J(1,1) = -x(1) + x(2)   ! ∂x/∂ξ
        mat_J(1,2) = -y(1) + y(2)   ! ∂y/∂ξ
        mat_J(2,1) = -x(1) + x(3)   ! ∂x/∂η
        mat_J(2,2) = -y(1) + y(3)   ! ∂y/∂η

        ! ── STEP2: 行列式・面積・逆行列の計算 ──
        ! detJは局所→全体座標変換の面積拡大率。
        ! 全体三角形の面積 = |detJ| × 局所三角形の面積(=1/2)
        ! absを取るのは、節点の並び順によってdetJが負になりうるため。
        detJ  = mat_J(1,1)*mat_J(2,2) - mat_J(1,2)*mat_J(2,1)
        val_m = 1.0d0 / detJ          ! 逆行列の係数
        area  = abs(detJ) * 0.5d0     ! 三角形の実面積（常に正）

        ! 2×2行列の逆行列の解析式：
        !   J⁻¹ = (1/detJ) * [  J22  -J12 ]
        !                     [ -J21   J11 ]
        mat_L(1,1) =  val_m * mat_J(2,2)
        mat_L(1,2) = -val_m * mat_J(1,2)
        mat_L(2,1) = -val_m * mat_J(2,1)
        mat_L(2,2) =  val_m * mat_J(1,1)

        ! ── STEP3: 形状関数の局所座標微分 v を設定 ──
        ! 1次三角形要素の形状関数：
        !   N1 = 1-ξ-η,  N2 = ξ,  N3 = η
        ! ξ,ηで偏微分すると定数になる（1次式なので微分後にξ,ηが残らない）：
        !   v(m,1) = ∂Nm/∂ξ,  v(m,2) = ∂Nm/∂η
        v(1,1) = -1.0d0;  v(1,2) = -1.0d0   ! N1の微分
        v(2,1) =  1.0d0;  v(2,2) =  0.0d0   ! N2の微分
        v(3,1) =  0.0d0;  v(3,2) =  1.0d0   ! N3の微分

        ! ── STEP4: 全体座標微分 res = J⁻¹ · v ──
        ! 連鎖律より v = J · res なので、両辺にJ⁻¹を掛けて res = J⁻¹ · v
        ! これにより局所座標微分を全体座標微分に変換する。
        !   res(i,1) = ∂Ni/∂x（x方向の勾配）
        !   res(i,2) = ∂Ni/∂y（y方向の勾配）
        do i = 1, 3
            res(i,1) = mat_L(1,1)*v(i,1) + mat_L(1,2)*v(i,2)   ! ∂Ni/∂x
            res(i,2) = mat_L(2,1)*v(i,1) + mat_L(2,2)*v(i,2)   ! ∂Ni/∂y
        end do

        ! ── STEP5: 要素剛性行列 k_mat の計算 ──
        ! k_mat(i,j) = ∫_Ωe (∇Ni · ∇Nj) dA
        ! ∇N が要素内で定数なので積分の外に出せる：
        !   = (∂Ni/∂x)(∂Nj/∂x) + (∂Ni/∂y)(∂Nj/∂y)) × area
        ! 対称行列になる（∇Ni·∇Nj = ∇Nj·∇Ni）。
        do i = 1, 3
            do j_global_nodes = 1, 3
                k_mat(i, j_global_nodes) = (res(i,1)*res(j_global_nodes,1) &
                                 + res(i,2)*res(j_global_nodes,2)) * area
            end do
        end do

    end subroutine calc_elem_stiffness


    !==========================================================================
    ! subroutine: coo_to_crs
    !
    ! 【目的】
    !   COO形式（タプルの羅列）をCRS形式（圧縮行格納）に変換する。
    !   アセンブリ時に重複して追加された同じ(row,col)のエントリを
    !   加算しながら圧縮し、メモリ効率の良いCRS形式を構築する。
    !
    ! 【変換の3ステップ】
    !   Step1: (row,col)の辞書順でソート → 重複が隣接するようにする
    !   Step2: 隣接比較しながら重複を加算・圧縮 → ユニークエントリ確定
    !   Step3: ソート済み行番号から row_ptr を構築
    !
    ! 【引数】
    !   coo_row, coo_col, coo_val : COO配列（重複込み）  [in]
    !   nnz_in                    : COOエントリ総数       [in]
    !   n                         : 行列のサイズ（節点数）[in]
    !   row_ptr, col_idx, values  : CRS配列               [out]
    !==========================================================================
   subroutine coo_to_crs(coo_row, coo_col, coo_val, nnz_in, &
                       n, row_ptr, col_idx, values)
    integer,          intent(in)  :: coo_row(:), coo_col(:), nnz_in, n
    double precision, intent(in)  :: coo_val(:)
    integer,          allocatable, intent(out) :: row_ptr(:), col_idx(:)
    double precision, allocatable, intent(out) :: values(:)

    integer          :: i, k, col, nnz
    integer          :: used_cols(20)
    integer          :: used_count
    integer,          allocatable :: pos(:)
    integer,          allocatable :: sorted_col(:)
    integer,          allocatable :: row_ptr_orig(:)
    double precision, allocatable :: sorted_val(:)
    double precision, allocatable :: dense_buffer(:)

    ! ── Step1: counting sort で row_ptr を構築 ──
    allocate(row_ptr(n+1))
    row_ptr = 0

    do k = 1, nnz_in
        row_ptr(coo_row(k) + 1) = row_ptr(coo_row(k) + 1) + 1
    end do

    row_ptr(1) = 1
    do i = 1, n
        row_ptr(i+1) = row_ptr(i+1) + row_ptr(i)
    end do

    ! ── Step2: COOエントリを行ごとに仕分け ──
    allocate(pos(n))
    allocate(sorted_col(nnz_in))
    allocate(sorted_val(nnz_in))

    pos = row_ptr(1:n)

    do k = 1, nnz_in
        i              = coo_row(k)
        sorted_col(pos(i)) = coo_col(k)
        sorted_val(pos(i)) = coo_val(k)
        pos(i)         = pos(i) + 1
    end do
    deallocate(pos)

    ! ── Step3: dense buffer で重複を加算・圧縮 ──
    allocate(row_ptr_orig(n+1))
    row_ptr_orig = row_ptr   ! 仕分け用row_ptrを保存

    allocate(dense_buffer(n))
    allocate(col_idx(nnz_in))
    allocate(values(nnz_in))
    dense_buffer = 0.0d0

    nnz        = 0
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
            col       = used_cols(k)
            nnz       = nnz + 1
            col_idx(nnz) = col
            values(nnz)  = dense_buffer(col)
            dense_buffer(col) = 0.0d0
        end do

        row_ptr(i+1) = nnz + 1
    end do

    deallocate(row_ptr_orig, dense_buffer, sorted_col, sorted_val)

end subroutine coo_to_crs
    !==========================================================================
    subroutine sort_coo(coo_row, coo_col, order, n)
        integer, intent(in)    :: coo_row(:), coo_col(:), n
        integer, intent(inout) :: order(:)
        integer :: i, j, tmp

        do i = 2, n
            tmp = order(i)   ! 今から挿入するエントリの元インデックス
            j   = i - 1      ! 既ソート部分の末尾から比較を始める

            ! tmp が指すエントリより大きいエントリを右に1つずつずらす
            ! 比較基準：まず行番号、同行なら列番号
            do while (j >= 1 .and. &
                      (coo_row(order(j)) > coo_row(tmp) .or. &
                      (coo_row(order(j)) == coo_row(tmp) .and. &
                       coo_col(order(j)) >  coo_col(tmp))))
                order(j+1) = order(j)   ! 1つ右にずらす
                j = j - 1
            end do
            order(j+1) = tmp   ! 適切な位置に挿入
        end do

    end subroutine sort_coo


    !==========================================================================
    ! subroutine: crs_subtract_col
    !
    ! 【目的】
    !   Dirichlet境界条件適用のStep A：右辺ベクトルの修正。
    !   境界節点 col_target の既知値 bc_val を使って
    !   F(row) -= K(row, col_target) * bc_val を計算する。
    !
    ! 【なぜこの処理が必要か】
    !   全体方程式 K u = F を既知節点(b)と未知節点(a)に分けると：
    !     Kaa ua = Fa - Kab ub
    !   「Kab ub を右辺に移項する」のがこの処理。
    !   密行列版では K_global(row, i) に直接アクセスできたが、
    !   CRS版では列方向の直接アクセスができないので行を走査して探す。
    !
    ! 【引数】
    !   row        : 修正する行番号           [in]
    !   col_target : 境界節点のグローバル番号 [in]
    !   bc_val     : 境界節点の既知値         [in]
    !   values, col_idx, row_ptr : CRS配列   [in]
    !   F_vec      : 右辺ベクトル（修正対象）[inout]
    !==========================================================================
    subroutine crs_subtract_col(row, col_target, bc_val, &
                                 values, col_idx, row_ptr, F_vec)
        integer,          intent(in)    :: row, col_target
        double precision, intent(in)    :: bc_val
        double precision, intent(in)    :: values(:)
        integer,          intent(in)    :: col_idx(:), row_ptr(:)
        double precision, intent(inout) :: F_vec(:)
        integer :: k

        ! row行のCRSエントリを走査して、列番号が col_target のものを探す
        ! row_ptr(row) ～ row_ptr(row+1)-1 が row行の非ゼロエントリの範囲
        do k = row_ptr(row), row_ptr(row+1)-1
            if (col_idx(k) == col_target) then
                ! 発見：右辺から K(row, col_target)*bc_val を引く
                F_vec(row) = F_vec(row) - values(k) * bc_val
                return   ! 同じ列は1つしかないので見つかったら即終了
            end if
        end do
        ! 見つからない場合：K(row, col_target) = 0 なので何もしない

    end subroutine crs_subtract_col


    !==========================================================================
    ! subroutine: apply_dirichlet_row
    !
    ! 【目的】
    !   Dirichlet境界条件適用のStep B：境界節点の行を単位行列化。
    !   node行のCRS値を「対角=1、その他=0」に書き換え、
    !   F(node) = bc_val とすることで
    !   「1 × u(node) = bc_val → u(node) = bc_val」の方程式にする。
    !
    ! 【なぜStep Aの後にやるか】
    !   Step Bでゼロにする前に、Step Aで他の行のFを修正しておく必要がある。
    !   順序が逆になると K(k, node) がすでにゼロになっており
    !   Step Aの修正が効かなくなる。
    !
    ! 【引数】
    !   node     : 境界節点のグローバル番号           [in]
    !   n        : 節点数（未使用だが対称性のために残す）[in]
    !   values   : CRS値配列（書き換え対象）           [inout]
    !   col_idx, row_ptr : CRS配列                    [in]
    !   F_vec    : 右辺ベクトル                        [inout]
    !   bc_val   : 境界節点の既知値                    [in]
    !==========================================================================
    subroutine apply_dirichlet_row(node, n, values, col_idx, row_ptr, F_vec, bc_val)
        integer,          intent(in)    :: node, n
        double precision, intent(inout) :: values(:), F_vec(:)
        integer,          intent(in)    :: col_idx(:), row_ptr(:)
        double precision, intent(in)    :: bc_val
        integer :: k

        ! node行のすべてのCRSエントリを走査する
        ! row_ptr(node) ～ row_ptr(node+1)-1 が node行の非ゼロエントリの範囲
        do k = row_ptr(node), row_ptr(node+1)-1
            if (col_idx(k) == node) then
                ! 対角成分：1にセット（u(node) の係数）
                values(k) = 1.0d0
            else
                ! 非対角成分：0にセット（他の節点との結合を切る）
                values(k) = 0.0d0
            end if
        end do

        ! 右辺を既知値にセット
        ! これで node行の方程式が「1 × u(node) = bc_val」になる
        F_vec(node) = bc_val

    end subroutine apply_dirichlet_row


    !==========================================================================
    ! subroutine: crs_matvec
    !
    ! 【目的】
    !   CRS形式の行列とベクトルの積 y = A·x を計算する。
    !   CG法の中で毎反復呼ばれる最重要演算。
    !
    ! 【CRS形式の利点】
    !   非ゼロ要素のみをループするので O(nnz) で計算できる。
    !   密行列の O(n²) と比べ、FEMのような疎行列では大幅に速い。
    !   （例：n=20万節点、nnz=140万 → 密行列比で約28倍速い）
    !
    ! 【引数】
    !   n                        : 行列・ベクトルのサイズ [in]
    !   values, col_idx, row_ptr : CRS配列               [in]
    !   x                        : 入力ベクトル           [in]
    !   y                        : 出力ベクトル y = A·x  [out]
    !==========================================================================
    subroutine crs_matvec(n, values, col_idx, row_ptr, x, y)
        integer,          intent(in)  :: n
        double precision, intent(in)  :: values(:), x(:)
        integer,          intent(in)  :: col_idx(:), row_ptr(:)
        double precision, intent(out) :: y(:)
        integer :: i, k

        do i = 1, n
            y(i) = 0.0d0   ! y(i)の初期化（前の計算結果が残らないようにする）

            ! i行の非ゼロエントリだけをループ
            ! row_ptr(i)   : i行の最初のエントリのインデックス
            ! row_ptr(i+1)-1 : i行の最後のエントリのインデックス
            do k = row_ptr(i), row_ptr(i+1)-1
                ! values(k) = A(i, col_idx(k))
                ! x(col_idx(k)) = x の col_idx(k) 列目
                y(i) = y(i) + values(k) * x(col_idx(k))
            end do
        end do

    end subroutine crs_matvec


    !==========================================================================
    ! subroutine: cg_solver
    !
    ! 【目的】
    !   共役勾配法（CG法）で連立方程式 K u = F を解く。
    !   Gauss消去法（O(n³)）の代わりに行列ベクトル積だけで解く反復解法。
    !   対称正定値行列にのみ使える。FEMの剛性行列はこれを満たす。
    !
    ! 【CG法のアルゴリズム概要】
    !   残差 r = b - A x を最小化する方向 p を更新しながら解に近づく。
    !   理論上は n ステップで厳密解に到達するが、
    !   実際は残差が閾値以下になった時点で打ち切る。
    !
    ! 【計算量】
    !   O(nnz × iter)：1反復あたり crs_matvec（O(nnz)）が1回
    !   Gauss消去法 O(n³) に対して、疎行列では格段に速い。
    !
    ! 【引数】
    !   n                        : 節点数               [in]
    !   values, col_idx, row_ptr : CRS行列              [in]
    !   b                        : 右辺ベクトル F       [in]
    !   x                        : 解ベクトル u（初期値から更新）[inout]
    !==========================================================================
    subroutine cg_solver(n, values, col_idx, row_ptr, b, x)
        integer,          intent(in)    :: n
        double precision, intent(in)    :: values(:), b(:)
        integer,          intent(in)    :: col_idx(:), row_ptr(:)
        double precision, intent(inout) :: x(:)

        double precision :: r(n), p(n), Ap(n)
        double precision :: alpha, beta, rr, rr_new
        integer :: iter, ii, kk

        ! ── 初期化 ──
        ! 初期残差 r = b - A*x を計算する
        ! x=0 なら r=b だが、Dirichlet節点の初期値が入っている場合は
        ! A*x を引く必要があるので crs_matvec を使う
        call crs_matvec(n, values, col_idx, row_ptr, x, Ap)
        r  = b - Ap    ! 初期残差
        p  = r         ! 初期探索方向（最初は残差と同じ方向）
        rr = dot_product(r, r)   ! r·r（収束判定と次ステップのαに使う）

        ! ── 反復ループ ──
        do iter = 1, 10000
            ! 収束判定：残差のノルム ||r|| が閾値以下なら終了
            if (sqrt(rr) < 1.0d-10) exit

            ! Ap = A·p（探索方向 p に行列を掛ける）
            ! これがCG法で唯一の行列演算→ CRS形式を活かす
            call crs_matvec(n, values, col_idx, row_ptr, p, Ap)

            ! ステップ幅 α の計算
            ! α = (r·r) / (p·Ap)
            ! この α で進むと残差が最も小さくなる（最適ステップ幅）
            alpha  = rr / dot_product(p, Ap)

            ! 解の更新：x = x + α·p
            x      = x + alpha * p

            ! 残差の更新：r = r - α·(A·p)
            ! 直接 b - A·x を計算し直すより数値誤差が小さい
            r      = r - alpha * Ap

            ! 新しい残差の内積
            rr_new = dot_product(r, r)

            ! 探索方向の更新係数 β = (r_new·r_new) / (r·r)
            ! β によって新しい p は古い p と共役（直交）になる
            ! これがCG法の「共役」の意味
            beta   = rr_new / rr

            ! 探索方向の更新：p = r + β·p
            ! 新しい残差方向に、前の探索方向を β 倍加えることで
            ! 既探索部分空間を再探索しない（効率化）
            p      = r + beta * p

            ! 次の反復のために更新
            rr     = rr_new
        end do

        ! 収束情報の出力
        print '(A,I5,A,ES10.3)', "  CG: iter =", iter-1, &
                                  "  residual =", sqrt(rr)

    end subroutine cg_solver


end program solve_laplace_fem