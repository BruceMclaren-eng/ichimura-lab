!=============================================================================
! ラプラス方程式の有限要素法 (FEM) 解法  ── 写経・学習用
!
! 【支配方程式】
!   ∇²u = ∂²u/∂x² + ∂²u/∂y² = 0    (0 ≤ x ≤ 1, 0 ≤ y ≤ 2)
!
! 【境界条件】
!   u(x, 0) = 0   下辺 (Dirichlet)
!   u(x, 2) = 1   上辺 (Dirichlet)
!   左辺・右辺    自由表面 (Neumann: ∂u/∂n = 0)
!
! 【FEMの基本的な考え方】
!   「真の解 u を全点で求めるのは無理なので、
!    有限個の節点値で近似解 ũ を作り、
!    残差 r = ∇²ũ の重み付き平均がゼロになるよう節点値を決める」
!
!   この考え方を「重み付き残差法 (Galerkin 法)」と呼ぶ。
!   ∫ Nj · r dΩ = 0  を全節点 j について解くと、
!   Green の定理 (部分積分) を経て
!   ∫ (∇Ni · ∇Nj) dΩ = 境界積分  という弱形式が導かれる。
!
! 【全体フロー】
!   1. メッシュ生成       : 節点座標 + 要素コネクティビティ
!   2. アセンブル         : 各要素の K^e (3×3) を全体 K (n×n) に足し込む
!   3. 境界条件の適用     : Dirichlet 節点を処理
!   4. 連立方程式の求解   : Gauss 消去法
!   5. 結果出力
!
! 【要素タイプ】  2次元 1次三角形要素 (CST: Constant Strain Triangle)
! 【求解手法】    Gauss の消去法
!=============================================================================

!　もともとの計算時間は1.715秒程度（コンパイル含む）。コンパイル以降は、0.365秒
!=============================================================================
! モジュール: fem_types
!
! 【目的】
!   プログラム全体で共通に使う型をまとめる。
!   x 座標と y 座標をペアで扱う構造体 pair を定義する。
!
! 【使い方】
!   type(pair) :: p
!   p%x = 1.0d0   ! x 成分へのアクセス
!   p%y = 2.0d0   ! y 成分へのアクセス
!=============================================================================
module fem_types
    implicit none

    ! 2次元座標 (x, y) を一組として扱う構造体
    type :: pair
        double precision :: x, y
    end type pair

end module fem_types


!=============================================================================
! メインプログラム: solve_laplace_fem
!=============================================================================
program solve_laplace_fem
    use fem_types
    implicit none

    !-------------------------------------------------------------------------
    ! 変数宣言
    !
    ! 【メッシュ情報】
    !   nx, ny   : x・y 方向の分割数
    !              例) nx=1, ny=2 のとき 2×3 の格子点 → 6 節点
    !   n_nodes  : 全節点数 = (nx+1)*(ny+1)
    !   n_elems  : 全要素数 = nx*ny*2  (各矩形を対角線で 2 三角形に分割)
    !   idx(3)   : 1 要素を構成する 3 節点のグローバル番号を一時保持する配列
    !
    ! 【動的配列】
    !   nodes    : 各節点の (x, y) 座標  [n_nodes]
    !   elems    : 要素コネクティビティ   [n_elems, 3]
    !              elems(e, 1:3) = 要素 e の 3 節点のグローバル番号
    !   K_global : 全体剛性行列           [n_nodes, n_nodes]
    !   F_vec    : 右辺ベクトル           [n_nodes]
    !   U_vec    : 解ベクトル (節点値)    [n_nodes]
    !   is_fixed : Dirichlet 境界フラグ   [n_nodes]
    !              .true. = その節点は値が既知 (境界条件で固定)
    !
    ! 【作業用変数】
    !   K_local  : 1 要素分の剛性行列 (3×3)
    !              アセンブル前の中間計算に使う
    !   bc_value : 境界節点に課す既知値
    !-------------------------------------------------------------------------
    integer :: nx, ny
    integer :: n_nodes, n_elems
    integer :: i, j, k, e
    integer :: idx(3)

    double precision :: x_min, x_max, y_min, y_max

    type(pair), allocatable :: nodes(:)
    integer,    allocatable :: elems(:,:)
    double precision, allocatable :: K_global(:,:)
    double precision, allocatable :: F_vec(:)
    double precision, allocatable :: U_vec(:)

    double precision :: K_local(3,3)
    double precision :: bc_value
    logical, allocatable :: is_fixed(:)


    !=========================================================================
    ! STEP 1: 解析パラメータの設定
    !
    ! 【領域】  x ∈ [0,1],  y ∈ [0,2]
    ! 【分割】  nx=1, ny=2 → Fig.1 のメッシュ (6 節点, 4 要素) と一致
    !
    ! 【節点番号の並び方 (nx=1, ny=2 の場合)】
    !
    !   節点5 ─ 節点6    y = 2
    !     |  ╲ 4 |
    !     | 3  ╲ |
    !   節点3 ─ 節点4    y = 1
    !     |  ╲ 2 |
    !     | 1  ╲ |
    !   節点1 ─ 節点2    y = 0
    !     x=0    x=1
    !
    !   左下から右へ、下段から上段へ順番に番号が付く。
    !=========================================================================
    x_min = 0.0d0;  x_max = 1.0d0
    y_min = 0.0d0;  y_max = 2.0d0
    nx = 100
    ny = 200

    n_nodes = (nx + 1) * (ny + 1)   ! = 6
    n_elems = nx * ny * 2            ! = 4


    !=========================================================================
    ! STEP 2: 配列のメモリ確保と初期化
    !
    ! Fortran の動的配列は allocate で確保してから使う。
    ! 確保直後の値は不定なので、明示的にゼロや .false. で初期化する。
    !=========================================================================
    allocate(nodes(n_nodes))
    allocate(elems(n_elems, 3))
    allocate(K_global(n_nodes, n_nodes))
    allocate(F_vec(n_nodes))
    allocate(U_vec(n_nodes))
    allocate(is_fixed(n_nodes))

    K_global = 0.0d0
    F_vec    = 0.0d0
    U_vec    = 0.0d0
    is_fixed = .false.


    !=========================================================================
    ! STEP 3: メッシュ生成
    !
    ! generate_mesh が nodes(:) と elems(:,:) を埋める。
    ! 詳細は contains 内のサブルーチンを参照。
    !=========================================================================
    call generate_mesh(nx, ny, x_min, x_max, y_min, y_max, nodes, elems)

    print *, "--- Mesh Generated: Nodes =", n_nodes, " Elements =", n_elems


    !=========================================================================
    ! STEP 4: 全体剛性行列のアセンブル
    !
    ! 【弱形式から導かれる剛性行列の定義】
    !   K[I,J] = ∫_Ω (∇N_I · ∇N_J) dΩ
    !
    !   ただし全体では N_I は全節点に定義された関数なので、
    !   要素をまたがない (要素内でのみ非ゼロ) 性質を利用して
    !   要素ごとに計算し足し合わせる (アセンブル)。
    !
    ! 【アセンブルの手順】
    !   各要素 e について:
    !     1) 要素 e を構成する 3 節点のグローバル番号 idx(1:3) を取得
    !     2) 3 節点の座標から K_local(3×3) を計算
    !     3) K_local(i,j) を K_global(idx(i), idx(j)) に加算
    !
    !   ローカル番号 i (1,2,3) → グローバル番号 idx(i)
    !   という変換がアセンブルの核心。
    !=========================================================================
    do e = 1, n_elems
        idx = elems(e, :)   ! この要素の 3 節点グローバル番号

        call calc_elem_stiffness( &
            nodes(idx(1)), nodes(idx(2)), nodes(idx(3)), K_local)

        ! K_local をグローバル番号の位置に足し込む
        do i = 1, 3
            do j = 1, 3
                K_global(idx(i), idx(j)) = &
                    K_global(idx(i), idx(j)) + K_local(i, j)
            end do
        end do
    end do


    !=========================================================================
    ! STEP 5: Dirichlet 境界条件の適用
    !
    ! 【境界条件】
    !   u(x, 0) = 0  → 節点 1, 2  (y ≈ 0)
    !   u(x, 2) = 1  → 節点 5, 6  (y ≈ 2)
    !
    ! 【2ステップ法の理由】
    !   全体方程式 K u = f を既知 (b) と未知 (a) に分けると:
    !   [ Kaa Kab ] [ ua ]   [ fa ]
    !   [ Kba Kbb ] [ ub ] = [ fb ]
    !
    !   未知節点の方程式:  Kaa ua = fa - Kab ub
    !   → 右辺に ub の寄与を移項する操作が Step A
    !   → Kbb の行・列を単位行列にする操作が Step B
    !   (対称行列の構造を保ったまま境界条件を組み込める)
    !
    ! 【Step A と Step B の順序が重要】
    !   Step B で K の列をゼロにする前に、
    !   Step A で F を修正しておかなければならない。
    !   順序を逆にすると K(k,i) がすでにゼロになっていて
    !   F の修正が効かない。
    !=========================================================================
    print *, "--- Applying Boundary Conditions ---"

    ! ── Step A: 既知節点値を右辺ベクトル F に移項 ──
    !
    ! 境界節点 i の既知値を ū とすると、
    ! 全行 k に対して  F(k) -= K(k,i) * ū  を計算する。
    ! これは「K ub の寄与を右辺に引く」操作。
    do i = 1, n_nodes
        bc_value = 0.0d0   ! 固定値の初期値 (Neumann 節点はここで終わる)

        if (abs(nodes(i)%y - 0.0d0) < 1.0d-6) then   ! 下辺: u = 0
            is_fixed(i) = .true.
            bc_value    = 0.0d0
        else if (abs(nodes(i)%y - 2.0d0) < 1.0d-6) then   ! 上辺: u = 1
            is_fixed(i) = .true.
            bc_value    = 1.0d0
        end if

        if (is_fixed(i)) then
            U_vec(i) = bc_value   ! 解配列に既知値をセット

            ! 全行に対して右辺を修正
            do k = 1, n_nodes
                F_vec(k) = F_vec(k) - K_global(k, i) * bc_value
            end do
        end if
    end do

    ! ── Step B: K_global の境界節点の行・列を単位行列化 ──
    !
    ! 境界節点 i の行・列をゼロにして対角を 1 にする。
    ! F(i) = ū にセットすることで、
    ! i 行の方程式が  1 * u(i) = ū  → u(i) = ū  になる。
    do i = 1, n_nodes
        if (is_fixed(i)) then
            K_global(i, :) = 0.0d0   ! i 行をゼロクリア
            K_global(:, i) = 0.0d0   ! i 列をゼロクリア (対称性維持)
            K_global(i, i) = 1.0d0   ! 対角を 1 に
            F_vec(i)       = U_vec(i)  ! 右辺に既知値をセット
        end if
    end do


    !=========================================================================
    ! STEP 6: 連立一次方程式 K u = F の求解
    !
    ! Gauss の消去法 (前進消去 + 後退代入) で解く。
    ! FEM の剛性行列は正定値対称なので、
    ! ピボット選択なしでも数値的に安定して動作する。
    !=========================================================================
    print *, "--- Solving Linear System (Gaussian Elimination) ---"
    call solve_linear_system(n_nodes, K_global, F_vec, U_vec)


    !=========================================================================
    ! STEP 7: 結果の出力
    !
    ! 【解析解との比較】
    !   この問題の解析解は u(x, y) = y/2
    !   (Neumann 条件から x 依存がなく、y 方向線形)
    !   → 節点 3 (y=1): u = 0.5
    !   → 節点 4 (y=1): u = 0.5
    !   線形三角形要素は線形関数を厳密に再現できるので
    !   FEM 解 = 解析解 が期待される。
    !=========================================================================
    print *, "========================================"
    print *, " FINAL RESULTS (Nodal Values)"
    print *, "========================================"
    print *, " Node ID |   X      Y    |   U(x,y)   "
    print *, "---------+---------------+------------"
    do i = 1, n_nodes
        write(*, '(I6, "   |", F6.2, " ", F6.2, "  |", F10.5)') &
            i, nodes(i)%x, nodes(i)%y, U_vec(i)
    end do

    ! メモリ解放 (allocate したものは deallocate で返す)
    deallocate(nodes, elems, K_global, F_vec, U_vec, is_fixed)


contains  ! ───────────────── 内部サブルーチン ─────────────────


    !==========================================================================
    ! サブルーチン: generate_mesh
    !
    ! 【目的】
    !   矩形領域を nx×ny の格子に分割し、各矩形を対角線で
    !   2 つの三角形要素に分けて nodes と elems を生成する。
    !
    ! 【節点番号の付け方】
    !   j (y 方向) の外ループ → i (x 方向) の内ループ の順で
    !   左下から右上へ連番を振る。
    !
    !   例) nx=2, ny=2 の場合:
    !     7─8─9    j=2
    !     │╲│╲│
    !     4─5─6    j=1
    !     │╲│╲│
    !     1─2─3    j=0
    !
    ! 【要素の分割】
    !   各矩形の 4 隅を n1(左下), n2(右下), n3(左上), n4(右上) として
    !   対角線 n1─n4 で 2 三角形に分割する:
    !
    !   n3─n4
    !   │ ╲│
    !   n1─n2
    !
    !   要素 (奇数): n1, n2, n4  ← 右下三角形
    !   要素 (偶数): n1, n4, n3  ← 左上三角形
    !
    ! 【引数】
    !   nx, ny         : x・y 方向の分割数                  [in]
    !   x_min .. y_max : 領域の範囲                          [in]
    !   nodes          : 節点座標配列 (サイズ n_nodes)       [out]
    !   elems          : コネクティビティ配列 (n_elems × 3)  [out]
    !==========================================================================
    subroutine generate_mesh(nx, ny, x_min, x_max, y_min, y_max, nodes, elems)
        integer,          intent(in)  :: nx, ny
        double precision, intent(in)  :: x_min, x_max, y_min, y_max
        type(pair),       intent(out) :: nodes(:)
        integer,          intent(out) :: elems(:,:)

        double precision :: dx, dy
        integer :: i, j, node_id, elem_id, n1, n2, n3, n4

        ! 格子間隔
        dx = (x_max - x_min) / dble(nx)
        dy = (y_max - y_min) / dble(ny)

        ! ── 節点座標の生成 ──
        ! y 方向 (外ループ) → x 方向 (内ループ) の順で番号を振る
        node_id = 0
        do j = 0, ny
            do i = 0, nx
                node_id = node_id + 1
                nodes(node_id)%x = x_min + dble(i) * dx
                nodes(node_id)%y = y_min + dble(j) * dy
            end do
        end do

        ! ── 要素コネクティビティの生成 ──
        ! 矩形セル (i,j) の 4 隅の節点番号を計算して 2 三角形に分割する
        !
        !   n3 = (j)行 (i)列     n4 = (j)行 (i+1)列
        !   n1 = (j-1)行 (i)列  n2 = (j-1)行 (i+1)列
        !
        elem_id = 0
        do j = 1, ny
            do i = 1, nx
                n1 = (j - 1) * (nx + 1) + i    ! 左下
                n2 = n1 + 1                      ! 右下
                n3 = n1 + (nx + 1)               ! 左上
                n4 = n3 + 1                      ! 右上

                ! 要素 1: 右下三角形 (n1, n2, n4)
                elem_id = elem_id + 1
                elems(elem_id, 1) = n1
                elems(elem_id, 2) = n2
                elems(elem_id, 3) = n4

                ! 要素 2: 左上三角形 (n1, n4, n3)
                elem_id = elem_id + 1
                elems(elem_id, 1) = n1
                elems(elem_id, 2) = n4
                elems(elem_id, 3) = n3
            end do
        end do

    end subroutine generate_mesh


    !==========================================================================
    ! サブルーチン: calc_elem_stiffness
    !
    ! 【目的】
    !   1次三角形要素の剛性行列 k_mat(3×3) を計算する。
    !
    ! 【理論的背景】
    !   弱形式から導かれる要素剛性行列の定義:
    !     k_mat(i,j) = ∫_Ωe (∇N_i · ∇N_j) dA
    !                = (∂N_i/∂x)(∂N_j/∂x) + (∂N_i/∂y)(∂N_j/∂y)) × A_e
    !
    !   1次三角形要素では ∇N_i が要素内で定数なので
    !   積分が「定数 × 面積」になり解析的に計算できる。
    !   (高次要素なら ∇N_i が x,y の関数になり Gauss 積分が必要になる)
    !
    ! 【局所座標と座標変換】
    !   形状関数は局所座標 (ξ, η) で定義される:
    !     N1(ξ,η) = 1 - ξ - η
    !     N2(ξ,η) = ξ
    !     N3(ξ,η) = η
    !
    !   全体座標 (x, y) との変換式:
    !     x = x1*N1 + x2*N2 + x3*N3
    !     y = y1*N1 + y2*N2 + y3*N3
    !
    !   この変換を ξ, η で微分するとヤコビアン J が出てくる:
    !     J = [ ∂x/∂ξ  ∂y/∂ξ ]   [ x2-x1  y2-y1 ]
    !         [ ∂x/∂η  ∂y/∂η ] = [ x3-x1  y3-y1 ]
    !
    ! 【計算手順 (5 ステップ)】
    !   1) 節点座標から J を構成
    !   2) det(J), 面積 A = |det J|/2, J の逆行列 J⁻¹ を計算
    !   3) 形状関数の局所座標微分 v を定義 (全要素共通の定数)
    !   4) 全体座標微分 res = J⁻¹ · v を計算 (連鎖律)
    !   5) k_mat(i,j) = (res_i · res_j) × A
    !
    ! 【引数】
    !   p1, p2, p3 : 三角形の 3 節点座標  [in]
    !   k_mat      : 要素剛性行列 3×3     [out]
    !==========================================================================
    subroutine calc_elem_stiffness(p1, p2, p3, k_mat)
        type(pair),       intent(in)  :: p1, p2, p3
        double precision, intent(out) :: k_mat(3,3)

        ! x(m), y(m) : m 番目のローカル節点の全体座標
        double precision :: x(3), y(3)

        ! mat_J(2,2) : ヤコビアン行列
        !   mat_J(1,1) = ∂x/∂ξ = x2-x1
        !   mat_J(1,2) = ∂y/∂ξ = y2-y1
        !   mat_J(2,1) = ∂x/∂η = x3-x1
        !   mat_J(2,2) = ∂y/∂η = y3-y1
        double precision :: mat_J(2,2)

        ! mat_L(2,2) : ヤコビアン逆行列 J⁻¹
        !   全体座標微分 ← 局所座標微分 の変換に使う
        double precision :: mat_L(2,2)

        ! detJ  : ヤコビアンの行列式 = (x2-x1)(y3-y1) - (y2-y1)(x3-x1)
        ! val_m : 1/detJ  (逆行列の係数)
        ! area  : 三角形の実面積 = |detJ| / 2
        double precision :: detJ, val_m, area

        ! v(m, k) : 形状関数 Nm の局所座標微分
        !   v(m, 1) = ∂Nm/∂ξ
        !   v(m, 2) = ∂Nm/∂η
        !   N1=1-ξ-η → [-1, -1]
        !   N2=ξ     → [ 1,  0]
        !   N3=η     → [ 0,  1]
        !   ※ Nm が ξ,η の 1 次式なので v は定数
        double precision :: v(3,2)

        ! res(m, k) : 形状関数 Nm の全体座標微分
        !   res(m, 1) = ∂Nm/∂x
        !   res(m, 2) = ∂Nm/∂y
        !   連鎖律より  res = J⁻¹ · v  で計算
        double precision :: res(3,2)

        integer :: i, j_idx

        ! 節点座標を配列にセット
        x(1) = p1%x;  y(1) = p1%y
        x(2) = p2%x;  y(2) = p2%y
        x(3) = p3%x;  y(3) = p3%y

        !----------------------------------------------------------------------
        ! STEP 1: ヤコビアン行列 J の構成
        !
        ! x = x1*N1 + x2*N2 + x3*N3 を ξ で微分:
        !   ∂x/∂ξ = -x1 + x2 = x2 - x1
        ! x を η で微分:
        !   ∂x/∂η = -x1 + x3 = x3 - x1
        ! y についても同様。
        !----------------------------------------------------------------------
        mat_J(1,1) = -x(1) + x(2)   ! x2 - x1  = ∂x/∂ξ
        mat_J(1,2) = -y(1) + y(2)   ! y2 - y1  = ∂y/∂ξ
        mat_J(2,1) = -x(1) + x(3)   ! x3 - x1  = ∂x/∂η
        mat_J(2,2) = -y(1) + y(3)   ! y3 - y1  = ∂y/∂η

        !----------------------------------------------------------------------
        ! STEP 2: det J, 面積 A, 逆行列 J⁻¹ の計算
        !
        ! 【行列式】
        !   2×2 行列の行列式: det J = J11*J22 - J12*J21
        !   det J の絶対値は、局所三角形を全体三角形に拡大するスケール。
        !
        ! 【面積】
        !   局所座標の基準三角形の面積 = 1/2
        !   全体三角形の面積 = |det J| * (局所の面積) = |det J| / 2
        !
        ! 【逆行列】
        !   2×2 の逆行列の解析式:
        !   J⁻¹ = (1/det J) * [  J22  -J12 ]
        !                      [ -J21   J11 ]
        !   J11 と J22 を入れ替え、J12 と J21 の符号を反転するだけ。
        !----------------------------------------------------------------------
        detJ  = mat_J(1,1)*mat_J(2,2) - mat_J(1,2)*mat_J(2,1)
        val_m = 1.0d0 / detJ
        area  = abs(detJ) * 0.5d0

        mat_L(1,1) =  val_m * mat_J(2,2)   !  J22 / det J
        mat_L(1,2) = -val_m * mat_J(1,2)   ! -J12 / det J
        mat_L(2,1) = -val_m * mat_J(2,1)   ! -J21 / det J
        mat_L(2,2) =  val_m * mat_J(1,1)   !  J11 / det J

        !----------------------------------------------------------------------
        ! STEP 3: 形状関数の局所座標微分 v を設定 (全要素共通)
        !
        ! N1 = 1 - ξ - η  を ξ で微分 → -1, η で微分 → -1
        ! N2 = ξ          を ξ で微分 →  1, η で微分 →  0
        ! N3 = η          を ξ で微分 →  0, η で微分 →  1
        !
        ! Nm が ξ,η の 1 次式 → 微分しても ξ,η が残らない → 定数になる
        ! この「v が定数」であることが、後で積分を面積×定数に簡略化できる根拠。
        !----------------------------------------------------------------------
        v(1,1) = -1.0d0;  v(1,2) = -1.0d0   ! ∂N1/∂ξ, ∂N1/∂η
        v(2,1) =  1.0d0;  v(2,2) =  0.0d0   ! ∂N2/∂ξ, ∂N2/∂η
        v(3,1) =  0.0d0;  v(3,2) =  1.0d0   ! ∂N3/∂ξ, ∂N3/∂η

        !----------------------------------------------------------------------
        ! STEP 4: 全体座標微分 res = J⁻¹ · v (連鎖律)
        !
        ! 【なぜ J⁻¹ を掛けるのか】
        !   連鎖律より:
        !     ∂N/∂ξ = (∂N/∂x)(∂x/∂ξ) + (∂N/∂y)(∂y/∂ξ)
        !     ∂N/∂η = (∂N/∂x)(∂x/∂η) + (∂N/∂y)(∂y/∂η)
        !
        !   行列形式で書くと:  v = J · res
        !   両辺に左から J⁻¹ を掛けると:  res = J⁻¹ · v
        !
        ! 【コードの実装】
        !   res(i,1) = ∂Ni/∂x = J⁻¹[1,1]*v[i,1] + J⁻¹[1,2]*v[i,2]
        !   res(i,2) = ∂Ni/∂y = J⁻¹[2,1]*v[i,1] + J⁻¹[2,2]*v[i,2]
        !----------------------------------------------------------------------
        do i = 1, 3
            res(i,1) = mat_L(1,1)*v(i,1) + mat_L(1,2)*v(i,2)   ! ∂Ni/∂x
            res(i,2) = mat_L(2,1)*v(i,1) + mat_L(2,2)*v(i,2)   ! ∂Ni/∂y
        end do

        !----------------------------------------------------------------------
        ! STEP 5: 要素剛性行列 k_mat の計算
        !
        ! 【定義】
        !   k_mat(i,j) = ∫_Ωe (∇Ni · ∇Nj) dA
        !
        ! 【なぜ定数×面積になるか】
        !   res (= ∇N) が要素内で定数なので積分の外に出せる:
        !     = (∇Ni · ∇Nj) × ∫_Ωe dA
        !     = (∇Ni · ∇Nj) × A
        !
        ! 【内積の展開】
        !   ∇Ni · ∇Nj = (∂Ni/∂x)(∂Nj/∂x) + (∂Ni/∂y)(∂Nj/∂y)
        !              = res(i,1)*res(j,1) + res(i,2)*res(j,2)
        !
        ! 【対称性】
        !   k_mat(i,j) = k_mat(j,i)  ← 内積は交換可能なので対称行列になる
        !----------------------------------------------------------------------
        do i = 1, 3
            do j_idx = 1, 3
                k_mat(i, j_idx) = ( res(i,1)*res(j_idx,1) &
                                  + res(i,2)*res(j_idx,2) ) * area
            end do
        end do

    end subroutine calc_elem_stiffness


    !==========================================================================
    ! サブルーチン: solve_linear_system
    !
    ! 【目的】
    !   Gauss の消去法で n×n 連立一次方程式 A x = b を解く。
    !
    ! 【アルゴリズム】
    !   1. 前進消去 (Forward Elimination)
    !      ピボット行 k を固定し、k より下の行から k 列成分を消去する。
    !      Row(i) = Row(i) - (A(i,k)/A(k,k)) × Row(k)
    !      → 上三角行列 U に変換される。
    !
    !   2. 後退代入 (Back Substitution)
    !      上三角行列になったので最終行から順に解を確定させる:
    !      x(i) = (b(i) - Σ_{j>i} A(i,j)*x(j)) / A(i,i)
    !
    ! 【なぜ安全か】
    !   FEM の剛性行列 (境界条件適用後) は正定値対称なので
    !   対角成分が必ず正になり、ゼロ除算が起きない。
    !
    ! 【注意】
    !   A_in, b_in をコピーして作業するので呼び出し元の値は変化しない。
    !   (Fortran では intent(in) でも配列の内容を書き換えようとすると
    !    コンパイルエラーになるが、ここは安全のためコピーして使う)
    !
    ! 【引数】
    !   n     : 方程式の数      [in]
    !   A_in  : 係数行列 n×n  [in]
    !   b_in  : 右辺ベクトル  [in]
    !   x_out : 解ベクトル    [out]
    !==========================================================================
    subroutine solve_linear_system(n, A_in, b_in, x_out)
        integer,          intent(in)  :: n
        double precision, intent(in)  :: A_in(n,n), b_in(n)
        double precision, intent(out) :: x_out(n)

        ! 作業用コピー (元の A_in, b_in を壊さない)
        double precision :: A(n,n), b(n)
        double precision :: factor, sum_val
        integer :: k, i, j

        A = A_in
        b = b_in

        !----------------------------------------------------------------------
        ! 1. 前進消去 (Forward Elimination)
        !
        ! k 行目をピボット行として固定し、
        ! k+1 行目以下から k 列の成分を消去していく。
        !
        ! 消去係数 factor = A(i,k) / A(k,k)
        ! → Row(i) = Row(i) - factor * Row(k)
        ! とすることで A(i,k) がゼロになる。
        !
        ! これを k=1 から k=n-1 まで繰り返すと
        ! 行列が上三角形になる。
        !----------------------------------------------------------------------
        do k = 1, n-1
            do i = k+1, n
                ! ゼロ割りチェック (FEM 行列では原則発生しない)
                if (abs(A(k,k)) < 1.0d-12) then
                    print *, "Error: Zero pivot at row", k
                    stop
                end if

                factor = A(i,k) / A(k,k)

                ! j = k から始めることで j < k は触らない (すでにゼロ)
                do j = k, n
                    A(i,j) = A(i,j) - factor * A(k,j)
                end do
                b(i) = b(i) - factor * b(k)
            end do
        end do

        !----------------------------------------------------------------------
        ! 2. 後退代入 (Back Substitution)
        !
        ! 上三角行列になったので、最終行 (n 行目) から解が確定する:
        !   x(n) = b(n) / A(n,n)
        !   x(n-1) = (b(n-1) - A(n-1,n)*x(n)) / A(n-1,n-1)
        !   ...
        !   x(i) = (b(i) - Σ_{j=i+1}^{n} A(i,j)*x(j)) / A(i,i)
        !----------------------------------------------------------------------
        x_out(n) = b(n) / A(n,n)

        do i = n-1, 1, -1
            sum_val = 0.0d0
            do j = i+1, n
                sum_val = sum_val + A(i,j) * x_out(j)
            end do
            x_out(i) = (b(i) - sum_val) / A(i,i)
        end do

    end subroutine solve_linear_system


end program solve_laplace_fem