module fem_types
    implicit none
    type :: pair
        double precision :: x,y
    end type pair

end module fem_types

program fem_practice
    use fem_types
    implicit none

    !【メッシュ情報】
    integer :: nx, ny ! x,y方向の分割数
    integer :: n_nodes ! 全ノード数
    integer :: n_elems　! 全要素数
    double precision :: x_min, x_max, y_min, y_max !領域の定義
    integer :: n_edges !メッシュの辺の数
    
    ! 【動的配列】
    type(pair), allocatable :: nodes(:)!ノードの座標
    integer, allocatable :: elems(:,:) !要素のノード（要素の番号,要素上の接点番号）
    double precision, allocatable :: K_global(:,:)!全体剛性行列
    double precision, allocatable ::F_vec(:)　!力ベクトル
    double precision, allocatable :: u_vec(:)!未知数ベクトル
    double precision, allocatable :: K_local(:,:)!要素剛性行列
    !===================================
    !【Step1:メッシュの作成】
    !===================================    
    ! 領域の定義
    x_min = 0.0d0
    x_max = 1.0d0
    y_min = 0.0d0
    y_max = 2.0d0
    !分割数の定義
    nx = 1
    ny = 2
    !全ノード数
    n_nodes = (nx + 1) * (ny + 1)
    !全要素数
    n_elems = nx * ny　* 2
    !メッシュの辺の数
    n_edges = 3


    !==================================
    !【Step2:配列の初期化】
    !==================================
    allocate(nodes(n_nodes))
    allocate(elems(n_elems,n_edges))
    allocate(K_global(n_nodes,n_nodes))
    allocate(F_vec(n_nodes))
    allocate(u_vec(n_nodes))
    allocate(K_local(n_edges,n_edges))


    subroutine generate_mesh(nx,ny,x_max,x_min,y_max,y_min,nodes,elems)


        dx = (x_max - x_min) / dble(nx)
        dy = (y_max - y_min) / dble(ny)

        node_id = 0
        do j = 0, ny
            do i = 0, nx
                node_id = node_id +1
                nodes(node_id)%x = x_min + dble(i) * dx
                nodes(node_id)%y = y_min + dble(i) * dy
            end do
        end do

