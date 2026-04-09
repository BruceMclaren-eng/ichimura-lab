program crs_practice
    implicit none
    
    integer :: n
    double precision, allocatable :: A(:,:)  
    real(8), allocatable :: values(:)
    integer, allocatable :: col_idx(:)
    integer, allocatable :: row_ptr(:)
    
    n=5
    allocate(A(n,n))
    A = 0.0d0
    A(1,1) = 1.0d0
    A(2,2) = 2.0d0
    A(3,3) = 3.0d0
    A(4,4) = 4.0d0
    A(5,5) = 5.0d0

    call dense_to_crs(A,n,values,col_idx,row_ptr)
    print *, "Values: ", values
    print *, "Column Indices: ", col_idx
    print *, "Row Pointers: ", row_ptr


contains

    subroutine dense_to_crs(A,n,values,col_idx,row_ptr)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: A(n,n)
        real(8), allocatable, intent(out) :: values(:)
        integer, allocatable, intent(out) :: col_idx(:)
        integer, allocatable, intent(out) :: row_ptr(:)
        integer :: i,j, nnz

        nnz = 0
        do i = 1,n
            do j = 1,n
                if(A(i,j) /= 0.0d0) then
                    nnz = nnz + 1
                end if
            end do
        end do

        allocate(values(nnz))
        allocate(col_idx(nnz))
        allocate(row_ptr(n + 1))

        nnz = 0
        row_ptr(1) = 1
        do i = 1, n
            do j = 1, n
                if(A(i,j) /= 0.0d0 ) then
                    nnz = nnz + 1
                    values(nnz) = A(i,j)
                    col_idx(nnz) = j
                end if
            end do
            row_ptr(i + 1) = nnz
        end do 

    end subroutine
end program crs_practice