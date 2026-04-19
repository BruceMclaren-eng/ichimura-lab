module matvec_cuda_module
    implicit none
    interface
        subroutine matvec_cuda_device(n, values, col_idx, row_ptr, x, y) &
            bind(C, name="matvec_cuda_device_")
            use iso_c_binding
            integer(c_int), intent(in)    :: n
            real(c_double), intent(in)    :: values(*)
            integer(c_int), intent(in)    :: col_idx(*)
            integer(c_int), intent(in)    :: row_ptr(*)
            real(c_double), intent(in)    :: x(*)
            real(c_double), intent(inout) :: y(*)
        end subroutine
    end interface
end module matvec_cuda_module