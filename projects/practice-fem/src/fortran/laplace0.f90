program laplace
!
! Copyright 2013 SofTek Kato
!
implicit none
integer(4) :: nn, nm
integer(4) :: n, m, iter, iter_max, i, j, k

parameter(nn=4096)
parameter(nm=4096)

real (kind=8) :: tol   = 1.0d-6 
real (kind=8) :: error = 1.0d0
real (kind=8), dimension (:,:) :: A (nn,nm)
real (kind=8), dimension (:,:) :: Anew(nn,nm)

! Elapsed time function initialized -- SofTek
integer ic, icr, icm
real*8 cpu0, cpu1, elapsed, t_ac
real*8   second
external second

call system_clock(ic,icr,icm)
t_ac= 1.0/real(icr,8)
print '(1x,a,e10.5)','Time measurement accuracy : ',t_ac
! 
n = nn
m = nm
iter_max = 1000

	A    = 0.0d0
	Anew = 0.0d0

	do j = 1, n
	  A(0, j)    = 1.0d0
	  Anew(0, j) = 1.0d0
	end do

    print *, "Jacobi relaxation Calculation:", n," x ", m, "mesh"

    cpu0 = second()
    iter = 0
    do while ( iter .le. iter_max-1 .and. error .gt. tol )

      error =0.d0
	  
	  do j = 2, m-1
		do i = 2, n-1
		  Anew(i,j) = 0.25 * ( A(i,j+1) + A(i,j-1) &
                             + A(i-1,j) + A(i+1,j) )	  
          error = max(error, abs(Anew(i,j) - A(i,j)))
        end do
	  end do 

	  do j = 2, m-1
		do i = 2, n-1
		  A(i,j) = Anew(i,j) 
        end do
	  end do 

      if ( mod (iter,100) == 0 ) print *, iter, error
      iter = iter + 1      

    end do 
    cpu1 = second()

! Printout Elapsed time (SofTek)

    elapsed = (cpu1 -cpu0) * t_ac
    print *, "total : ",elapsed, "sec"

end program

function second() result(rtime)
implicit none
!
      integer :: ic,ir,im
      real(8) :: rtime
!
      call system_clock(ic,ir,im)
!
      rtime= real(ic,8)
!
      return
end function