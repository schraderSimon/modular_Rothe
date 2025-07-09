program ho_chain
   !
   ! 1-D harmonic-oscillator spectrum in a Gaussian μ-grid basis
   !
   use, intrinsic :: iso_fortran_env, only : dp => real64
   use ho_kernels                   ! S , x , x2 , H2 H already here
   implicit none
   !–––– grid definition ––––
   real(dp), parameter :: width = 5_dp          !  -5 … +5
   integer,  parameter :: n  = 21                 !  -5 … +5 in steps of 0.1
   real(dp), parameter :: a  = 1.0_dp/sqrt(2.0_dp) !  ground-state width, ω=1
   real(dp)           :: mus(n), as(n),bs(n),ps(n)  ! μ-grid, a, b, p
   real(dp), allocatable :: p(:,:)                 ! (n,4) parameter table
   !–––– matrices & eigen stuff ––––
   complex(dp), allocatable :: H2mat(:,:), Hmat(:,:), Smat(:,:), Amat(:,:), Bmat(:,:)   ! full (n,n)
   real(dp),    allocatable :: eval(:) , eval2(:)               ! eigenvalues
   !––– LAPACK work arrays (will be allocated after a size query) –––
   complex(dp), allocatable :: work(:)
   real(dp),    allocatable :: rwork(:)
   integer                    :: lwork, info
   integer                    :: i
   !–––– fill the μ grid ––––
   real(dp) :: delta
    delta = 2.0_dp*width / (n-1)   ! 10.0 = 5.0 - (-5.0)
   do i = 1, n
      mus(i) = -width + delta*(i-1)
      as(i) = a + (rand()-0.5_dp) * 0.1_dp
      bs(i) = 0_dp+(rand()-0.5_dp) * 1_dp  ! b-skew, 0.01 < b < 0.11
      ps(i) = 0_dp+(rand()-0.5_dp) * 1_dp                ! p = 0.0
   end do
   !–––– build parameter table p(:,1:4) = (a, b, μ, p) ––––
   allocate(p(n,4))
   p(:,1) =  as                 ! width      a
   p(:,2) =  bs            ! b-skew     b
   p(:,3) =  mus               ! centre     μ
   p(:,4) =  ps            ! momentum   p
! add a second eigenvalue array

! …

!–––– assemble matrices ––––
Hmat  = H (p)        ! ⟨g_i|Ĥ|g_j⟩
H2mat = H2(p)        ! ⟨g_i|Ĥ²|g_j⟩
Smat  = S (p)        ! ⟨g_i|g_j⟩

allocate(eval(n), eval2(n))

!========================== 1.  diagonalise Ĥ ==============================
do i=1, n
   print *, Hmat(i,:) ! 
end do
lwork  = -1
allocate(work(1), rwork(1))
call zhegv(1, 'N', 'U', n, Hmat, n, Smat, n, eval,  work, lwork, rwork, info)
lwork = int(real(work(1)))
deallocate(work, rwork)
allocate(work(lwork), rwork(3*n-2))

call zhegv(1, 'N', 'U', n, Hmat, n, Smat, n, eval,  work, lwork, rwork, info)
if (info /= 0) stop 'ZHEGV for H failed'

deallocate(work, rwork)         ! reuse workspace for the next solve

!========================== 2.  diagonalise Ĥ² =============================
Smat= S(p)        ! ⟨g_i|g_j⟩
do i=1, n
   print *, H2mat(i,:) ! 
end do
lwork  = -1
allocate(work(1), rwork(1))
call zhegv(1, 'N', 'U', n, H2mat, n, Smat, n, eval2, work, lwork, rwork, info)
lwork = int(real(work(1)))
deallocate(work, rwork)
allocate(work(lwork), rwork(3*n-2))

call zhegv(1, 'N', 'U', n, H2mat, n, Smat, n, eval2, work, lwork, rwork, info)
if (info /= 0) stop 'ZHEGV for H² failed'

!–––– print results ––––
print '(A)', 'Eigenvalues of Ĥ      :'
do i = 1, n
   print '(I3,2X,F12.8)', i-1, eval (i)
end do

print '(A)', 'Eigenvalues of Ĥ squared:'
do i = 1, n
   print '(I3,2X,F12.8)', i-1, eval2(i)
end do
end program ho_chain