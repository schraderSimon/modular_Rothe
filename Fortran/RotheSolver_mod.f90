module RotheSolver_mod
   use, intrinsic :: iso_fortran_env, only : dp => real64
   implicit none
   private

   real(dp), parameter :: FD_eps = 1.0e-6_dp   ! finite-difference step for gradient
   complex(dp), parameter :: Im =(0.0_dp,1.0_dp)
   abstract interface
      function matrix_func(p) result(Smat)
         import dp
         real(dp), intent(in)  :: p(:,:)
         complex(dp), allocatable :: Smat(:,:)
      end function matrix_func
      function matrix_func_t(p, t) result(mat)
         import dp
         real(dp), intent(in)  :: p(:,:), t
         complex(dp), allocatable :: mat(:,:)
      end function matrix_func_t
   end interface

   public :: RotheSolver
   type :: RotheSolver
      real(dp),  allocatable :: p(:,:)          ! (n, Nb)  nonlinear
      complex(dp), allocatable :: c(:)          ! (n)      linear
      real(dp) :: dt
      character(len=:), allocatable :: name
      procedure(matrix_func), pointer, nopass :: S  => null()
      procedure(matrix_func_t), pointer, nopass :: H  => null()
      procedure(matrix_func_t), pointer, nopass :: H2 => null()
      contains
      procedure :: init
      procedure :: calc_Rothe_error   => RS_error
      procedure :: calc_Rothe_gradient => RS_gradient_numerical
      procedure :: timeEvolve => RS_timeEvolve
   end type RotheSolver
contains
subroutine init(self, dt, name,S_cb, H_cb, H2_cb)
   class(RotheSolver), intent(inout) :: self
   real(dp),          intent(in)    :: dt
   character(*),      intent(in)    :: name
   procedure(matrix_func) , pointer       :: S_cb
   procedure(matrix_func_t) , pointer       :: H_cb, H2_cb

   self%dt   = dt
   self%name = name
   self%S  => S_cb
   self%H  => H_cb
   self%H2 => H2_cb
end subroutine init

function RS_error(self, t, p_old,c_old,p_new,c_new) result(err)
   class(RotheSolver), intent(in)    :: self
   real(dp),           intent(in)    :: t
   real(dp),           intent(in)    :: p_new(:,:), p_old(:,:)     ! (n,Nb)
   complex(dp),        intent(in)    :: c_old(:)                    ! (n)
   complex(dp), intent(out)   :: c_new (:)   ! (n) and reshaped to (n,1)
   real(dp)             :: err

   integer :: n,Nb, ngo, ngn, k, info,j,lwork
   real(dp), allocatable :: pcat(:,:), p_plus(:,:), p_minus(:,:)
   integer, allocatable :: ipiv(:)
   complex(dp), allocatable :: work(:)
   complex(dp), allocatable :: Sfull(:,:), Hfull(:,:), H2full(:,:)
   complex(dp), allocatable :: S_tilde_full(:,:), rho_mat(:,:), rho_vec(:)
   complex(dp), allocatable ::  Sigma_nn(:,:), &
                                Sigma_oo(:,:), eye(:,:), c_new_reshaped(:,:)
   real(dp), parameter :: lbmd = 1.0e-16_dp ! This should really be an argument
   complex(dp) :: overlap, proj

   n  = size(p_old,1);   Nb = size(p_old,2)
   ngo = n
   ngn = size(p_new,1)

   ! concatenate
   allocate(pcat(ngo+ngn,Nb))
   pcat(1:ngo,1:Nb)      = p_old
   pcat(ngo+1:ngo+ngn,1:Nb) = p_new
   allocate(Sfull(ngo+ngn, ngo+ngn), &
            Hfull(ngo+ngn, ngo+ngn), &
            H2full(ngo+ngn, ngo+ngn), &
            S_tilde_full(ngo+ngn, ngo+ngn), &
            rho_mat(ngo+ngn, ngo+ngn), &
            rho_vec(ngn))
   Sfull  = self%S (pcat)
   Hfull  = self%H (pcat,t+self%dt*0.5_dp)  ! t+dt/2
   H2full = self%H2(pcat,t+self%dt*0.5_dp)

   S_tilde_full      = Sfull + (self%dt**2/4.0_dp)*H2full
   allocate(Sigma_nn(ngn, ngn), Sigma_oo(ngo, ngo))

   Sigma_nn = S_tilde_full(ngo+1:, ngo+1:)
   Sigma_oo = S_tilde_full(:ngo,   :ngo   )
   rho_mat        = Sfull - (self%dt**2/4.0_dp)*H2full + Im*self%dt*Hfull
   rho_vec      = matmul(conjg(transpose(rho_mat(:ngo, ngo+1:))), c_old)

   allocate(eye(ngn,ngn)); eye = 0.0_dp
   do k= 1, ngn
      eye(k,k) = 1.0_dp
   end do

   Sigma_nn = Sigma_nn + lbmd*eye
   Sigma_nn = 0.5_dp*(Sigma_nn + conjg(transpose(Sigma_nn)))   ! enforce Hermiticity
   lwork = -1
   allocate(ipiv(ngn))
   allocate(work(1))
   allocate(c_new_reshaped(ngn,1))  ! ensure column vector
   c_new_reshaped(:,1)= rho_vec  ! ensure column vector

   call zhesv('U', ngn, 1, Sigma_nn, ngn, ipiv, c_new_reshaped, ngn, work, lwork, info)
   lwork = int(real(work(1)))          ! optimal size returned in WORK(1)
   deallocate(work);  allocate(work(lwork))
   call zhesv('U', ngn, 1, Sigma_nn, ngn, ipiv, c_new_reshaped, ngn, work, lwork, info)
   if (info /= 0) stop 'zhesv failed'

   overlap = dot_product((c_old), matmul(Sigma_oo, c_old))
   c_new   = c_new_reshaped(:,1)  ! restore original shape
   proj    = dot_product((rho_vec),    c_new)
   err     = real(overlap - proj)
end function RS_error
function RS_gradient_numerical(self, t, p_old, c_old, p_new) result(grad)
   class(RotheSolver), intent(in)    :: self
   real(dp),           intent(in)    :: t
   real(dp),           intent(in)    :: p_new(:,:), p_old(:,:)     ! (n,Nb)
   complex(dp),        intent(in)    :: c_old(:)                    ! (n)
   complex(dp), allocatable   :: c_new (:)
   real(dp), allocatable :: grad(:,:)

   integer :: n, Nb, k,j
   real(dp), allocatable :: p_trial_plus(:,:), p_trial_minus(:,:)
   real(dp) :: step,  rplus,rminus

   n  = size(p_old,1);   Nb = size(p_old,2)

   allocate(grad(n,Nb))
   allocate(c_new(n))
   allocate(p_trial_plus(n,Nb))
   allocate(p_trial_minus(n,Nb))
   step = FD_eps
   do k=1,Nb
      do j=1, n
         p_trial_plus = p_new
         p_trial_minus = p_new
         p_trial_plus(j,k) = p_trial_plus(j,k) + step
         p_trial_minus(j,k) = p_trial_minus(j,k) - step
         rplus= RS_error(self, t, p_old, c_old, p_trial_plus, c_new)
         rminus= RS_error(self, t, p_old, c_old, p_trial_minus, c_new)
         grad(j,k) = (rplus-rminus)/(2*step)
      end do
   end do
   deallocate(p_trial_plus, p_trial_minus)
end function RS_gradient_numerical
!─────────────────────────────────────────────────────────────────────────────
subroutine RS_timeEvolve(self)
    class(RotheSolver), intent(in)    :: self
end subroutine RS_timeEvolve
end module RotheSolver_mod
