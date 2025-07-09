
program test_eval_RotheError
   use, intrinsic :: iso_fortran_env, only : dp => real64
   use RotheSolver_mod, only : RotheSolver   
   use ho_kernels      , only : S, H, H2        
   use optimizers
   use stdlib_io_npy, only: save_npy
   implicit none

   abstract interface
      function matf(p) result(mat)
         import dp
         real(dp), intent(in)  :: p(:,:)
         complex(dp), allocatable :: mat(:,:)
      end function matf
      function matft(p, t) result(mat)
         import dp
         real(dp), intent(in)  :: p(:,:), t
         complex(dp), allocatable :: mat(:,:)
      end function matft
   end interface


   procedure(matf), pointer :: Sp  => S   
   procedure(matft), pointer :: Hp  => H
   procedure(matft), pointer :: H2p => H2

   integer, parameter :: n = 4, Nb = 4
   integer :: i, n_steps,niter
   real(dp)            :: p_old(n,Nb), p_new(n,Nb), a,b,mu,q
   complex(dp)         :: c_old(n)
   complex(dp), allocatable :: c_new(:)
   complex(dp), allocatable :: c_out(:,:)
   real(dp), allocatable :: p_out(:,:,:)
   real(dp), allocatable :: errors(:), times(:)
   complex(dp), allocatable :: gradient(:,:)
   type(RotheSolver)  :: solver
   real(dp) :: dt = 0.01_dp, t0 = 0.0_dp, err, err_old

   a=1.0_dp/sqrt(3.0_dp)
   b=0.1_dp
   mu=1.0_dp
   q= 0.0_dp
   p_old = 0.0_dp
   p_old(:,1) = a  
   p_old(:,2) = b  
   p_old(1,3) = mu                 
   p_old(2,3) = -mu                
   p_old(3,3) = 2.0_dp*mu
   p_old(4,3) = -2.0_dp*mu            
   p_old(:,4) = q       
   p_new = 0_dp
   do i=1, n
      print *, "p_old(", i, ") = ", p_old(i,:)
   end do
   p_new = p_old + p_new  ! add perturbation to p_old
   c_old(:) = 1_dp  ! initial coefficients
   allocate(c_new(n))
   call solver%init(dt, "HOscillator_1D", Sp, Hp, H2p)
   print *, "Initiated solver"
   err_old = solver%calc_Rothe_error(t0, p_old, c_old, p_old, c_new)
   err = solver%calc_Rothe_error(t0, p_old, c_old, p_new, c_new)
   print '(A,F12.6)', "Initial Rothe error = ", err
   gradient = solver%calc_Rothe_gradient(t0, p_old, c_old, p_new) 
   print *, "Gradient:"
   do i=1,n
      print *, "gradient(", i, ") = ", gradient(i,:)
   end do
   
   call optimize_gradientDescent(solver, t0, p_old, c_old,           &
                                 p_old,p_new, err,c_new, tol=1e-14_dp, max_iter=500, &
                                 learning_rate=1e-0_dp,niter=niter, verbose=.false.)
   print '(A,F12.6)', "Optimised Rothe error = ", err
   print *, "Optimised parameters:"
   do i=1, n
      print *, "p_old(", i, ") = ", p_new(i,:)
   end do
   n_steps=1000
   allocate(p_out(n,Nb,n_steps+1), c_out(n,n_steps+1), errors(n_steps+1),times(n_steps+1))
   call propagate_Nsteps(solver, t0, p_old, c_old,n_steps, p_out,c_out, errors,times, tol=1e-14_dp, &
                                    max_iter=500, learning_rate=1e0_dp, verbose=.false.)
   call save_npy("output/nonlinear_data.npy", p_out)
   caLL save_npy("output/coefficients.npy", c_out)
   call save_npy("output/errors.npy", errors)
   call save_npy("output/times.npy", times)
   print *, "Final error at step ", n_steps, " = ", sum(sqrt(errors))
end program test_eval_RotheError
