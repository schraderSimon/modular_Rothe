module optimizers
   use, intrinsic :: iso_fortran_env, only: dp => real64
contains
subroutine optimize_gradientDescent(solver , time, p_old , c_old,          &
                                    p_start,p_opt  , err_opt, c_new, tol,               &
                                    max_iter, learning_rate, niter,verbose)
   use, intrinsic :: iso_fortran_env , only : dp => real64
   use RotheSolver_mod               , only : RotheSolver
   implicit none
   class(RotheSolver), intent(inout) :: solver        ! already initialised
   complex(dp)         , intent(in)    :: time             ! start time
   real(dp)         , intent(in)    :: p_old(:,:)    ! (n,Nb) initial guess
   real(dp)         , intent(in)    :: p_start(:,:)  ! (n,Nb) initial guess for next step
   complex(dp)      , intent(in)    :: c_old(:)       ! previous coeffs
   real(dp)         , intent(out)   :: p_opt(:,:)     ! best parameters
   complex(dp)      , intent(out)   :: c_new(:)       ! coefficients at end
   real(dp)         , intent(out)   :: err_opt        ! final error
   integer, intent(out) :: niter
   real(dp), intent(in)   :: tol          
   integer , intent(in)   :: max_iter       ! default 200
   real(dp), intent(in)   :: learning_rate     
   logical , intent(in)   :: verbose       


   integer                :: n_par1, n_par2
   real(dp)               :: err_now, err_new, gnorm2, lr
   real(dp), allocatable  :: p_now(:,:), p_new(:,:), g_now(:,:), g_new(:,:)
   logical :: talk
   n_par1 = size(p_old,1);  n_par2 = size(p_old,2)
   allocate(p_now(n_par1,n_par2), p_new(n_par1,n_par2), g_now(n_par1,n_par2),g_new(n_par1,n_par2))
   lr=learning_rate
   p_now  = p_start
   err_now = solver%calc_Rothe_error(time, p_old, c_old, p_now, c_new)
   g_now = real( solver%calc_Rothe_gradient(time, p_old, c_old, p_now) )
   talk = verbose

   do niter = 1, max_iter
      p_new = p_now - lr*g_now
      err_new= solver%calc_Rothe_error(time, p_old, c_old, p_new, c_new)
      
      if (err_new>err_now) then
         if (talk) write(*,'("Step rejected, err_new > err_now")')
         ! step was worse, so reduce learning rate
         lr = lr * 0.8_dp
         if (lr < 1.0e-8_dp) then
            write(*,'("Learning rate too small, stopping")')
            exit
         end if
         cycle  ! skip to next iteration
      end if
      err_now = err_new
      g_new = real( solver%calc_Rothe_gradient(time, p_old, c_old, p_new) )
      gnorm2 = sum( g_new**2 )
      p_now  = p_new
      g_now=g_new
      if (talk) then
         write(*,'(I4,2ES14.5)') niter, err_now, gnorm2
      end if
      if (err_now < tol) exit  ! small error
      
   end do

   p_opt  = p_now
   err_opt = err_now
   deallocate(p_now, p_new, g_now, g_new)
   if (talk) write(*,'("Finished after",I4," iterations; err =",ES14.5)') niter, err_opt
end subroutine optimize_gradientDescent
subroutine optimize_LBFGS (solver , time, p_old , c_old,           &
                           p_start,p_opt ,  err_opt, c_new, tol,   &
                           max_iter, learning_rate, niter, verbose)
   use, intrinsic :: iso_fortran_env , only : dp => real64
   use RotheSolver_mod               , only : RotheSolver
   implicit none
   class(RotheSolver) , intent(inout) :: solver
   complex(dp)           , intent(in)    :: time
   real(dp)           , intent(in)    :: p_old(:,:)     
   complex(dp)        , intent(in)    :: c_old(:)
   real(dp)           , intent(in)    :: p_start(:,:)    ! initial guess
   real(dp)           , intent(out)   :: p_opt(:,:)      ! final parameters
   complex(dp)        , intent(out)   :: c_new(:)        ! coeffs at optimum
   real(dp)           , intent(out)   :: err_opt         ! f(x*)
   real(dp)           , intent(in)    :: tol             ! <— unused, exists for consistency
   integer            , intent(in) , optional :: max_iter
   real(dp)           , intent(in) , optional :: learning_rate  ! <— unused, exists for consistency
   logical            , intent(in) :: verbose
   integer            , intent(out) :: niter
   integer  :: n_par1, n_par2, n, m, iflag, maxit, iter
   logical  :: talk, diagco
   integer  :: iprint(2)
   real(dp) :: eps, xtol, gnorm
   real(dp) , allocatable :: x(:), g(:), diag(:), w(:)
   real(dp) , allocatable :: gmat(:,:), p_work(:,:)
   complex(dp), allocatable :: c_tmp(:)

   integer :: mp, lp !Stuff for lbfgs
   real(dp) :: gtol, stpmin, stpmax !Stuff for lbfs
   EXTERNAL LB2
   common /lb3/ mp, lp, gtol, stpmin, stpmax ! These are global parameters, i.e. LBFGS.f has access to them
   interface
      subroutine lbfgs(nv, mv, x_v, f, g_v, diagco_v, diag_v, iprint_v, gtolerance, x_tol, ww, iflag_lbfgs)
         use, intrinsic :: iso_fortran_env , only : dp => real64
         integer, intent(in)    :: nv, mv
         real(dp), intent(inout):: x_v(nv), f
         real(dp), intent(inout):: g_v(nv)
         logical , intent(in)   :: diagco_v
         real(dp), intent(inout):: diag_v(nv)
         integer, intent(in)    :: iprint_v(2)
         real(dp), intent(in)   :: gtolerance, x_tol
         real(dp), intent(inout):: ww(*)
         integer, intent(inout) :: iflag_lbfgs
      end subroutine lbfgs
   end interface
   n_par1 = size(p_old,1);   n_par2 = size(p_old,2)
   n      = n_par1 * n_par2
   m      = 30                       ! Can be set to very large value to restore full BFGS
   maxit  = merge(max_iter, 1000, present(max_iter))
   talk   = verbose

   iprint = -1                       ! LBFGS print control. -1 means no printing
   if (talk) iprint = [1,0]       ! print every iteration

   xtol = 3e-8_dp   !For finite difference gradient, this is a good value
   eps = 3e-8_dp  !For finite difference gradient, this is a good value
   diagco = .true. ! User uses an initial diagonal scaling
  
   allocate(x(n), g(n), diag(n), w(n*(2*m+1)+2*m))
   allocate(gmat(n_par1,n_par2), p_work(n_par1,n_par2))
   allocate(c_tmp(size(c_old)))

   x = reshape(p_start, [n]) ! p_start is (n_par1,n_par2) matrix, x is now a vector of length n

   gtol   = 0.9_dp ! Smaller value -> Better line search, but more function evaluations. Should be between 0.1 and 0.9
   stpmin = 1.0e-20_dp
   stpmax = 1.0e+20_dp


   iflag = 0
   iter  = 0
   p_work = reshape(x, shape(p_work))
   err_opt = solver%calc_Rothe_error(time, p_old, c_old, p_work, c_tmp)
   gmat    = real( solver%calc_Rothe_gradient(time, p_old, c_old, p_work) )
   g= reshape(gmat, [n]) ! Flatten gradient into g(:)

   gnorm= sqrt(sum(g**2)) ! initial gradient norm
   print *, "Initial gradient norm = ", gnorm
   diag(1:n) = 1.0_dp/(abs(g)+1e-10_dp) ! initial diagonal scaling, avoid division by zero

   do
      p_work = reshape(x, shape(p_work))

      err_opt = solver%calc_Rothe_error(time, p_old, c_old, p_work, c_tmp)
      gmat    = real( solver%calc_Rothe_gradient(time, p_old, c_old, p_work) )
      g       = reshape(gmat, [n])

      call lbfgs(n, m, x, err_opt, g, diagco, diag, iprint, eps, xtol, w, iflag)
      diagco= .false. ! reset diagco to false for next iteration
      iter = iter + 1

      if (iflag <= 0) exit          ! < 0 : error, 0 : converged
      if (iter >= maxit) then
         if (talk) write(*,'("LBFGS – reached max_iter (",I0,")")') maxit
         exit
      end if
   end do
   niter = iter
   p_opt = reshape(x, shape(p_opt))
   c_new = c_tmp

   if (talk) &
      write(*,'("LBFGS finished after ",I0," iterations,  f = ",ES12.5)') &
            niter, err_opt

   deallocate(x, g, diag, w, gmat, p_work, c_tmp)
end subroutine optimize_LBFGS

subroutine propagate_Nsteps(solver, t0, p_init, c_init,n_steps, p_out,c_out, errors,times, tol, &
                                    max_iter, learning_rate, verbose)
   use, intrinsic :: iso_fortran_env , only : dp => real64
   use RotheSolver_mod               , only : RotheSolver
   implicit none
   class(RotheSolver), intent(inout) :: solver        
   complex(dp)         , intent(in)    :: t0             ! start time
   real(dp)         , intent(in)    :: p_init(:,:)    ! (n,Nb) initial guess
   complex(dp)      , intent(in)    :: c_init(:)       ! previous coeffs
   integer          , intent(in)    :: n_steps    ! number of steps to propagate
   real(dp)         , intent(in)    :: tol            ! default 1e-8
   real(dp)         , allocatable,intent(out)   :: p_out(:,:,:)     ! best parameters at each step
   complex(dp)      , allocatable, intent(out)   :: c_out(:,:)       ! coefficients at each step
   real(dp)         ,allocatable, intent(out)   :: errors(:)     ! errors at each step
   real(dp), allocatable, intent(out) :: times(:) ! times at each step
   
   integer , intent(in), optional   :: max_iter       
   real(dp), intent(in), optional   :: learning_rate       
   logical , intent(in), optional   :: verbose       
   integer                :: n_par1, n_par2, iteration, niter, max_iterations
   real(dp)               :: err_new, dx, lr
   complex(dp) :: time
   logical :: talk
   real(dp), allocatable  :: p_old(:,:), p_new(:,:), p_startguess(:,:)
   complex(dp), allocatable :: c_new(:), c_old(:)
   if (present(max_iter)) then
      max_iterations = max_iter
   else
      max_iterations = 200
   end if
   if (present(learning_rate)) then
      lr=learning_rate
   else
      lr=1
   end if
   if (present(verbose)) then
      talk=verbose
   else
      talk=.false.
   end if 
   n_par1 = size(p_init,1)
   n_par2 = size(p_init,2)
   allocate(p_old(n_par1, n_par2), p_new(n_par1, n_par2), &
            p_startguess(n_par1, n_par2))
   allocate(c_new(size(c_init)))
   allocate(c_old(size(c_init)))
   allocate(p_out(n_par1, n_par2, n_steps+1))
   allocate(c_out(size(c_init), n_steps+1))
   allocate(errors(n_steps+1),times(n_steps+1))
   time = t0
   p_old = p_init
   c_old = c_init
   p_out(:,:,1) = p_old
   c_out(:,1) = c_old
   errors(1) = 0_dp  ! Initial error is zero (Could in principle be the variance of the ground state for initial eigenstates)
   times(1) = time  ! Initial time
   p_startguess = p_init
   dx=1.0_dp
   do iteration = 1, n_steps
      call optimize_LBFGS(solver, time, p_old, c_old, p_startguess,p_new,err_new, c_new, tol, &
                                    max_iterations, lr, niter,talk)
      p_old = p_new
      c_old = c_new
      p_startguess = (1.0_dp+dx)*p_new - dx*p_out(:,:,iteration)  ! use previous change as initial guess
      c_out(:,iteration+1) = c_new
      p_out(:,:,iteration+1) = p_new
      print *, c_new
      errors(iteration+1) = err_new
      times(iteration+1) = time
      time =time+ solver%dt  ! increment time by dt
      print *, "Error at step ", iteration, " = ", err_new, "niter: ", niter, "Time: ", time
      
   end do

end subroutine propagate_Nsteps
end module optimizers