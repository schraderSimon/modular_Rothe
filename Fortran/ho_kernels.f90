module ho_kernels
   use, intrinsic :: iso_fortran_env, only : dp => real64
   implicit none
   private
   public :: S,  H, H2, dS
   real(dp), parameter         :: pi = 3.141592653589793238462643383279502884197169399375_dp
   complex(dp), parameter      :: Im = (0.0_dp, 1.0_dp)
!! Parameters you don’t want to thread through every call
   real(dp), parameter :: omega = 1.0_dp

contains
!────────────── overlap S ──────────────
   pure subroutine moments(bA_i, bA_j, mu_i, mu_j, p_i, p_j, &
                           alpha, alpha2, alpha3,alpha4,alpha5,alpha6, beta, beta2,beta4)
      complex(dp), intent(in)  :: bA_i, bA_j
      real(dp),    intent(in)  :: mu_i, mu_j, p_i, p_j
      complex(dp), intent(out) :: alpha, alpha2, alpha3, beta, beta2, beta4, alpha4, alpha5, alpha6
      alpha  = conjg(bA_i) + bA_j
      alpha2=alpha*alpha
      alpha3=alpha2*alpha
      alpha2 = alpha*alpha
      alpha3 = alpha2*alpha
      alpha4 = alpha3*alpha
      alpha5 = alpha4*alpha
      alpha6 = alpha5*alpha
      beta  = 2.0_dp*(conjg(bA_i)*mu_i + bA_j*mu_j) + Im*(p_j - p_i)
      beta2  = beta*beta
      beta4  = beta2*beta2
   end subroutine moments

   function S(p) result(Smat)
      real(dp),    intent(in)     :: p(:,:)          ! (n,4), the parameters of the polynomial
      complex(dp), allocatable    :: Smat(:,:)       ! (n,n), the overlap matrix
      real(dp), allocatable       :: normalizations(:) ! (n), the normalization factors
      integer                     :: n, i, j!, Nparam
      complex(dp)                 :: N_i, N_j, alpha, beta, gamma, Z, bA_i, bA_j
      
      
      real(dp)                    :: a_i, a_j, b_i, b_j, mu_i, mu_j, p_i, p_j
      
      n = size(p,1)
      !Nparam= size(p,2) ! Number of parameters per Gaussian
      allocate(Smat(n,n))
      allocate (normalizations(n))
      do i = 1, n
         a_i = p(i,1) ! a_i
         b_i = p(i,2) ! b_i
         mu_i = p(i,3) ! mu_i
         p_i = p(i,4) ! p_i
         ! Compute normalization factor for each Gaussian
         normalizations(i) = (2.0_dp * a_i**2 / pi)**0.25_dp
      end do

      do j = 1, n
         a_j = p(j,1) ! a_j
         b_j = p(j,2) ! b_j
         mu_j = p(j,3) ! mu_j
         p_j = p(j,4) ! p_j
         N_j= cmplx(normalizations(j),0.0_dp, kind=dp) ! Normalization factor for Gaussian j
         bA_j = cmplx(a_j**2, b_j, kind=dp)

         do i = 1, j
            N_i= cmplx(normalizations(i),0.0_dp, kind=dp) ! Normalization factor for Gaussian i
            a_i = p(i,1) ! a_i
            b_i = p(i,2) ! b_i
            mu_i = p(i,3) ! mu_i
            p_i = p(i,4) ! p_i
            
            bA_i = cmplx(a_i**2, b_i, kind=dp)
            alpha = conjg(bA_i) + bA_j
            beta = 2.0_dp*(conjg(bA_i)*mu_i + bA_j*mu_j) + Im*(p_j - p_i)
            gamma = -(conjg(bA_i)*mu_i**2 + bA_j*mu_j**2) + Im*(p_i*mu_i - p_j*mu_j)
            Z = N_i * N_j * sqrt(pi/alpha) * exp(gamma + beta**2/(4.0_dp*alpha))
            Smat(i,j) = Z
            Smat(j,i) = conjg(Z) ! Ensure symmetry
         end do
      end do
      deallocate(normalizations)
   end function S
   function xn(p,order) result(mat)
      real(dp),    intent(in)  :: p(:,:)       ! (n,4)
      integer,     intent(in)  :: order        ! Order of the polynomial
      complex(dp), allocatable :: mat(:,:)   ! (n,n)
      integer :: n, i, j
      complex(dp) :: bA_i,bA_j, alpha,beta, fact
      complex(dp) :: alpha2,alpha3,alpha4,alpha5,alpha6, beta2,beta4
      real(dp) :: a_i,b_i,mu_i,p_i, a_j,b_j,mu_j,p_j

      n = size(p,1)
      allocate(mat(n,n))
      select case (order)
         case (1)
            do j=1,n
               a_j = p(j,1); b_j = p(j,2); mu_j = p(j,3); p_j = p(j,4)
               bA_j = cmplx(a_j**2,b_j,kind=dp)
               do i=1,j
                  a_i = p(i,1); b_i = p(i,2); mu_i = p(i,3); p_i = p(i,4)
                  bA_i = cmplx(a_i**2,b_i,kind=dp);      

                  call moments(bA_i, bA_j, mu_i, mu_j, p_i, p_j, &
                           alpha, alpha2, alpha3,alpha4,alpha5,alpha6, beta, beta2,beta4)

                  fact = beta / (2.0_dp*alpha)
                  mat(i,j) = fact
                  mat(j,i) = conjg(mat(i,j))
               end do
          end do
         case (2)
            do j=1,n
               a_j = p(j,1); b_j = p(j,2); mu_j = p(j,3); p_j = p(j,4)
               bA_j = cmplx(a_j**2,b_j,kind=dp)
               do i=1,j
                  a_i = p(i,1); b_i = p(i,2); mu_i = p(i,3); p_i = p(i,4)
                  bA_i = cmplx(a_i**2,b_i,kind=dp);      

                  call moments(bA_i, bA_j, mu_i, mu_j, p_i, p_j, &
                           alpha, alpha2, alpha3,alpha4,alpha5,alpha6, beta, beta2,beta4)

                  fact = beta2/(4.0_dp*alpha2) + 1.0_dp/(2.0_dp*alpha)
                  mat(i,j) =  fact
                  mat(j,i) = conjg(mat(i,j))
               end do
            end do
         case (3)
            do j=1,n
               a_j = p(j,1); b_j = p(j,2); mu_j = p(j,3); p_j = p(j,4)
               bA_j = cmplx(a_j**2,b_j,kind=dp)
               do i=1,j
                  a_i = p(i,1); b_i = p(i,2); mu_i = p(i,3); p_i = p(i,4)
                  bA_i = cmplx(a_i**2,b_i,kind=dp);      

                  call moments(bA_i, bA_j, mu_i, mu_j, p_i, p_j, &
                           alpha, alpha2, alpha3,alpha4,alpha5,alpha6, beta, beta2,beta4)

                  fact = 3.0_dp/4.0_dp * beta/(alpha2) + 1.0_dp/8.0_dp* (beta2*beta/alpha3)
                  mat(i,j) =  fact
                  mat(j,i) = conjg(mat(i,j))
               end do
            end do

         case (4)
            do j=1,n
               a_j = p(j,1); b_j = p(j,2); mu_j = p(j,3); p_j = p(j,4)
               bA_j = cmplx(a_j**2,b_j,kind=dp)
               do i=1,j
                  a_i = p(i,1); b_i = p(i,2); mu_i = p(i,3); p_i = p(i,4)
                  bA_i = cmplx(a_i**2,b_i,kind=dp);      

                  call moments(bA_i, bA_j, mu_i, mu_j, p_i, p_j, &
                           alpha, alpha2, alpha3,alpha4,alpha5,alpha6, beta, beta2,beta4)
                  fact = 3.0_dp/(4.0_dp*alpha2) + 3.0_dp*beta2/(4.0_dp*alpha3) + beta4/(16.0_dp*alpha4)
                  mat(i,j) = fact
                  mat(j,i) = conjg(mat(i,j))
               end do
            end do
         case (5)
            do j=1,n
               a_j = p(j,1); b_j = p(j,2); mu_j = p(j,3); p_j = p(j,4)
               bA_j = cmplx(a_j**2,b_j,kind=dp)
               do i=1,j
                  a_i = p(i,1); b_i = p(i,2); mu_i = p(i,3); p_i = p(i,4)
                  bA_i = cmplx(a_i**2,b_i,kind=dp);      

                  call moments(bA_i, bA_j, mu_i, mu_j, p_i, p_j, &
                           alpha, alpha2, alpha3,alpha4,alpha5,alpha6, beta, beta2,beta4)

                  fact = 15.0_dp*beta/(8.0_dp*alpha3) + 5.0_dp*beta2*beta/(8.0_dp*alpha4) + beta**5/(32.0_dp*alpha5)
                  mat(i,j) = fact
                  mat(j,i) = conjg(mat(i,j))
               end do
            end do
         case (6)
            do j=1,n
               a_j = p(j,1); b_j = p(j,2); mu_j = p(j,3); p_j = p(j,4)
               bA_j = cmplx(a_j**2,b_j,kind=dp)
               do i=1,j
                  a_i = p(i,1); b_i = p(i,2); mu_i = p(i,3); p_i = p(i,4)
                  bA_i = cmplx(a_i**2,b_i,kind=dp);      

                  call moments(bA_i, bA_j, mu_i, mu_j, p_i, p_j, &
                           alpha, alpha2, alpha3,alpha4,alpha5,alpha6, beta, beta2,beta4)

                  fact = 15.0_dp/(8.0_dp*alpha3) + 45.0_dp*beta2/(16.0_dp*alpha4) + 15.0_dp*beta4/(32.0_dp*alpha5) + beta2*beta4/(64.0_dp*alpha6)

                  mat(i,j) = fact
                  mat(j,i) = conjg(mat(i,j))
               end do
            end do
         case default 
            print *, "Error: Unsupported order for xn function. Supported orders are 1 to 6."
            stop
      end select
   end function xn

   function H(p,t) result(Hmat)
      use, intrinsic :: iso_fortran_env, only : dp => real64
      implicit none
      real(dp),    intent(in)  :: p(:,:)
      complex(dp), intent(in) :: t
      complex(dp), allocatable :: Hmat(:,:)
      !– helper matrices already coded –
      complex(dp), allocatable :: Smat(:,:), Xmat(:,:), X2mat(:,:)
      !– loop indices & scratch –
      integer           :: n, i, j
      real(dp)          :: a_j, b_j, mu_j, p_j
      complex(dp)       :: bA_j                      !  bA_j = a_j² + i b_j

      n = size(p,1)
      allocate(Smat(n,n), Xmat(n,n), X2mat(n,n), Hmat(n,n))

      Smat = S (p)      ! ⟨g_i|g_j⟩
      Xmat = xn (p,1)*Smat      ! ⟨g_i|x|g_j⟩
      X2mat = xn(p,2)*Smat     ! ⟨g_i|x²|g_j⟩

      do j = 1, n                      ! column → parameters of |g_j⟩
         a_j  = p(j,1);  b_j = p(j,2)
         mu_j = p(j,3);  p_j = p(j,4)
         bA_j  = cmplx(a_j**2, b_j, kind=dp)

         do i = 1, n                   ! row → ⟨g_i|
            !──────── kinetic part  –½⟨∂²⟩  (see derivation in analysis) ────────
            Hmat(i,j) =  -2.0_dp*bA_j**2 * ( X2mat(i,j) - 2.0_dp*mu_j*Xmat(i,j) + mu_j**2*Smat(i,j) )   &
                        + 2.0_dp*Im*bA_j*p_j * ( Xmat(i,j) - mu_j*Smat(i,j) )                          &
                        + ( 0.5_dp*p_j**2 + bA_j ) * Smat(i,j)                                         &
                        !──────── potential part  ½ ω² ⟨x²⟩ ────────
                        + 0.5_dp*omega**2 * X2mat(i,j)
         end do
      end do
      !– ensure symmetry –
      Hmat = 0.5_dp * (Hmat + conjg(transpose(Hmat)))
   deallocate(Smat, Xmat, X2mat)
   end function H


   function H2(p,t) result(H2mat)
      use, intrinsic :: iso_fortran_env, only : dp => real64
      implicit none
      real(dp),    intent(in)  :: p(:,:)
      complex(dp), intent(in) :: t
      complex(dp), allocatable :: H2mat(:,:)
      !– helper matrices already coded –
      complex(dp), allocatable :: Smat(:,:), Xmat(:,:), X2mat(:,:), X3mat(:,:), X4mat(:,:)
      !– loop indices & scratch –
      integer           :: n, i, j
      real(dp)          :: a_j, b_j, mu_j, p_j, a_i, b_i, mu_i, p_i
      complex(dp)       :: bA_j, bA_i                      !  bA_j = a_j² + i b_j
      complex(dp)       :: bA_i_c,c(0:10), in(0:39)
      n = size(p,1)
      allocate(Smat(n,n), Xmat(n,n), X2mat(n,n), X3mat(n,n), X4mat(n,n), H2mat(n,n))
      Smat = S (p)      ! ⟨g_i|g_j⟩
      Xmat = xn (p,1)!*Smat      ! ⟨g_i|x|g_j⟩
      X2mat = xn(p,2)!*Smat     ! ⟨g_i|x²|g_j⟩
      X3mat = xn(p,3)!*Smat
      X4mat = xn(p,4)!*Smat
      do j= 1, n
         a_j  = p(j,1);  b_j = p(j,2)
         mu_j = p(j,3);  p_j = p(j,4)
         bA_j  = cmplx(a_j**2, b_j, kind=dp)
         do i = 1,n
            a_i  = p(i,1);  b_i = p(i,2)
            mu_i = p(i,3);  p_i = p(i,4)
            bA_i  = cmplx(a_i**2, b_i, kind=dp)
            bA_i_c= conjg(bA_i)
            ! Kinetic squared
            ! This uses the fact that <gi|∂⁴|gj> = <∂²gi|∂²gj> (i.e. it is relatively simple to compute) 
            in(0) = bA_i_c*bA_j
            in(1) = p_j**2
            in(2) = bA_i_c*in(1)
            in(3) = p_i**2
            in(4) = bA_j*in(3)
            in(5) = Im*p_i
            in(6) = 8.0_dp*in(0)
            in(7) = in(5)*in(6)
            in(8) = Im*p_j
            in(9) = 16.0_dp*in(0)*p_i*p_j
            in(10) = in(9)*mu_i
            in(11) = mu_j**2
            in(12) = bA_j**2
            in(13) = 8.0_dp*bA_i_c*in(12)
            in(14) = mu_i**2
            in(15) = bA_i_c**2
            in(16) = 8.0_dp*bA_j*in(15)
            in(17) = 4.0_dp*in(2)*in(5)
            in(18) = 4.0_dp*in(1)*in(15)
            in(19) = 4.0_dp*in(12)*in(3)
            in(20) = 16.0_dp*in(12)
            in(21) = bA_i_c*in(20)*in(5)
            in(22) = 16.0_dp*bA_j*in(15)*in(8)
            in(23) = in(14)*in(22)
            in(24) = in(15)*in(20)
            in(25) = in(14)*in(24)
            in(26) = in(15)*mu_i
            in(27) = in(26)*mu_j
            in(28) = 32.0_dp*bA_j*in(8)
            in(29) = 32.0_dp*in(12)
            in(30) = in(26)*in(29)
            in(31) = in(29)*mu_j
            in(32) = in(15)*in(31)
            in(33) = in(21)*mu_i
            in(34) = in(8)*mu_j
            in(35) = in(12)*mu_j
            in(36) = 16.0_dp*bA_i_c
            in(37) = in(15)*mu_i
            in(38) = Im*bA_j*p_j
            in(39) = Im*p_i


            c(0) = 4.0_dp*in(0) + in(1)*in(3) + in(10)*mu_j - in(11)*in(13) - in(11)*in(19) + in(11)*in(25) - in(11)*in(33) &
                   - in(14)*in(16) - in(14)*in(18) + in(17)*mu_i + 2.0_dp*in(2) + in(23)*mu_j - 4.0_dp*in(34)*in(4) - in(34)*in(6) + 2.0_dp*in(4) + in(7)*mu_i
            c(1) = 32.0_dp*bA_i_c*in(35)*in(39)*mu_i + 8.0_dp*bA_i_c*in(38) + 16.0_dp*bA_j*in(37) + 8.0_dp*in(1)*in(37) - in(10) &
                   + in(11)*in(12)*in(36)*in(39) - in(11)*in(30) - in(14)*in(32) - in(17) - in(23) - in(27)*in(28) + 8.0_dp*in(3)*in(35) + 4.0_dp*in(3)*in(38) + in(35)*in(36) - in(7) - in(9)*mu_j
            c(2) = -bA_i_c*in(31)*in(5) + in(11)*in(24) + 64.0_dp*in(12)*in(27) - in(13) - in(16) - in(18) - in(19) + in(22)*mu_j + in(25) + in(26)*in(28) - in(33) + in(9)
            c(3) = 16.0_dp*Im*bA_i_c*in(12)*p_i - in(22) - in(30) - in(32)
            c(4) = in(24)


            H2mat(i,j) = 0.25_dp * (c(0) + c(1)*Xmat(i,j) + c(2)*X2mat(i,j) + c(3)*X3mat(i,j) + c(4)*X4mat(i,j))

            ! Potential squared
            H2mat(i,j) = H2mat(i,j) + 0.25_dp * omega**4 *X4mat(i,j) !Very easy
            ! Kinetic times potential
            c(5) = 4_dp*bA_i_c**2*mu_i**2 - 4_dp*Im*bA_i_c*mu_i*p_i - 2_dp*bA_i_c - p_i**2
            c(6) = 4_dp*bA_i_c*(-2_dp*bA_i_c*mu_i + Im*p_i)
            c(7) = 4_dp*bA_i_c**2
            H2mat(i,j) = H2mat(i,j) - 0.25_dp * (c(5)*X2mat(i,j) + c(6)*X3mat(i,j) + c(7)*X4mat(i,j)) * omega**2 
            ! Potential times kinetic
            c(8) = 4_dp*bA_j**2*mu_j**2 + 4_dp*Im*bA_j*mu_j*p_j - 2_dp*bA_j - p_j**2
            c(9) = 4_dp*bA_j*(-2.0_dp*bA_j*mu_j - Im*p_j)
            c(10) = 4_dp*bA_j**2
            H2mat(i,j) = H2mat(i,j) - 0.25_dp * (c(8)*X2mat(i,j) + c(9)*X3mat(i,j) + c(10)*X4mat(i,j)) * omega**2 
         end do
      end do
      H2mat=H2mat*Smat ! Don't forget that all matrix elements need to be multiplied with the overlap!
      deallocate(Smat, Xmat, X2mat, X3mat, X4mat)
   end function H2


   function dS(p,Smat_in) result(dSmat)
      use, intrinsic :: iso_fortran_env, only : dp => real64
      implicit none
      real(dp),    intent(in)  :: p(:,:)              ! (n,4)
      complex(dp), intent(in), optional, target :: Smat_in(:,:)      ! (n,n)
      complex(dp), pointer          :: Smat(:,:)
      complex(dp), allocatable :: dSmat(:,:,:)        ! (n,n,4)

      integer :: n, i, j
      real(dp) :: a_i,b_i,mu_i,p_i,  a_j,b_j,mu_j,p_j
      complex(dp) :: bA_i,bA_j, alpha,beta, Sij, dlog(4)
      n = size(p,1)
      if (present(Smat_in)) then
         Smat => Smat_in             ! alloc/copy handled automatically
      else
         allocate(Smat(n,n))
         Smat = S(p)
      end if
      allocate(dSmat(n,n,4))
      do i = 1, n
         a_i = p(i,1);  b_i = p(i,2);  mu_i = p(i,3);  p_i = p(i,4)
         bA_i = cmplx(a_i**2, b_i, kind=dp)

         do j = 1, n
            a_j = p(j,1);  b_j = p(j,2);  mu_j = p(j,3);  p_j = p(j,4)
            bA_j = cmplx(a_j**2, b_j, kind=dp)

            alpha = conjg(bA_i) + bA_j
            beta  = 2.0_dp*(conjg(bA_i)*mu_i + bA_j*mu_j) + Im*(p_j - p_i)

            Sij = Smat(i,j)

            dlog(1) =  0.5_dp/a_j - a_j/alpha + 2.0_dp*a_j*mu_j*beta/alpha    &
                     - 2.0_dp*a_j*mu_j**2 - a_j*beta**2/(2.0_dp*alpha**2)

            dlog(2) =  Im*(-mu_j**2 - 0.5_dp/alpha + mu_j*beta/alpha          &
                           - beta**2/(4.0_dp*alpha**2))

            dlog(3) =  bA_j*(beta/alpha - 2.0_dp*mu_j) - Im*p_j
            dlog(4) =  Im*(beta/(2.0_dp*alpha) - mu_j)

            dSmat(i,j,:) = Sij * dlog
         end do
      end do

      ! add bra-side contribution on the diagonal (matches the Python version)
      do i = 1, n
         dSmat(i,i,:) = dSmat(i,i,:) + conjg(dSmat(i,i,:))
      end do
      if (.not. present(Smat_in))  deallocate(Smat)
   end function dS


   function dXn(p, order, Smat_in, dSmat_in) result(dXmat)
      use, intrinsic :: iso_fortran_env, only : dp => real64
      implicit none

      ! ---- dummy arguments -------------------------------------------------
      complex(dp), allocatable       :: dXmat(:,:,:), dprefac(:,:,:) ! (n,n,4)
      real(dp), intent(in)           :: p(:,:)            ! (n,4)
      integer, intent(in)           :: order            ! order of the derivative
      complex(dp), intent(in), optional, target :: Smat_in(:,:)      ! (n,n)
      complex(dp), intent(in), optional, target :: dSmat_in(:,:,:)   ! (n,n,4)

      ! ---- local working arrays -------------------------------------------
      complex(dp), pointer          :: Smat(:,:), dSmat(:,:,:)

      integer :: n, i, j
      real(dp) :: a_i,b_i,mu_i,p_i,  a_j,b_j,mu_j,p_j
      complex(dp) :: bA_i,bA_j, alpha, beta, alpha2, beta2,alpha3, alpha4, beta4, alpha5, alpha6

      ! ---------------------------------------------------------------------
      n = size(p,1)

      ! --- S matrix ---------------------------------------------------------
      if (present(Smat_in)) then
         Smat => Smat_in             ! alloc/copy handled automatically
      else
         allocate(Smat(n,n))
         Smat = S(p)
      end if

      ! --- dS matrix --------------------------------------------------------
      if (present(dSmat_in)) then
         dSmat => dSmat_in
      else
         allocate(dSmat(n,n,4))
         dSmat = dS(p)
      end if

      ! --- arrays we *always* need -----------------------------------------
      allocate(dXmat(n,n,4), dprefac(n,n,4))
      select case (order)
         case(1)
            do  j= 1, n
               a_j = p(j,1);  b_j = p(j,2);  mu_j = p(j,3);  p_j = p(j,4)
               bA_j = cmplx(a_j**2, b_j, kind=dp)

               do i= 1, n
                  a_i = p(i,1);  b_i = p(i,2);  mu_i = p(i,3);  p_i = p(i,4)
                  bA_i = cmplx(a_i**2, b_i, kind=dp)

                  call moments(bA_i, bA_j, mu_i, mu_j, p_i, p_j, &
                           alpha, alpha2, alpha3,alpha4,alpha5,alpha6, beta, beta2,beta4)
                  dXmat(i,j,:)=beta / (2.0_dp*alpha)*dSmat(i,j,:) ! prefac * dS_ij/dpj
                  ! Calculate prefac derivatives
                  dprefac(i,j,1) = a_j*(2.0_dp*alpha*mu_j - beta)/alpha**2
                  dprefac(i,j,2)  =  Im*(alpha*mu_j - 0.5_dp*beta)/alpha**2
                  dprefac(i,j,3)  =  bA_j/alpha
                  dprefac(i,j,4)  =  0.5_dp*Im/alpha

               end do
            end do
         case(2)
            do  j= 1, n
               a_j = p(j,1);  b_j = p(j,2);  mu_j = p(j,3);  p_j = p(j,4)
               bA_j = cmplx(a_j**2, b_j, kind=dp)

               do i= 1, n
                  a_i = p(i,1);  b_i = p(i,2);  mu_i = p(i,3);  p_i = p(i,4)
                  bA_i = cmplx(a_i**2, b_i, kind=dp)

                  call moments(bA_i, bA_j, mu_i, mu_j, p_i, p_j, &
                           alpha, alpha2, alpha3,alpha4,alpha5,alpha6, beta, beta2,beta4)

                  dprefac(i,j,1)  =  a_j*(2.0_dp*alpha*beta*mu_j - alpha - beta2)/alpha3
                  dprefac(i,j,2)  =  Im*(2.0_dp*alpha*beta*mu_j - alpha - beta2)/(2.0_dp*alpha3)
                  dprefac(i,j,3)  =  beta*bA_j/alpha2
                  dprefac(i,j,4)  =  Im*beta/(2.0_dp*alpha2)
                  dXmat(i,j,:)=(beta2/(4.0_dp*alpha2) + 1.0_dp/(2.0_dp*alpha))*dSmat(i,j,:) ! prefac * dS_ij/dpj

               end do
            end do
         case(3)
            do  j= 1, n
               a_j = p(j,1);  b_j = p(j,2);  mu_j = p(j,3);  p_j = p(j,4)
               bA_j = cmplx(a_j**2, b_j, kind=dp)

               do i= 1, n
                  a_i = p(i,1);  b_i = p(i,2);  mu_i = p(i,3);  p_i = p(i,4)
                  bA_i = cmplx(a_i**2, b_i, kind=dp)

                  call moments(bA_i, bA_j, mu_i, mu_j, p_i, p_j, &
                           alpha, alpha2, alpha3,alpha4,alpha5,alpha6, beta, beta2,beta4)
                  dprefac(i,j,1)  =  3.0_dp*a_j*(2.0_dp*alpha*mu_j*(2.0_dp*alpha + beta2) - beta*(4.0_dp*alpha + beta2))/(4.0_dp*alpha4)
                  dprefac(i,j,2)  =  3.0_dp*Im*(2.0_dp*alpha*mu_j*(2.0_dp*alpha + beta2) - beta*(4.0_dp*alpha + beta2))/(8.0_dp*alpha4)
                  dprefac(i,j,3)  =  3.0_dp*bA_j*(2.0_dp*alpha + beta2)/(4.0_dp*alpha3)
                  dprefac(i,j,4)  =  3.0_dp*Im*(2.0_dp*alpha + beta2)/(8.0_dp*alpha3)
                  dXmat(i,j,:) = (3.0_dp/4.0_dp * beta/(alpha**2) + 1.0_dp/8.0_dp* (beta/alpha)**3)*dSmat(i,j,:) ! prefac * dS_ij/dpj

               end do
            end do
         case(4)
            do  j= 1, n
               a_j = p(j,1);  b_j = p(j,2);  mu_j = p(j,3);  p_j = p(j,4)
               bA_j = cmplx(a_j**2, b_j, kind=dp)
               do i= 1, n
                  a_i = p(i,1);  b_i = p(i,2);  mu_i = p(i,3);  p_i = p(i,4)
                  bA_i = cmplx(a_i**2, b_i, kind=dp)

                  call moments(bA_i, bA_j, mu_i, mu_j, p_i, p_j, &
                           alpha, alpha2, alpha3,alpha4,alpha5,alpha6, beta, beta2,beta4)

                  dprefac(i,j,1)  =  -a_j*(6.0_dp*alpha2 + 9.0_dp*alpha*beta2 - 2.0_dp*alpha*beta*mu_j*(6.0_dp*alpha + beta2) + beta4)/(2.0_dp*alpha5)
                  dprefac(i,j,2)  =  Im*(-6.0_dp*alpha2 - 9.0_dp*alpha*beta2 + 2.0_dp*alpha*beta*mu_j*(6.0_dp*alpha + beta2) - beta4)/(4.0_dp*alpha5)
                  dprefac(i,j,3)  =  beta*bA_j*(6.0_dp*alpha + beta2)/(2.0_dp*alpha4)
                  dprefac(i,j,4)  =  Im*beta*(6.0_dp*alpha + beta2)/(4.0_dp*alpha4)
                  dXmat(i,j,:) = (12.0_dp*alpha2 + 12.0_dp*alpha*beta2 + beta4)/(16.0_dp*alpha4)*dSmat(i,j,:) ! prefac * dS_ij/dpj
               end do
            end do
         case default
            write(*,*) "dXn: order not implemented"
            stop 1
      end select
      
      do i = 1,n
         dprefac(i,i,:)=dprefac(i,i,:) + conjg(dprefac(i,i,:))
      end do
      do j = 1,n
         do i =1,n
            dXmat(i,j,:)=dXmat(i,j,:)+dprefac(i,j,:)*Smat(i,j)
         end do
      end do
      if (.not. present(Smat_in))  deallocate(Smat)
      if (.not. present(dSmat_in)) deallocate(dSmat)
      deallocate(dprefac)   ! always ours   end function dX2
   end function dXn
  function dH(p,t) result(dHmat)
      use, intrinsic :: iso_fortran_env, only : dp => real64
      implicit none
      real(dp),    intent(in)  :: p(:,:)
      complex(dp), intent(in) :: t
      complex(dp), allocatable :: dHmat(:,:,:) ! (n,n,4)
      complex(dp), allocatable :: dSmat(:,:,:), dXmat(:,:,:), dX2mat(:,:,:)
      complex(dp), allocatable :: Smat(:,:), Xmat(:,:), X2mat(:,:)
      integer :: n, i, j

      n = size(p,1)
      allocate(Smat(n,n), dSmat(n,n,4))
      Smat= S(p) ! (n,n)
      dSmat = dS(p,Smat) ! 
      Xmat = xn(p,1)*Smat ! (n,n)
      X2mat = xn(p,2)*Smat ! (n,n)
      allocate(dXmat(n,n,4), dX2mat(n,n,4))
      dXmat = dXn(p,1,Smat,dSmat) ! (n,n,4)
      dX2mat = dXn(p,2,Smat,dSmat) !
   deallocate(Smat,Xmat,X2mat,dSmat,dXmat)
   end function dH
end module ho_kernels
