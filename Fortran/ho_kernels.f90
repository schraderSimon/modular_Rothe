module ho_kernels
   use, intrinsic :: iso_fortran_env, only : dp => real64
   implicit none
   private
   public :: S, x, x2, H, H2, x3, x4
   real(dp), parameter         :: pi = 3.141592653589793238462643383279502884197169399375_dp
   complex(dp), parameter      :: Im = (0.0_dp, 1.0_dp)
!! Parameters you don’t want to thread through every call
   real(dp), parameter :: omega = 1.0_dp

contains
!────────────── overlap S ──────────────
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
        normalizations(i) = (2.0_dp * a_i**2.0_dp / pi)**0.25_dp
    end do

    do i = 1, n
        do j = i, n
            N_i= cmplx(normalizations(i),0.0_dp, kind=dp) ! Normalization factor for Gaussian i
            N_j= cmplx(normalizations(j),0.0_dp, kind=dp) ! Normalization factor for Gaussian j
            a_i = p(i,1) ! a_i
            b_i = p(i,2) ! b_i
            mu_i = p(i,3) ! mu_i
            p_i = p(i,4) ! p_i
            a_j = p(j,1) ! a_j
            b_j = p(j,2) ! b_j
            mu_j = p(j,3) ! mu_j
            p_j = p(j,4) ! p_j
            bA_i = cmplx(a_i**2.0_dp, b_i, kind=dp)
            bA_j = cmplx(a_j**2.0_dp, b_j, kind=dp)
            alpha = conjg(bA_i) + bA_j
            beta = 2.0_dp*(conjg(bA_i)*mu_i + bA_j*mu_j) + Im*(p_j - p_i)
            gamma = -(conjg(bA_i)*mu_i**2.0_dp + bA_j*mu_j**2.0_dp) + Im*(p_i*mu_i - p_j*mu_j)
            Z = N_i * N_j * sqrt(pi/alpha) * exp(gamma + beta**2.0_dp/(4.0_dp*alpha))
            Smat(i,j) = Z
            Smat(j,i) = conjg(Z) ! Ensure symmetry
        end do
    end do
    deallocate(normalizations)
      !– ensure symmetry –
   Smat = 0.5_dp * (Smat + conjg(transpose(Smat)))
    end function S
    function x(p) result(Xmat)
      real(dp),    intent(in)  :: p(:,:)       ! (n,4)
      complex(dp), allocatable :: Xmat(:,:)   ! (n,n)
      complex(dp), allocatable :: Smat(:,:)
      integer :: n, i, j
      complex(dp) :: bA_i,bA_j, alpha,beta, fact
      real(dp) :: a_i,b_i,mu_i,p_i, a_j,b_j,mu_j,p_j

      n = size(p,1)
      allocate(Xmat(n,n))

      do i=1,n
         do j=i,n
            a_i = p(i,1); b_i = p(i,2); mu_i = p(i,3); p_i = p(i,4)
            a_j = p(j,1); b_j = p(j,2); mu_j = p(j,3); p_j = p(j,4)

            bA_i = cmplx(a_i**2,b_i,kind=dp);      bA_j = cmplx(a_j**2,b_j,kind=dp)

            alpha = conjg(bA_i) + bA_j
            beta  = 2.0_dp*(conjg(bA_i)*mu_i + bA_j*mu_j) + Im*(p_j-p_i)

            fact = beta / (2.0_dp*alpha)
            Xmat(i,j) = fact
            Xmat(j,i) = conjg(Xmat(i,j))
         end do
      end do

   end function x


!────────────────── matrix ⟨x²⟩ ──────────────────
    function x2(p) result(X2mat)
      real(dp),    intent(in)  :: p(:,:)       ! (n,4)
      complex(dp), allocatable :: X2mat(:,:)   ! (n,n)
      complex(dp), allocatable :: Smat(:,:)
      integer :: n, i, j
      complex(dp) :: bA_i,bA_j, alpha,beta, fact
      real(dp) :: a_i,b_i,mu_i,p_i, a_j,b_j,mu_j,p_j

      n = size(p,1)
      allocate(X2mat(n,n))
      do i=1,n
         do j=i,n
            a_i = p(i,1); b_i = p(i,2); mu_i = p(i,3); p_i = p(i,4)
            a_j = p(j,1); b_j = p(j,2); mu_j = p(j,3); p_j = p(j,4)

            bA_i = cmplx(a_i**2,b_i,kind=dp);      bA_j = cmplx(a_j**2,b_j,kind=dp)

            alpha = conjg(bA_i) + bA_j
            beta  = 2.0_dp*(conjg(bA_i)*mu_i + bA_j*mu_j) + Im*(p_j-p_i)
            !gamma = -(conjg(bA_i)*mu_i**2 + bA_j*mu_j**2) + Im*(p_i*mu_i - p_j*mu_j)

            fact = beta**2/(4.0_dp*alpha**2) + 1.0_dp/(2.0_dp*alpha)
            X2mat(i,j) =  fact
            X2mat(j,i) = conjg(X2mat(i,j))
         end do
      end do

   end function x2
function x3(p) result(X3mat)
      real(dp),    intent(in)  :: p(:,:)       ! (n,4)
      complex(dp), allocatable :: X3mat(:,:)   ! (n,n)
      complex(dp), allocatable :: Smat(:,:)
      integer :: n, i, j
      complex(dp) :: bA_i,bA_j, alpha,beta, fact
      real(dp) :: a_i,b_i,mu_i,p_i, a_j,b_j,mu_j,p_j

      n = size(p,1)
      allocate(X3mat(n,n))

      do i=1,n
         do j=i,n
            a_i = p(i,1); b_i = p(i,2); mu_i = p(i,3); p_i = p(i,4)
            a_j = p(j,1); b_j = p(j,2); mu_j = p(j,3); p_j = p(j,4)

            bA_i = cmplx(a_i**2,b_i,kind=dp);      bA_j = cmplx(a_j**2,b_j,kind=dp)

            alpha = conjg(bA_i) + bA_j
            beta  = 2.0_dp*(conjg(bA_i)*mu_i + bA_j*mu_j) + Im*(p_j-p_i)
            !gamma = -(conjg(bA_i)*mu_i**2 + bA_j*mu_j**2) + Im*(p_i*mu_i - p_j*mu_j)

            fact = 3.0_dp/4.0_dp * beta/(alpha**2) + 1.0_dp/8.0_dp* (beta/alpha)**3
            X3mat(i,j) =  fact
            X3mat(j,i) = conjg(X3mat(i,j))
         end do
      end do

   end function x3
function x4(p) result(X4mat)
      real(dp),    intent(in)  :: p(:,:)       ! (n,4)
      complex(dp), allocatable :: X4mat(:,:)   ! (n,n)
      complex(dp), allocatable :: Smat(:,:)
      integer :: n, i, j
      complex(dp) :: bA_i,bA_j, alpha,beta, fact
      real(dp) :: a_i,b_i,mu_i,p_i, a_j,b_j,mu_j,p_j

      n = size(p,1)
      allocate(X4mat(n,n))
      do i=1,n
         do j=i,n
            a_i = p(i,1); b_i = p(i,2); mu_i = p(i,3); p_i = p(i,4)
            a_j = p(j,1); b_j = p(j,2); mu_j = p(j,3); p_j = p(j,4)

            bA_i = cmplx(a_i**2,b_i,kind=dp);      bA_j = cmplx(a_j**2,b_j,kind=dp)

            alpha = conjg(bA_i) + bA_j
            beta  = 2.0_dp*(conjg(bA_i)*mu_i + bA_j*mu_j) + Im*(p_j-p_i)
            !gamma = -(conjg(bA_i)*mu_i**2 + bA_j*mu_j**2) + Im*(p_i*mu_i - p_j*mu_j)

            fact = 1/(2.0_dp*alpha**2)+beta**2/(2.0_dp*alpha**3) + (1/(2.0_dp*alpha)+(beta/(2.0_dp*alpha))**2)**2
            X4mat(i,j) = fact
            X4mat(j,i) = conjg(X4mat(i,j))
         end do
      end do

   end function x4
   !────────────────── Hamiltonian H = –½∇² + ½ ω² x² ──────────────────
   function H(p,t) result(Hmat)
      use, intrinsic :: iso_fortran_env, only : dp => real64
      implicit none
      real(dp),    intent(in)  :: p(:,:), t
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
      Xmat = x (p)*Smat      ! ⟨g_i|x|g_j⟩
      X2mat = x2(p)*Smat     ! ⟨g_i|x²|g_j⟩

      do j = 1, n                      ! column → parameters of |g_j⟩
         a_j  = p(j,1);  b_j = p(j,2)
         mu_j = p(j,3);  p_j = p(j,4)
         bA_j  = cmplx(a_j**2, b_j, kind=dp)

         do i = 1, n                   ! row → ⟨g_i|
            !──────── kinetic part  –½⟨∂²⟩  (see derivation in analysis) ────────
            Hmat(i,j) =  -2.0_dp*bA_j**2 * ( X2mat(i,j) - 2.0_dp*mu_j*Xmat(i,j) + mu_j**2.0_dp*Smat(i,j) )   &
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
      real(dp),    intent(in)  :: p(:,:),t
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
      Xmat = x (p)!*Smat      ! ⟨g_i|x|g_j⟩
      X2mat = x2(p)!*Smat     ! ⟨g_i|x²|g_j⟩
      X3mat = x3(p)!*Smat
      X4mat = x4(p)!*Smat
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
            c(5) = 4_dp*bA_i_c**2.0_dp*mu_i**2 - 4_dp*Im*bA_i_c*mu_i*p_i - 2_dp*bA_i_c - p_i**2
            c(6) = 4_dp*bA_i_c*(-2_dp*bA_i_c*mu_i + Im*p_i)
            c(7) = 4_dp*bA_i_c**2
            H2mat(i,j) = H2mat(i,j) - 0.25_dp * (c(5)*X2mat(i,j) + c(6)*X3mat(i,j) + c(7)*X4mat(i,j)) * omega**2 
            ! Potential times kinetic
            c(8) = 4_dp*bA_j**2.0_dp*mu_j**2 + 4_dp*Im*bA_j*mu_j*p_j - 2_dp*bA_j - p_j**2
            c(9) = 4_dp*bA_j*(-2.0_dp*bA_j*mu_j - Im*p_j)
            c(10) = 4_dp*bA_j**2
            H2mat(i,j) = H2mat(i,j) - 0.25_dp * (c(8)*X2mat(i,j) + c(9)*X3mat(i,j) + c(10)*X4mat(i,j)) * omega**2 
         end do
      end do
      H2mat=H2mat*Smat ! Don't forget that all matrix elements need to be multiplied with the overlap!
      deallocate(Smat, Xmat, X2mat, X3mat, X4mat)
   end function H2
end module ho_kernels