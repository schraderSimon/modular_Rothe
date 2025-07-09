module lapack_interfaces
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none
    interface
        subroutine zgetrf(m, n, a, lda, ipiv, info) !Compute LU factorization of a complex matrix
            integer, intent(in) :: m, n, lda
            integer, intent(out) :: ipiv(*)
            integer, intent(out) :: info
            complex*16, intent(inout) :: a(lda, *)
        end subroutine zgetrf

        subroutine zgetri(n, a, lda, ipiv, work, lwork, info) !Compute the inverse of a complex matrix using LU factorization
            integer, intent(in) :: n, lda, lwork
            integer, intent(in) :: ipiv(*)
            integer, intent(out) :: info
            complex*16, intent(inout) :: a(lda, *)
            complex*16, intent(inout) :: work(lwork)
        end subroutine zgetri

        subroutine zhesv(uplo, n, nrhs, a, lda, ipiv, b, ldb, info) !Solve a system of linear equations with a complex Hermitian matrix
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            integer, intent(out) :: ipiv(*)
            integer, intent(out) :: info
            complex*16, intent(inout) :: a(lda, *)
            complex*16, intent(inout) :: b(ldb, *)
        end subroutine zhesv
        subroutine zposv(uplo, n, nrhs, a, lda, b, ldb, info)
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            integer, intent(out) :: info
            complex*16, intent(inout) :: a(lda, *)
            complex*16, intent(inout) :: b(ldb, *)
        end subroutine zposv    
    end interface
end module lapack_interfaces
module matrix_operations
    use lapack_interfaces
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
contains
    subroutine invert(A,Ainv,check)
        implicit none
        ! Invert the matrix A and store the result in B
        complex(dp), intent(in) :: A(:,:)
        complex(dp), intent(out) :: Ainv(:,:)
        complex(dp), allocatable :: A_times_Ainv(:,:)
        logical, intent(in) :: check ! Wether or not to check the result
        
        integer :: info,n,i,j
        real(dp) :: error ! Size of the matrix for checking
        complex(dp) :: element
        
        complex(dp), allocatable :: work(:) ! Work array for zgetri
        integer, allocatable :: ipiv(:) ! ipiv is now allocatable
        
        n = size(A, 1) ! Assuming A is square
        allocate(ipiv(n)) ! Allocate pivot array of size n
        allocate(work(n))
        allocate(A_times_Ainv(n, n)) ! Allocate the result matrix
        ! Perform LU decomposition
        Ainv = A ! Initialize Ainv with A
        call zgetrf(n, n, Ainv, n, ipiv, info) ! Perform LU factorization, storing the result in Ainv
        call zgetri(n, Ainv, n, ipiv, work, n, info)
        deallocate(ipiv) ! Deallocate pivot array
        deallocate(work) ! Deallocate work array
        A_times_Ainv= matmul(A, Ainv) ! Compute the product of A and Ainv
        if (check) then
            do i=1, n
                do j=1, n
                    if (i==j) then
                        element= A_times_Ainv(i,j) - (1.0_dp, 0.0_dp) ! Check diagonal elements
                        error= abs(element)
                    else
                        element= A_times_Ainv(i,j)
                        error= abs(element)
                    end if
                    !print * , A_times_Ainv(i,j), 'Element (', i, ',', j, ')'
                    if (error > 1.0e-6_dp) then ! Check if the error is within a tolerance
                        !print *, 'Error in element (', i, ',', j, '): ', error
                    end if
                end do
            end do
        end if
        
        deallocate(A_times_Ainv) ! Deallocate the result matrix
        print *, 'Matrix inversion completed successfully.'
    end subroutine invert
    function inv(A) result(Ainv)
        implicit none
        ! Invert the matrix A and store the result in B
        complex(dp), intent(in) :: A(:,:)
        complex(dp), allocatable :: Ainv(:,:)
        logical :: check

        check = .false. ! Set to true to check the result
        allocate(Ainv(size(A, 1), size(A, 2))) ! Allocate space for the inverse matrix
        call invert(A, Ainv, check) ! Call the main inversion subroutine
    end function inv
    function solve_posdef(A,b) result(c)
        implicit none
        integer ::  n, nrhs, info
        complex(dp), intent(in):: A(:,:)
        complex(dp), intent(in):: b(:)
        complex(dp), allocatable :: c(:)
        complex(dp), allocatable :: Acopy(:,:), bcopy(:,:)
        n= size(A,1)
        nrhs = 1
        allocate(c(n))
        allocate(bcopy(n,1))
        bcopy(:,1)=b
        Acopy=A
        call zposv('U', n, nrhs, Acopy, n, bcopy, n, info)
        if (info /= 0) then
            print *, 'zposv failed with info =', info
        stop 1
    end if
        c=bcopy(:,1)
    end function solve_posdef
end module matrix_operations
program test_invert
    use matrix_operations
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none
    integer :: i,j, n
    complex(dp), allocatable :: A(:,:), Ainv(:,:), Ainv2(:,:), vector(:), sol1(:), sol2(:), ip1, ip2
    real(dp), allocatable :: rpart(:), ipart(:) ! Real and imaginary parts for random numbers
    real(dp) :: t_start, t_end, t_start_s,t_end_s
    logical :: check
    character(len=32) :: arg
    integer :: istat

    call get_command_argument(1, arg)
    read(arg, *, iostat=istat) n
    if (istat /= 0 .or. n < 1) then
        print *, 'Usage: ./program n'
        stop 1
    end if
    allocate(A(n, n))
    do i = 1, n
        do j = 1, n
            if (i /= j) then
                A(i,j) = cmplx(1.0_dp/(i+j-1), 0.0_dp, kind=dp)  
            else 
                A(i,j)= cmplx(1.0_dp, 0.0_dp, kind=dp) ! Diagonal elements
            end if
            A(i,j) = A(i,j) + cmplx(0.0_dp, 1.0_dp, kind=dp) * real(i+j, kind=dp) * 0.1_dp ! Adding small imaginary part
        end do
    end do
    A=matmul(A,conjg(transpose(A)))
    ! Allocate space for the inverse matrix
    allocate(Ainv(n, n))
    allocate(Ainv2(n, n)) ! Allocate space for the inverse matrix from the inv subroutine
    allocate(rpart(n), ipart(n), vector(n), sol1(n), sol2(n)) ! Allocate space for random parts and solutions
    check = .true. ! Set to true to check the result
    !call cpu_time(t_start)
    !call invert(A, Ainv, check)
    !call cpu_time(t_end)
    !print *, 'Matrix A:'
    !print *, A
    !print *, 'Inverse Matrix Ainv:'
    !print *, Ainv
    !print*, 'Check result:'
    !print*, matmul(A, Ainv) ! Print the product of A and Ainv to verify the inversion
    !Ainv2= inv(A) ! Call the inv subroutine to test the interface
    !print *, 'Inverse Matrix Ainv from inv subroutine:'
    !print *, Ainv2 ! Print the inverse matrix from the inv subroutine
    !print *, 'Check result from inv subroutine:'
    !print *, matmul(A, Ainv2) ! Print the product of A and Ainv
    
    call random_number(rpart)
    ipart=0
    vector=cmplx(rpart,ipart,dp)
    ! Solve linear equation using invert
    sol1=matmul(Ainv,vector)
    ! from lapack import solve_posdef
    call cpu_time(t_start_s)
    sol2=solve_posdef(A,vector)
    call cpu_time(t_end_s)
    !print *, "Solution using inverse"
    !print *,sol1
    !print *, "Solution using solve"
    !print *, sol2
    ip1=sum(conjg(sol1)*sol1)
    ip2=sum(conjg(sol2)*sol2)
    print *, "Relative norm difference"
    print *, (ip1-ip2)/ip1
    print *, "Time taken for inversion"
    print *, t_end-t_start
    print *, "Time taken for solve"
    print *, t_end_s-t_start_s
    deallocate(A, Ainv, Ainv2, rpart, ipart, vector, sol1, sol2) ! Deallocate all allocated arrays
end program test_invert