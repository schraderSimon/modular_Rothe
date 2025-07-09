program vector_matrix_product
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none
    integer, parameter :: size = 3 ! Fixed matrix and vector size
    integer :: i, j
    complex(dp) :: cv(0:size-1), cm(0:size-1, 0:size-1), cresult(0:size-1)
    complex(dp) :: z1, z2, z3
    real(dp) :: rpart(size), ipart(size)
    real(dp) :: rmat(0:size-1, 0:size-1), imat(0:size-1, 0:size-1)
    call random_number(rpart) ! Generate random numbers for real part
    call random_number(ipart) ! Generate random numbers for imaginary part
    call random_number(rmat) ! Generate random numbers for real part of matrix
    call random_number(imat) ! Generate random numbers for imaginary part of matrix
    
    ! Alternatively, just assign values to the matrix
    do concurrent(i = 0:size-1)
        do concurrent(j=0:size-1)
            rmat(i,j) = real(i + j, kind=dp) ! Example values for real part
            imat(i,j) = real(i - j, kind=dp) ! Example values for imaginary part
        end do
    end do
    cv = cmplx(rpart, ipart, kind=dp)
    cm = cmplx(rmat, imat, kind=dp)
    cresult = matmul(cm, cv)
    z1= (1.0_dp,2.0_dp)
    z2 = (3.0_dp, 4.0_dp)
    z3 = z1 + z2
    print *, 'z1:', z1
    print *, 'z2:', z2
    print *, 'z3:', z3
    print *, 'Matrix cm:'
    do i = 0, size-1
        write(*,'(3("(",F8.3,"+",F8.3,"j) "))') cm(i,:)
    end do
    print *, 'Vector cv:'
    do i = 0, size-1
        write(*, '(F10.4)') cv(i)
    end do
    print *, 'Result of matrix-vector multiplication cresult:'
    do i = 0, size-1
        write(*, '(F10.4)') cresult(i)
    end do

end program vector_matrix_product