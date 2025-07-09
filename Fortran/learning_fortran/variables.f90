program variables
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none
    real(dp) :: x
    character :: name
    integer :: i
    ! Initialize variables
    i = 10
    x = 3.14_dp
    name = 'A'
    ! Print variables
    print *, 'Integer i:', i
    print *, 'Real x:', x
    print *, 'Character name:', name
    ! Sum of i and x
    print *, 'Sum of i and x:', i + x
end program variables