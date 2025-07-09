
program ho_test
   use, intrinsic :: iso_fortran_env, only : dp => real64
   use ho_kernels
   implicit none
   real(dp), parameter :: a = 1.0_dp / sqrt(2.0_dp)   ! ground-state width, ω = 1
   real(dp), allocatable :: p(:,:)
   complex(dp) :: Hmat(1,1), Smat(1,1)
   real(dp)    :: E

   allocate(p(1,4))
   p = 0.0_dp                ! (a, b, μ, p)
   p(1,1) = a                ! a = 1/√2

   Hmat = H(p)
   Smat = S(p)
   print '(A, F10.7)', 'normalization 〈S〉 = ',real(Smat(1,1))
   E    = real( Hmat(1,1) / Smat(1,1) )   ! expectation value

   print '(A, F10.7)', 'Energy 〈H〉 = ', E
   if (abs(E - 0.5_dp) < 1.0e-10_dp) then
      print *, '  Test passed round-state energy is 0.5'
   else
      print *, '  Test failed something is off'
   end if
end program ho_test