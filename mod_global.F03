!-------------------------------------------------------------------------------
! module mod_global :
!     Global variables
!
module mod_global
    use, intrinsic            :: iso_c_binding
    implicit none
    integer(c_int), parameter :: DIM=3
    integer(c_int), parameter :: cd = c_double ! used for literal constants
end module mod_global
