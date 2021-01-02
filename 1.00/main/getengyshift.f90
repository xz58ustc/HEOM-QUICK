subroutine getengyshift(tt, ialf, ispin, eshift)
use matmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent(in)  :: ialf, ispin     ! 1 for left, 2 for right
real*8,  intent(in)  :: tt
real*8,  intent(out) :: eshift
!
if (tt >= t_off) then
   eshift = 0.d0
   return 
end if
!
if (fieldtype .eq. 0) then
   eshift = amps(ialf,ispin) * (1.d0 - dexp(-tt / tchar(ialf,ispin)) )
!
else if (fieldtype .eq. 1) then
   eshift = amps(ialf,ispin) * dsin(2.0 * pi * tt / tchar(ialf,ispin))
!
else if (fieldtype .eq. -1) then
   eshift = 0.d0
else
   write(6,*)'getengyshift: error! unknown fieldtype =', fieldtype
   stop
end if
!
end subroutine getengyshift
