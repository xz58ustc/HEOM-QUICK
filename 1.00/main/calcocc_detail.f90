subroutine calcocc_detail
use matmod
use tmpmatmod
use sparsemod
implicit none
include '../include/sizes'
include '../include/common'
!
integer    :: ni, nj, nk, nnz
integer*8  :: lni
complex*16 :: ctmp1
!
! population_ms = tr( c^\dag_{ms} c_{ms} \rho )
!
pocc(1:norbs,1:nspin) = 0.d0
if (.not. lsparse) then
   do nk=1,nspin
      do ni=1,norbs
         call zamsmm1('l', 'n', 'n', ni, nk, cunity, rho(1,1,1), nrho, czero, cmtmp1, nrho)
         call zams_trace_ab('t', ni, nk, cmtmp1, nrho, ctmp1)
         pocc(ni,nk) = dble(ctmp1)
      end do
   end do
else
   nnz = nnz_spa(1)
   lni = ind_spa(1)
   do nk=1,nspin
      do ni=1,norbs
         call zamsmm_coo('l', 'n', 'n', ni, nk, cunity, rho_spa(lni), nnz, irow_spa(lni), icol_spa(lni), &
                         czero, cmtmp1, nrho, cmtmp2, nrho)
         call zams_trace_ab('t', ni, nk, cmtmp1, nrho, ctmp1)
         pocc(ni,nk) = dble(ctmp1)
      end do
   end do
end if
!
return
end subroutine calcocc_detail
