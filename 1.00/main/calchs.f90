subroutine calchs(rhoinp)
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
complex*16, intent(in) :: rhoinp(nrho,*)
!
integer :: ni, nj, nk, iorbs1, ispin, ispin2
real*8  :: dtmp1, dtmp2, dtmp3, dtmp4, dtmp5
real*8  :: delta0, delta1, tksymm, tkasym, edotp
logical :: lherm
!
namelist / para1 / eup, edown, uu, tupdn, fixdot
namelist / para2 / edot, fixdot
namelist / para3 / engy01, engy02, u12, t12, phi12, fixdot
namelist / para4 / e1up, e1down, e2up, e2down, u12, t12, uu1, uu2, j12, fixdot
namelist / sys3level / engy01, engy02, engy03, t12, t23, t13, fixdot
namelist / para_hubbard / lhubbard, edot, vcoup, uu, lhf, lspinsymm, lorbsymm, lspindiff, eup, edown, fixdot
namelist / hamil_sys / lreadhs, lfixhs, lgenhs
!
integer :: itimes
data itimes /0/
!
! single-site, spin-splitted
! 
! H_S = eup * c^\dag_up * c_up + edown * c^\dag_down * c_down +
!        uu * c^\dag_up * c_up * c^\dag_down * c_down
!
itimes = itimes + 1
!
if (itimes > 1) then
   if (lfixhs) return
end if
!
if (itimes .eq. 1) then
   lhubbard = .false.
   lhf      = .false. 
   lspinsymm = .false.
   lorbsymm  = .false.
   lspindiff = .false.
   vcoup    = 0.d0
   edot     = 0.d0
   eup      = 0.d0 
   edown    = 0.d0
   uu       = 0.d0
   fixdot   = .true.
   rewind(5)
   read(5, para_hubbard, end=51) 
   51 continue
!
   if (lhubbard .and. lspindiff) then
      if (nspin .ne. 2 .or. lspinsymm .eq. .true.) then   
         write(6,*)'calchs: error! intrinsic input error for lspindiff=T case '
         stop
      end if
   end if   
!
   lreadhs = .false.
   lfixhs  = .false.
   lgenhs  = .false.
   rewind(5)
   read(5, hamil_sys, end=52) 
   52 continue
   if (lfixhs) then
      fixdot = .true.
   end if
   if (lreadhs) then
      open(unit=40, file='hamil_sys.data', form='unformatted', status='unknown')
      rewind(40)
      read(40) ((hs(ni,nj), ni=1,nrho), nj=1,nrho)
      close(40)
      write(6,*)'calchs: <hamil_sys.data> read successfully'
      goto 999
   end if
end if
!
hs = czero
!
if (lgenhs) then
   call genhs
   goto 999
end if
!
if (lhubbard) then
   if (itimes .eq. 1) then
      write(6,*)'calchs: Hubbard model Hamiltonian is used, parameters are '
      write(6,*)'calchs: lhf    = ', lhf
      write(6,*)'calchs: vcoup  = ', vcoup
      write(6,*)'calchs: uu     = ', uu
      write(6,*)'calchs: fixdot = ', fixdot
      write(6,*)'calchs: lspinsymm = ', lspinsymm
      write(6,*)'calchs: lorbsymm  = ', lorbsymm
      write(6,*)'calchs: lspindiff = ', lspindiff
      if (lspindiff) then
         write(6,*)'calchs: eup    = ', eup 
         write(6,*)'calchs: edown  = ', edown
      else
         write(6,*)'calchs: edot   = ', edot
      end if
      call flush(6)
      edot  = edot / hbar
      eup   = eup  / hbar
      edown = edown / hbar
      vcoup = vcoup / hbar
      uu    = uu / hbar
   end if
!
   if (lhf) then  
!
     call calc_rsdm(rhoinp(1,1), nrho, norbs, nspin, rsdm)
!
! symmetrize spin
!     
     if (lspinsymm .and. nspin > 1) then
        rsdm(1:norbs,1:norbs,1) = (rsdm(1:norbs,1:norbs,1) + rsdm(1:norbs,1:norbs,2)) * 5.d-1
        rsdm(1:norbs,1:norbs,2) =  rsdm(1:norbs,1:norbs,1)
     end if
!
     if (lorbsymm) then
        if (norbs .eq. 2) then
           do ispin=1,nspin
              rsdm(1,1,ispin) = (rsdm(1,1,ispin) + rsdm(2,2,ispin)) * 5.d-1
              rsdm(2,2,ispin) = rsdm(1,1,ispin)
           end do
        end if
     end if
!
! diagonal term
!
     cmtmp2(1:nrho, 1:nrho) = czero
     do ispin=1,nspin
        do iorbs1=1,norbs 
           dtmp1 = edot 
           if (nspin > 1) then
              !dtmp1 = dtmp1 + uu * 5.d-1 * rsdm(iorbs1,iorbs1,nspin-ispin+1)
              dtmp1 = dtmp1 + uu * rsdm(iorbs1,iorbs1,nspin-ispin+1)
           end if
           cmtmp1(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, iorbs1, ispin), 0.d0)
           call zamsmm1('l', 'c', 'n', iorbs1, ispin, dcmplx(dtmp1,0.d0), cmtmp1, nrho, cunity, cmtmp2, nrho)
        end do
     end do
     hs(1:nrho, 1:nrho) = hs(1:nrho, 1:nrho) + cmtmp2(1:nrho, 1:nrho)
!
! off-diagonal term
!
     cmtmp2(1:nrho, 1:nrho) = czero
     do ispin=1,nspin
        do iorbs1=1,norbs-1
           cmtmp1(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, iorbs1+1, ispin), 0.d0)
           call zamsmm('l', 'c', 'n', iorbs1, ispin, dcmplx(vcoup,0.d0), cmtmp1, nrho, cunity, cmtmp2, nrho)
        end do
     end do
     do nj=1,nrho
        do ni=1,nrho
           hs(ni,nj) = hs(ni,nj) + cmtmp2(ni,nj) + dconjg(cmtmp2(nj,ni))
        end do
     end do      
!
   else
   !   
   !  edot * a_{ms}^dag a_ms
   ! 
      cmtmp2(1:nrho, 1:nrho) = czero
      if (lspindiff) then
         do iorbs1=1,norbs
            cmtmp1(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, iorbs1, 1), 0.d0)
            call zgemm('c', 'n', nrho, nrho, nrho, dcmplx(eup, 0.d0), cmtmp1, nrho, cmtmp1, nrho, &
                       cunity, cmtmp2, nrho)         
            cmtmp1(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, iorbs1, 2), 0.d0)
            call zgemm('c', 'n', nrho, nrho, nrho, dcmplx(edown, 0.d0), cmtmp1, nrho, cmtmp1, nrho, &
                       cunity, cmtmp2, nrho)         
         end do
      else 
         do iorbs1=1,norbs
            do ispin=1,nspin
               cmtmp1(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, iorbs1, ispin), 0.d0)
               call zgemm('c', 'n', nrho, nrho, nrho, dcmplx(edot, 0.d0), cmtmp1, nrho, cmtmp1, nrho, &
                          cunity, cmtmp2, nrho)         
            end do
         end do
      end if
      hs(1:nrho, 1:nrho) = hs(1:nrho, 1:nrho) + cmtmp2(1:nrho, 1:nrho)
   !
   ! vcoup * a_{ms}^dag a_{m+1 s}^dag + H.c.
   !      
      cmtmp3(1:nrho, 1:nrho) = czero
      do iorbs1=1,norbs-1 
         do ispin=1,nspin
            cmtmp1(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, iorbs1, ispin), 0.d0)
            cmtmp2(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, iorbs1+1, ispin), 0.d0)
            call zgemm('c', 'n', nrho, nrho, nrho, dcmplx(vcoup, 0.d0), cmtmp1, nrho, cmtmp2, nrho, &
                       cunity, cmtmp3, nrho)     
         end do
      end do   
      do nj=1,nrho
         do ni=1,nrho
            hs(ni,nj) = hs(ni,nj) + cmtmp3(ni,nj) + dconjg(cmtmp3(nj,ni))
         end do
      end do      
   !
   ! uu / 2 * a_{ms}^dag a_{ms'}^dag a_{ms') a_{ms}
   !
      cmtmp5(1:nrho, 1:nrho) = czero
      do iorbs1=1,norbs
         do ispin=1,nspin  
            cmtmp1(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, iorbs1, ispin), 0.d0)      
            do ispin2=1,nspin
               cmtmp2(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, iorbs1, ispin2), 0.d0)            
               call zgemm('c', 'c', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp2, nrho, &
                          czero, cmtmp3, nrho)  ! a_{ms}^dag a_{ms'}^dag
               call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp3, nrho, cmtmp2, nrho, &
                          czero, cmtmp4, nrho)  ! a_{ms}^dag a_{ms'}^dag a_{ms'}
               call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp4, nrho, cmtmp1, nrho, &
                          czero, cmtmp3, nrho)  ! a_{ms}^dag a_{ms'}^dag a_{ms'} a_{ms}
               cmtmp5(1:nrho, 1:nrho) = cmtmp5(1:nrho, 1:nrho) + cmtmp3(1:nrho, 1:nrho)
            end do
         end do
      end do
      hs(1:nrho, 1:nrho) = hs(1:nrho, 1:nrho) + dcmplx(5.d-1 * uu, 0.d0) * cmtmp5(1:nrho, 1:nrho)
   end if
   goto 999
end if
!
!--------------------- for (norbs=1, nspin=1) ----------------
if (nspin .eq. 1 .and. norbs .eq. 1) then
  edot = 0.d0
  fixdot = .true.
  rewind(5)
  read(5, para2, end=91)
  91 continue
  epara1 = edot
  if (lprths) then
     write(6,*)
     write(6,*)' parameters for Hamiltonian '
     write(6,*)' edot = ', edot
  end if
  edot = edot / hbar
!   
  cmtmp1(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, 1, 1), 0.d0)
  call zgemm('c', 'n', nrho, nrho, nrho, dcmplx(edot, 0.d0), cmtmp1, nrho, cmtmp1, nrho, &
              czero, cmtmp2, nrho)
  hs = hs + cmtmp2
end if
!---------------------end of area for (norbs=1, nspin=1) ----------------
!
if (nspin .eq. 2 .and. norbs .eq. 1) then
  eup   = 0.d0
  edown = 0.d0
  uu    = 0.d0
  tupdn = 0.d0
  fixdot = .true.
  rewind(5)
  read(5, para1, end=101)
  101 continue
  epara1 = eup
  epara2 = edown
  epara3 = uu
!
  delta0 = max(dlwidth(1,1), dlwidth(1,2)) * 5.d-1   ! choose the larger coupling component to estimate T_K
  dtmp3  = min(eup, edown)
  dtmp4  = uu / delta0
!
! symmetric Anderson model Kondo scale
!
  tksymm = dsqrt(uu * delta0 * 5.d-1) * dexp(                   &
           -pi * dtmp4 / 8.d0 + pi / 2.d0 / dtmp4 )
!
! Bethe Ansatz (exact) Kondo energy scale
!
  tkondo  = dsqrt(uu * delta0 * 5.d-1) * dexp(                  &
            pi * dtmp3 * (dtmp3 + uu) / (uu * delta0 * 2.d0 ))
!
  if (lprths) then
     write(6,*)
     write(6,*)' parameters for Hamiltonian '
     write(6,111)eup, edown, uu
     if (dabs(tupdn) .ge. dnano) then
       write(6,120)tupdn
     end if
  end if
  111 format('   eup = ', f14.8, 2x, ' edown = ', f14.8, 2x, ' U = ', f14.8, 2x)
  120 format(' tupdn = ', f14.8, 2x, ' spin-flip mechanism is invoked ')
!  
  if (lprths) then
     write(6,*)
     write(6,*)' Bethe Ansatz Kondo temperature, T_K (energy unit) = ', tkondo
  end if
!
  dtmp1 = pi * (-dtmp3) * (dtmp3 + uu) / (uu * delta0 * 2.d0)
  dtmp2 = dsqrt(uu * delta0 * 5.d-1) * dexp(-dtmp1 + pi**2 / 16.d0 / dtmp1) 
  if (lprths) then
     write(6,*)'   With the additional exponent, T_K (energy unit) = ', dtmp2
     if (dabs(uu + 2.d0 * eup) .le. dpico .and. dabs(uu + 2.d0 * edown) .le. dpico) then
        write(6,*)
        write(6,*)' Symmetric Anderson Kondo temperature              = ', tksymm
     end if
     call flush(6)
  end if
!
! Kondo temperature for asymmetric Anderson model
!
  delta1 = min(eup, edown) + uu * 5.d-1
!
  dtmp2 = delta0  * 2.d0 / pi  * (1.d0 / dabs(delta1 - uu * 5.d-1) + 1.d0 / dabs(delta1 + uu * 5.d-1)) 
  tkondo = 0.182d0 * uu * dsqrt(dtmp2) * dexp(-1.d0 / dtmp2) 
  if (lprths) then
     write(6,*)
     write(6,*)' Haldane''s expression for Kondo temperature        = ', tkondo
  end if
  tkondo = tkondo * dexp(pi**2 / 16.d0 * dtmp2)
  if (lprths) then
     write(6,*)'   With the additional exponent, T_K (energy unit) = ', tkondo
     call flush(6)
  end if
!
!  if (dabs(delta1) .ge. dnano) then
     if (lprths) then
        write(6,*)'calchs: deviation from p-h symmetry = ', delta1 / uu
     end if
     dtmp1 = min(eup, edown)
     dtmp2 = 1.d-4
     nk    = 1000 
     dtmp3 = dtmp1  ! renormalized dot level
     do ni=1,nk
        dtmp4 = dtmp3
        dtmp3 = dtmp1 - delta0 / pi * dlog(-uu / dtmp3)
        if (dabs(dtmp3 - dtmp4) .le. dtmp2) then 
           if (lprths) write(6,*)'calchs: renormalized dot energy = ', dtmp3
           exit
        end if
     end do
     if (dabs(dtmp3 - dtmp4) .gt. dtmp2) then
        if (lprths) write(6,*)'calchs: failed to find renormalized edot '
        goto 998
     end if
     edotp = dtmp3
     dtmp1 = delta0 / (2.d0 * pi) * (1.d0 / dabs(delta1 - uu * 5.d-1) - 1.d0 / dabs(delta1 + uu * 5.d-1)) 
     dtmp2 = delta0  * 2.d0 / pi  * (1.d0 / dabs(delta1 - uu * 5.d-1) + 1.d0 / dabs(delta1 + uu * 5.d-1)) 
     dtmp3 = dtmp2 / (1.d0 + (pi * dtmp1)**2)
     tkasym = 0.182d0 * dabs(edotp) * dsqrt(dtmp3) * dexp(-1.d0 / dtmp3) 
     if (lprths) then
        write(6,*)
        write(6,*)' Asymmetric Anderson Kondo temperature             = ', tkasym
     end if
!
!     tkasym = tkasym * dexp(pi**2/ 16.d0 * dtmp3)
!     write(6,*)'   With the additional exponent, T_K (energy unit) = ', tkasym
!
     call flush(6)
!  end if
  998 continue
!
!
  eup   = eup   / hbar
  edown = edown / hbar
  uu    = uu    / hbar
  tupdn = tupdn / hbar
!
  cmtmp1(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, 1, 1), 0.d0)  ! c_up
  cmtmp2(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, 1, 2), 0.d0)  ! c_down
!
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp1, nrho, &
             czero, cmtmp3, nrho)
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp2, nrho, cmtmp2, nrho, &
             czero, cmtmp4, nrho)
  hs = hs + eup * cmtmp3 + edown * cmtmp4
!
  call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp3, nrho, cmtmp4, nrho, &
             czero, cmtmp5, nrho)
  hs = hs + uu * cmtmp5
!
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp2, nrho, &
             czero, cmtmp6, nrho)
  do ni=1,nrho
    do nj=1,nrho
      hs(ni,nj) = hs(ni,nj) + (cmtmp6(ni,nj) + dconjg(cmtmp6(nj,ni))) * tupdn
    end do
  end do
end if
!
!--------------------- for (norbs=2, nspin=1) ----------------
if (nspin .eq. 1 .and. norbs .eq. 2) then
  engy01 = 0.d0   ! if doubledot, coupled to left  lead
  engy02 = 0.d0   ! if doubledot, coupled to right lead (if nalf .ne. 1)
  u12    = 0.d0
  t12    = 0.d0
  phi12  = 0.d0
  fixdot = .true.
  rewind(5)
  read(5, para3, end=114)
  114 continue
  if (lprths) then
     write(6,*)
     write(6,*)' parameters for Hamiltonian '
     write(6,115)engy01, engy02, u12, t12, phi12
  end if
  115 format(' e01 = ', f14.8, 2x, ' e02 = ', f14.8, 2x,  &
             ' u12 = ', f14.8, 2x, ' t12 = ', f14.8, ' phi12 = ', f14.8)
  call flush(6)
!
  engy01 = engy01 / hbar
  engy02 = engy02 / hbar
  u12    = u12 / hbar
  t12    = t12 / hbar
  cmtmp1(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, 1, 1), 0.d0) ! c_1
  cmtmp2(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, 2, 1), 0.d0) ! c_2
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp1, nrho, &
             czero, cmtmp3, nrho)
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp2, nrho, cmtmp2, nrho, &
             czero, cmtmp4, nrho)
  hs = hs + engy01 * cmtmp3 + engy02 * cmtmp4
  call zgemm('c', 'n', nrho, nrho, nrho, cdexp( eye*phi12), cmtmp1, nrho, cmtmp2, nrho, &
             czero, cmtmp5, nrho)  ! c1^dag c2
  call zgemm('c', 'n', nrho, nrho, nrho, cdexp(-eye*phi12), cmtmp2, nrho, cmtmp1, nrho, &
             cunity, cmtmp5, nrho) ! c2^dag c1
  hs = hs + t12 * cmtmp5 
  call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp3, nrho, cmtmp4, nrho, &
             czero, cmtmp5, nrho)
  hs = hs + u12 * cmtmp5
!
  if (forcesteady .and. dfieldtype .eq. 2) then
    hs = hs + dedot1 * cmtmp3 +  dedot2 * cmtmp4
  end if
!
end if
!
!--------------------- for (norbs=2, nspin=2) ----------------
if (nspin .eq. 2 .and. norbs .eq. 2) then
  e1up   = 0.d0
  e2up   = 0.d0
  e1down = 0.d0
  e2down = 0.d0
  u12    = 0.d0
  t12    = 0.d0
  uu1    = 0.d0
  uu2    = 0.d0
  j12    = 0.d0
  fixdot = .true.
  rewind(5)
  read(5, para4, end=116)
  116 continue
  if (lprths) then
     write(6,*)
     write(6,*)' parameters for Hamiltonian '
     write(6,117)e1up, e1down, e2up, e2down
     write(6,625)u12, t12, uu1, uu2
  end if
  117 format('  e1up = ', f14.8, 2x, ' e1down = ', f14.8, 2x, ' e2up = ', f14.8, 2x, ' e2down = ', f14.8, 2x)
  625 format('   u12 = ', f14.8, 2x, '    t12 = ', f14.8, 2x, '  uu1 = ', f14.8, 2x, '    uu2 = ', f14.8, 2x)
  call flush(6)
!
! Heisenberg model  H' = J S1.S2 (dot product), J < 0 implies ferromagnetic coupling
!
  if (lspin3d) then
     if (lprths .and. dabs(j12) .gt. dpico) then
        write(6,*)'calchs: Heisenberg coupling (J12 * S1.S2) is invoked ' 
        if (j12 .gt. 0.d0) then
           write(6,*)'calchs: antiferromagnetic coupling J12 = ', j12
        else 
           write(6,*)'calchs: ferromagnetic coupling J12 = ', j12
        end if
     end if
     call flush(6)
  end if
!
  e1up   = e1up   / hbar
  e1down = e1down / hbar
  e2up   = e2up   / hbar
  e2down = e2down / hbar
  u12    = u12    / hbar
  t12    = t12    / hbar
  uu1    = uu1    / hbar
  uu2    = uu2    / hbar
  j12    = j12    / hbar
  cmtmp1(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, 1, 1), 0.d0) ! c1up
  cmtmp2(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, 1, 2), 0.d0) ! c1down
  cmtmp3(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, 2, 1), 0.d0) ! c2up
  cmtmp4(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, 2, 2), 0.d0) ! c2down
!
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp1, nrho, &
             czero, cmtmp5, nrho)                                     ! n1up
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp2, nrho, cmtmp2, nrho, &
             czero, cmtmp6, nrho)                                     ! n1down
  hs = hs + e1up * cmtmp5 + e1down * cmtmp6
!
  call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp5, nrho, cmtmp6, nrho, &
             czero, cmtmp7, nrho)                                     ! n1up*n1down
  hs = hs + uu1 * cmtmp7
  cmtmp7 = cmtmp5 + cmtmp6                                            ! n1
!
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp3, nrho, cmtmp3, nrho, &
             czero, cmtmp5, nrho)                                     ! n2up
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp4, nrho, cmtmp4, nrho, &
             czero, cmtmp6, nrho)                                     ! n2down
  hs = hs + e2up * cmtmp5 + e2down * cmtmp6
!
  call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp5, nrho, cmtmp6, nrho, &
             czero, cmtmp8, nrho)                                     ! n2up*n2down
  hs = hs + uu2 * cmtmp8
  cmtmp8 = cmtmp5 + cmtmp6                                            ! n2
!
  call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp7, nrho, cmtmp8, nrho, &
             czero, cmtmp5, nrho)                                     ! n1*n2
  hs = hs + u12 * cmtmp5                                                
!
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp3, nrho, &
             czero, cmtmp5, nrho)
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp3, nrho, cmtmp1, nrho, &
             cunity, cmtmp5, nrho)
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp2, nrho, cmtmp4, nrho, &
             cunity, cmtmp5, nrho)
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp4, nrho, cmtmp2, nrho, &
             cunity, cmtmp5, nrho)
  hs = hs + t12 * cmtmp5
!
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp1, nrho, &
             czero, cmtmp5, nrho)                                     ! n1up
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp2, nrho, cmtmp2, nrho, &
             czero, cmtmp6, nrho)                                     ! n1down
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp3, nrho, cmtmp3, nrho, &
             czero, cmtmp7, nrho)                                     ! n2up
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp4, nrho, cmtmp4, nrho, &
             czero, cmtmp8, nrho)                                     ! n2down
!
  if (forcesteady .and. dfieldtype .eq. 2) then
    hs = hs + dedot1 * (cmtmp5 + cmtmp6) + dedot2 * (cmtmp7 + cmtmp8)
  end if
!
  if (lspin3d) then
! Heisenberg term
     if (dabs(j12) .gt. dpico) then
        call zgemm('n', 'n', nrho, nrho, nrho, cunity, sopr(1,1,1,1), nrho,  &
                   sopr(1,1,1,2), nrho, czero,  cmtmp1, nrho)
        call zgemm('n', 'n', nrho, nrho, nrho, cunity, sopr(1,1,2,1), nrho,  &
                   sopr(1,1,2,2), nrho, cunity, cmtmp1, nrho)
        call zgemm('n', 'n', nrho, nrho, nrho, cunity, sopr(1,1,3,1), nrho,  &
                   sopr(1,1,3,2), nrho, cunity, cmtmp1, nrho)
! check Hermicity     
        call checkhermicity(cmtmp1, nrho, cmtmp2, nrho, lherm, dtmp1)
        if (.not. lherm) then
           write(6,*)'calchs: error! check Hermicity fails 1 ', dtmp1
           stop
        end if
        hs(1:nrho,1:nrho) = hs(1:nrho,1:nrho) + j12 * cmtmp1(1:nrho,1:nrho)
     end if
! interaction with magnetic field
     if (lbfield3d) then
        cmtmp1(1:nrho,1:nrho) = -dbf3d(1,1) * sopr(1:nrho,1:nrho,1,1)  &
                                -dbf3d(2,1) * sopr(1:nrho,1:nrho,2,1)  &
                                -dbf3d(3,1) * sopr(1:nrho,1:nrho,3,1)  &  
                                -dbf3d(1,2) * sopr(1:nrho,1:nrho,1,2)  &
                                -dbf3d(2,2) * sopr(1:nrho,1:nrho,2,2)  &
                                -dbf3d(3,2) * sopr(1:nrho,1:nrho,3,2)     
        call checkhermicity(cmtmp1, nrho, cmtmp2, nrho, lherm, dtmp1)
        if (.not. lherm) then
           write(6,*)'calchs: error! check Hermicity fails 2 ', dtmp1
           stop
        end if
        hs(1:nrho,1:nrho) = hs(1:nrho,1:nrho) + cmtmp1(1:nrho,1:nrho)
     end if
! zero-field splitting
     if (lzfs) then
        cmtmp1(1:nrho,1:nrho) = sopr(1:nrho,1:nrho,1,1) + sopr(1:nrho,1:nrho,1,2)  ! Sx
        call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp1, nrho, czero, cmtmp2, nrho)  ! Sx^2
        hs(1:nrho,1:nrho) = hs(1:nrho,1:nrho) + d_xx * cmtmp2(1:nrho,1:nrho)
!
        cmtmp1(1:nrho,1:nrho) = sopr(1:nrho,1:nrho,2,1) + sopr(1:nrho,1:nrho,2,2)  ! Sy
        call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp1, nrho, czero, cmtmp2, nrho)  ! Sy^2
        hs(1:nrho,1:nrho) = hs(1:nrho,1:nrho) + d_yy * cmtmp2(1:nrho,1:nrho)
!
        cmtmp1(1:nrho,1:nrho) = sopr(1:nrho,1:nrho,3,1) + sopr(1:nrho,1:nrho,3,2)  ! Sz
        call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp1, nrho, czero, cmtmp2, nrho)  ! Sz^2
        hs(1:nrho,1:nrho) = hs(1:nrho,1:nrho) + d_zz * cmtmp2(1:nrho,1:nrho)
     end if
  end if
end if
!
!--------------------- for (norbs=3, nspin=1) ----------------
if (norbs .eq. 3 .and. nspin .eq. 1) then
  engy01 = 0.d0  
  engy02 = 0.d0  
  engy03 = 0.d0
  t12    = 0.d0
  t23    = 0.d0
  t13    = 0.d0
  fixdot = .true.
  rewind(5)
  read(5, sys3level, end=118)
  118 continue
  if (lprths) then
     write(6,*)
     write(6,*)' parameters for Hamiltonian '
     write(6,119)engy01, engy02, engy03, t12, t23, t13
  end if
  119 format(' e01 = ', f14.8, 2x, ' e02 = ', f14.8, 2x, ' e03 = ', f14.8,  &
             ' t12 = ', f14.8, 2x, ' t23 = ', f14.8, 2x, ' t13 = ', f14.8) 
  call flush(6)
!
  engy01 = engy01 / hbar
  engy02 = engy02 / hbar
  engy03 = engy03 / hbar
  t12    = t12 / hbar
  t23    = t23 / hbar
  t13    = t13 / hbar
  cmtmp1(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, 1, 1), 0.d0) ! c_1
  cmtmp2(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, 2, 1), 0.d0) ! c_2
  cmtmp3(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, 3, 1), 0.d0) ! c_3
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp1, nrho, &
             czero, cmtmp4, nrho)
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp2, nrho, cmtmp2, nrho, &
             czero, cmtmp5, nrho)
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp3, nrho, cmtmp3, nrho, &
             czero, cmtmp6, nrho)
  hs = hs + engy01 * cmtmp4 + engy02 * cmtmp5 + engy03 * cmtmp6
!
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp2, nrho, &
             czero, cmtmp4, nrho)   ! c_1^dag c_2
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp2, nrho, cmtmp1, nrho, &
             cunity, cmtmp4, nrho)  ! c_2^dag c_1
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp2, nrho, cmtmp3, nrho, &
             czero, cmtmp5, nrho)   ! c_2^dag c_3
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp3, nrho, cmtmp2, nrho, &
             cunity, cmtmp5, nrho)  ! c_3^dag c_2
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp3, nrho, &
             czero, cmtmp6, nrho)   ! c_1^dag c_3
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp3, nrho, cmtmp1, nrho, &
             cunity, cmtmp6, nrho)  ! c_3^dag c_1
  hs = hs + t12 * cmtmp4 + t23 * cmtmp5 + t13 * cmtmp6
end if
!
999 continue
!
if (lprths .and. .not. lhf) then
   dmtmp1(1:nrho,1:nrho) = dble(hs(1:nrho,1:nrho)) * hbar
   write(6,*)
   write(6,*)'calchs: real part of equilibrium Hamiltonian '
   call amatout(nrho, nrho, dmtmp1)
   if (fixdot) then
     write(6,*)'calchs: fixdot=TRUE The system Hamiltonian is fixed during time evolution '
   end if
end if
call flush(6)
!
lprths = .false.   ! by default, this subroutine executes in silent mode
!
return
end subroutine calchs
