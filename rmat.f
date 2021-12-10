      module rmat
      implicit none
      real*8,dimension(:),allocatable :: nfc,ngc,nfcp,ngcp ! used for coul90
      contains


      subroutine scatt2b()
      use channels
      use lagrange_mesh_single_channel
      use mesh
      use coulfunc
      use input
      use constants
      use precision
      use pot
      implicit none
      real*8 :: eta, kr
      integer :: ifail
      integer :: nch, l
      real*8 :: s,ls,j
      real*8 :: hm
      complex*16,dimension(1:nr) :: wf
      complex*16 :: smat
      integer :: ir
      integer :: ie
      real*8,dimension(0:lmax) :: cph  !Coulomb phase-shift
      complex*16,dimension(1:nr) :: X_vector


      !initial Lagrange mesh
      call initial_lagrange_func(rmax)
      call T_and_Bloch(mu)
      allocate(nfc(0:lmax),ngc(0:lmax),nfcp(0:lmax),ngcp(0:lmax))

      ! compute the boundary conditions
      eta=za*zb*e2*mu/hbarc/hbarc/k
      kr=k*rmax
      call coul90(kr,eta,zero,lmax,nfc,ngc,nfcp,ngcp,0,ifail)
      if (ifail/=0) then
      write(*,*) 'coul90: ifail=',ifail; stop
      endif

      ! compute the Coulomb phase-shift
      call coulph(eta,cph,lmax)

      do nch=1, alpha2b%nchmax
        l=alpha2b%l(nch)
        s=alpha2b%s(nch)
        j=alpha2b%j(nch)
        ls=0.5_dpreal*(j*(j+1)-l*(l+1)-s*(s+1))
         call potr(za*zb,ls)


C         if(nch==1)  call R_matrix(l,mu,ecm,v,cph(l),ngc(l),ngcp(l),nfc(l),nfcp(l),smat,wf,X_vector)
         call R_matrix_iov(l,mu,ecm,v,cph(l),ngc(l),ngcp(l),nfc(l),nfcp(l),smat,wf)


C        call rmat_inho(nr,rmax,v,ecm,eta,hm,l,smat,wf)

        write(*,100) l,s,j,real(smat),aimag(smat)
        write(32,101)l,s,j
        write(42,101)l,s,j
        write(52,101)l,s,j
        do ir=1,nr
           write(32,*) rr(ir), real(wf(ir)),aimag(wf(ir))
           write(42,*)rr(ir), abs(wf(ir))
           write(52,*) rr(ir), real(wf(ir))
        end do



      end do


      deallocate(nfc,ngc,nfcp,ngcp)


100   format('l=',I3,' s=',f3.1, ' j=', f7.1, ' s-mat=(',2f10.6,')')
101   format('&l=',I3,' s=',f3.1, ' j=', f5.1)

C     end do   ! add loop

      end subroutine




      end module
