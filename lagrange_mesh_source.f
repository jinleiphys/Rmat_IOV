      module  lagrange_mesh_single_channel
      use mesh
      implicit none
      ! note the delta function structure fi(xj) = \lambda_i^{-1/2} \delta_{ij}
      real*8,dimension(:),allocatable :: lagrange_func
      real*8,dimension(:,:), allocatable :: kinetic_matrix
      real*8,dimension(:,:),allocatable :: l_barrier_matrix
      complex*16,dimension(:,:),allocatable :: V_matrix
      real*8,dimension(:),allocatable ::  mesh_rr, mesh_rw
      real*8 :: rmax_rmatrix
      contains

c***********************************************************************
       subroutine initial_lagrange_func (rmatch_rmatrix)
       !
       ! before calling this subroutine please initial rmax_rmatrix and nr
       !
       ! this subroutine is used to initial the Lagrange function
       ! for the purpose to reconstruct the scattering wave function

       !!!!!! please gives the matching radius
c***********************************************************************
       use gauss
       use mesh
       use precision
       implicit none
       integer :: ir
       real*8 :: rmatch_rmatrix
       rmax_rmatrix=rmatch_rmatrix

       if(allocated(lagrange_func)) deallocate(lagrange_func)
       if(allocated(mesh_rr)) deallocate(mesh_rr)
       if(allocated(mesh_rw)) deallocate(mesh_rw)

       allocate(lagrange_func(1:nr))
       allocate(mesh_rr(1:nr))
       allocate(mesh_rw(1:nr))

       call gauleg(nr,0.0_dpreal,1.0_dpreal,mesh_rr,mesh_rw)

       do ir=1, nr
         lagrange_func(ir) = 1/ sqrt(rmax_rmatrix * mesh_rw(ir))
       end do

       end subroutine
c-----------------------------------------------------------------------
c***********************************************************************
       subroutine T_and_Bloch(mu)
       ! this subroutine is used to compute the kinetic and Bloch operator
       ! in the matrix form
c***********************************************************************
       use precision
       use constants
       implicit none
       integer :: ir, irp
       real*8 :: f1,f2,f3
       real*8 ::  xi, xj
       real*8 :: mu


       if(allocated(kinetic_matrix)) deallocate(kinetic_matrix)
       allocate(kinetic_matrix(1:nr, 1:nr))

       do ir=1 , nr
         do irp=1, nr

            f1=0.0_dpreal
            f2=0.0_dpreal
            f3=0.0_dpreal
            xi= mesh_rr(ir)
            xj= mesh_rr(irp)

            if(ir==irp) then
               f1= hbarc**2 / (6.0_dpreal * rmax_rmatrix**2 * xi * (1-xi) * mu)
               f2= 4.0_dpreal * nr * (nr+1) + 3.0
               f3= (1.0_dpreal - 6.0_dpreal *xi )/( xi * (1-xi) )
               kinetic_matrix(ir,irp)= f1 * ( f2 + f3 )
            else
               f1= (-1.0_dpreal)**(ir+irp) * hbarc**2 /2.0_dpreal/rmax_rmatrix/rmax_rmatrix/
     +             sqrt( xi * xj * (1-xi) * (1-xj) ) / mu
                f2= nr * (nr+1)  + 1.0 + ( xi + xj - 2.0_dpreal * xi * xj )/( (xi-xj)**2 )
              f3= 1.0_dpreal/(1.0_dpreal-xi)  + 1.0_dpreal/(1.0_dpreal -xj)


               kinetic_matrix(ir,irp)= f1* (f2-f3)
            end if
         end do
       end do



       end subroutine
c-----------------------------------------------------------------------
c***********************************************************************
       subroutine centrifugal_barrier (l,mu)
       ! this subroutine is used to compute centrifugal_barrier in
       ! lagrange mesh basis
       ! V_l = \frac{\hbarc^2 l (l+1)}{2 * mu * r}
       ! input : mu  in MeV
       !         l
c***********************************************************************
       use constants
       use precision
       implicit none
       integer :: l , ir
       real*8 :: mu ,xi
       if(allocated(l_barrier_matrix)) deallocate(l_barrier_matrix)
       allocate(l_barrier_matrix(1:nr, 1:nr))
       l_barrier_matrix=0.0_dpreal

       do ir=1 , nr
         xi = mesh_rr(ir)
         l_barrier_matrix(ir,ir) = ( hbarc**2 * l * (l+1) )/ ( 2.0_dpreal * mu * rmax_rmatrix * rmax_rmatrix * xi * xi  )
       end do

       end subroutine
c-----------------------------------------------------------------------
c***********************************************************************
       subroutine lagrange_V(Vpot)
       ! this subroutine is used to compute the V-matrix in the Lagrange
       ! mesh basis
       ! input V is given in the mesh point in the step of hcm
       ! this subroutine interpolate the Vpot to give the correct mesh sets
c***********************************************************************
       use interpolation
       use precision
       implicit none
       complex*16,dimension(1:nr),intent(in) :: vpot
       integer :: ir
       real*8 ::  r

       if(allocated(V_matrix)) deallocate(V_matrix)
       allocate(V_matrix(1:nr,1:nr))
       V_matrix=0.0_dpreal

       do ir = 1, nr
         r = rmax_rmatrix*  mesh_rr(ir)
         V_matrix(ir,ir) = vpot(ir)
       end do


       end subroutine
c-----------------------------------------------------------------------
c***********************************************************************
       subroutine R_matrix(l,mu,ecm,vpot,cph,gc,gcp,fc,fcp,Smatrix,scattwf)
       ! this subroutine is used to compute the single channel R-matrix
c***********************************************************************
       use precision
       use constants
       implicit none
       integer :: l
       real*8 :: mu, ecm ,k
       real*8 :: gc,gcp,fc,fcp
       complex*16,dimension(1:nr) :: vpot
       complex*16,dimension(1:nr,1:nr) :: C_minus_E
       complex*16,dimension(1:nr) ::  B_vector , X_vector
       integer :: ir,INFO,NRHS
       integer,dimension(1:nr) :: IPIV
       real*8 :: N
       complex*16 ::  Rmatrix , Zmatrix , Smatrix, Zmatrix_O, Zmatrix_I
       complex*16 :: hc ,hc1 !H(+), H(-)
       complex*16 :: hcp,hcp1 ! derivatives of H(+),H(-)
       complex*16,dimension(1:nr) :: scattwf
       complex*16 :: wf_bound
       real*8 ::  cph
       integer :: irp
       real*8 :: delta


       k=sqrt(2.*mu*ecm/(hbarc**2))
       hc=cmplx(gc,fc,kind=8)
       hc1=cmplx(gc,-fc,kind=8)
       hcp=cmplx(gcp,fcp,kind=8)
       hcp1=cmplx(gcp,-fcp,kind=8)

       call lagrange_V(Vpot)
       call centrifugal_barrier(l,mu)

       C_minus_E=0.0_dpreal
       B_vector=0.0_dpreal   ! fn(a) in sofia's notes
       C_minus_E = kinetic_matrix + l_barrier_matrix + V_matrix

       do ir=1, nr
          C_minus_E(ir,ir) =  C_minus_E(ir,ir) - ecm
          B_vector(ir) = (-1.0_dpreal)**ir / sqrt( rmax_rmatrix* mesh_rr(ir) * ( 1.0-mesh_rr(ir) ) )   ! calculate f_n^(a) in the notes
       end do

       NRHS=1
       X_vector= B_vector



       call ZGESV( nr, NRHS, C_minus_E, nr, IPIV, X_vector, nr, INFO )
       If(INFO/=0) stop "error in calling ZGESV"


       do ir=1,nr
         write(201,*) mesh_rr(ir),abs(X_vector(ir))

       end do

       write(201,*) "&"


       ! normalize factor
       N = hbarc**2 / (2.0_dpreal * mu * rmax_rmatrix)

       Rmatrix=0.0_dpreal
       do ir=1,nr

            Rmatrix = Rmatrix + B_vector(ir) * X_vector(ir)
       end do

       Rmatrix=Rmatrix * N

       Zmatrix= (hc-k*rmax_rmatrix*Rmatrix*hcp) / (k*rmax_rmatrix)

       Zmatrix_O=Zmatrix
       Zmatrix_I= (hc1-k*rmax_rmatrix*Rmatrix*hcp1) / (k*rmax_rmatrix)

C      Smatrix= CONJG(Zmatrix) / Zmatrix

        Smatrix= Zmatrix_I / Zmatrix_O

       delta=0.5d0*ACOS(real(Smatrix))
       if (aimag(Smatrix)<0) delta=-delta+pi
       DELTA=DELTA*180.0d0/pi

C       write(*,*) "l=",l,"Smatrix=", Smatrix


       ! compute the scattering wave function  ! still testing

       wf_bound = 0.5* iu * (hc1-hc*Smatrix) * exp(iu*cph)
       do ir=1,nr

        scattwf(ir)=lagrange_func(ir)*X_vector(ir)*wf_bound *N / Rmatrix

       write(99,*) mesh_rr(ir)* rmax_rmatrix , real(scattwf(ir)), aimag(scattwf(ir))
       end do

       write(99,*) "&"


       end subroutine
c-----------------------------------------------------------------------


c***********************************************************************
       subroutine R_matrix_iov(l,mu,ecm,vpot,cph,gc,gcp,fc,fcp,Smatrix,scattwf)
       ! this subroutine is used to compute the single channel R-matrix
c***********************************************************************
       use precision
       use constants
       implicit none
       integer :: l
       real*8 :: mu, ecm ,k
       real*8 :: gc,gcp,fc,fcp
       complex*16,dimension(1:nr) :: vpot
       complex*16,dimension(1:nr,1:nr) :: C_minus_E
       complex*16,dimension(1:nr) ::  B_vector,X_vector,X_vector_store,xx
       complex*16,dimension(1:maxiter,1:maxiter) :: A_matrix
       complex*16,dimension(1:maxiter,1:maxiter) :: A1_matrix
       complex*16,target,dimension(1:nr,1:maxiter+1) :: eta_bar
       complex*16,dimension(1:nr) :: eta1, eta_tilde,f0
       complex*16,dimension(:),pointer :: eta_bar1
       integer :: ir,INFO,NRHS
       complex*16,dimension(1:maxiter) :: etapsi,af
       integer, dimension(1:maxiter) :: IPIV
       real*8 :: N
       complex*16 ::  Rmatrix , Zmatrix , Smatrix, Zmatrix_O, Zmatrix_I
       complex*16 :: hc ,hc1 !H(+), H(-)
       complex*16 :: hcp,hcp1 ! derivatives of H(+),H(-)
       complex*16,dimension(1:nr) :: scattwf
       complex*16 :: wf_bound
       real*8 ::  cph
       integer :: irp
       real*8 :: delta
       complex*16 :: dotprot,mod,mod1
       real*8 :: abs_eta_tilde,abs_eta_tilde1
       integer :: i,j


       k=sqrt(2.*mu*ecm/(hbarc**2))
       hc=cmplx(gc,fc,kind=8)
       hc1=cmplx(gc,-fc,kind=8)
       hcp=cmplx(gcp,fcp,kind=8)
       hcp1=cmplx(gcp,-fcp,kind=8)

       call lagrange_V(Vpot)
       call centrifugal_barrier(l,mu)

       C_minus_E=0.0_dpreal
       B_vector=0.0_dpreal   ! fn(a) in sofia's notes
       C_minus_E = kinetic_matrix + l_barrier_matrix + V_matrix

       do ir=1, nr
          C_minus_E(ir,ir) =  C_minus_E(ir,ir) - ecm
          B_vector(ir) = (-1.0_dpreal)**ir / sqrt( rmax_rmatrix* mesh_rr(ir) * ( 1.0-mesh_rr(ir) ) )   ! calculate f_n^(a) in the notes
       end do

       NRHS=1
C       X_vector= B_vector


! begin IOV method
      eta_bar(:,1)=B_vector-matmul(C_minus_E,B_vector)
C      eta_bar(:,1)=1.0_dpreal
      eta_bar1=>eta_bar(1:nr,1)
      call COMDOTPRODUCT(eta_bar1,eta_bar1,dotprot)
      abs_eta_tilde=abs(dotprot)
      abs_eta_tilde=sqrt(abs_eta_tilde)
      eta_bar1=eta_bar1/abs_eta_tilde
      etapsi=0.0_dpreal

      do i=1, maxiter
        eta_bar1=>eta_bar(1:nr,i)
        call COMDOTPRODUCT(eta_bar1,B_vector,etapsi(i))
        ! calculate eta_{i+1}
        eta1=matmul(C_minus_E,eta_bar1)

        do j=1,nr
           write(102,*)mesh_rr(j),abs(eta1(j))
        end do
        write(102,*)"&"


        !Orthogonalize eta_{i+1} => eta_tilde
        eta_tilde=eta1
C        f0=0.0_dpreal
        do j=1, i
          eta_bar1=>eta_bar(1:nr,j)
          call COMDOTPRODUCT(eta_bar1,eta1,A_matrix(j,i))
C          f0=f0+A_matrix(j,i)* eta_bar1
          eta_tilde=eta_tilde-A_matrix(j,i)* eta_bar1
        end do
C        eta_tilde=eta1-f0


        call COMDOTPRODUCT(eta_tilde,eta_tilde,dotprot)
        abs_eta_tilde=abs(dotprot)
        abs_eta_tilde=sqrt(abs_eta_tilde)
        eta_bar1=>eta_bar(1:nr,i+1)
        eta_bar1=eta_tilde/abs_eta_tilde

        !!!!store value of abs_eta_tilde
        if (i/=1) A_matrix(i,i-1) = abs_eta_tilde1
        abs_eta_tilde1=abs_eta_tilde

        if (i>1) then
         A1_matrix=A_matrix
         af=etapsi ! a coefficient
         call ZGESV(i,NRHS,A1_matrix,maxiter,IPIV,af,maxiter,INFO )
         If(INFO/=0) stop "error in calling ZGESV"
         X_vector=matmul(eta_bar,af)
         call COMDOTPRODUCT(X_vector,X_vector,mod)
        end if

C        write(*,*) "i=",i,"mod=",mod

C        do j=1,nr
C           write(103,*)mesh_rr(j),abs(eta_bar(j,i))
C        end do
C        write(103,*)"&"
C        stop
C         write(*,*) "abs(mod)=",abs(mod), "abs(mod1)=",abs(mod1)

C        if (i>1) write(*,*) "abs(abs(mod)-abs(mod1))=",abs(abs(mod)-abs(mod1))
        if (i/=1 .and. abs(abs(mod)-abs(mod1)) < 1e-9) exit

        mod1=mod
       end do
        write(*,*) "i=",i

C        xx=X_vector

! finish IOV method



       ! normalize factor
       N = hbarc**2 / (2.0_dpreal * mu * rmax_rmatrix)

       Rmatrix=0.0_dpreal
       do ir=1,nr

            Rmatrix = Rmatrix + B_vector(ir) * X_vector(ir)
       end do

       Rmatrix=Rmatrix * N

       Zmatrix= (hc-k*rmax_rmatrix*Rmatrix*hcp) / (k*rmax_rmatrix)

       Zmatrix_O=Zmatrix
       Zmatrix_I= (hc1-k*rmax_rmatrix*Rmatrix*hcp1) / (k*rmax_rmatrix)

C      Smatrix= CONJG(Zmatrix) / Zmatrix

        Smatrix= Zmatrix_I / Zmatrix_O

       delta=0.5d0*ACOS(real(Smatrix))
       if (aimag(Smatrix)<0) delta=-delta+pi
       DELTA=DELTA*180.0d0/pi

C       write(*,*) "l=",l,"Smatrix=", Smatrix


       ! compute the scattering wave function  ! still testing

       wf_bound = 0.5* iu * (hc1-hc*Smatrix) * exp(iu*cph)
       do ir=1,nr

        scattwf(ir)=lagrange_func(ir)*X_vector(ir)*wf_bound *N / Rmatrix

       write(98,*) mesh_rr(ir)* rmax_rmatrix , real(scattwf(ir)), aimag(scattwf(ir))
       end do

       write(98,*) "&"


C       stop


       end subroutine
c-------------------------------------------------------------------------------
       subroutine COMDOTPRODUCT(eta_bar1,w,product)
 ! this subroutine is used to calculate the dot product of eta_bar1^T W
        implicit none
        complex*16, dimension(1:nr) :: eta_bar1,w
        complex*16 :: product
        integer :: i
        product=0.0d0
        do i=1,nr
           product=product+CONJG(eta_bar1(i))*w(i)
        end do
       end subroutine
c-------------------------------------------------------------------------------


      end module
