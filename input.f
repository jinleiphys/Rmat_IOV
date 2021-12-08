      module input
      implicit none
      real*8 :: Ecm, mu ! energy in CM frame and reduce mass : unit in MeV
      real*8 ::  k !  momentum
C     integer :: method ! method=1 lagrange mesh method, /=1 numerov method
      real*8 :: massa, massb !mass munber of interaction partciles
      real*8 :: za, zb ! charge of interaction partciles

      logical,dimension(9999) :: written

      contains

      subroutine readin()
      use mesh
      use precision
      use constants
      use channels
      use pot
      use gauss
      implicit none
      real*8 :: Elab



      namelist /global/ Elab,hcm,rmax,lmin,lmax,nr,maxiter
      namelist /systems/ massa, massb, za,zb,ja,jb

      namelist /potential/ ptype,a1,a2,rc,uv,av,rv,uw,aw,rw,
     +                     vsov,rsov,asov,vsow,rsow,asow,vd,avd,rvd,wd,awd,rwd


      hcm=0.05;rmax=40;lmin=0;lmax=0;nr=60;maxiter=20
      read(5,nml=global)
      irmax=nint(rmax/hcm)




      massa=1.0; massb=1.0; za=0.0;zb=0.0;ja=0.0;jb=0.0
      read(5,nml=systems)
      mu=massa*massb*amu/(massa+massb)
      write(*,*) "mu=",mu
      ecm=elab*massb/(massa+massb)
      write(*,*) "ecm=",ecm

      k=sqrt(2.0_dpreal*mu*ecm)/hbarc
      call alpha_2b()


       ptype=1;a1=0.0d0;a2=massb;rc=1.3d0
       uv=0.0d0;av=0.1d0;rv=0.0d0
       uw=0.0d0;aw=0.1d0;rw=0.0d0
       vsov=0.0d0;rsov=0.0d0;asov=0.1d0
       vsow=0.0d0;rsow=0.0d0;asow=0.1d0
       vd=0.0d0;avd=0.1d0;rvd=0.0d0
       wd=0.0d0;awd=0.1d0;rwd=0.0d0
       read(5,nml=potential)




      end subroutine





      subroutine fkind()
       character*40 flkind(9999)
       integer writf(9999),nwrit,i
       flkind(:) = ' '

       flkind(8)='channels coupling index'
       written(8)=.TRUE.


       flkind(30)='potential used in the calcuation'
       written(30)=.TRUE.



       flkind(32)='solutions of r-matrix'
       written(32)=.TRUE.



       flkind(42)='module solutions of r-matrix'
       written(42)=.TRUE.


       nwrit = 0
       do i=1,9999
        if(written(i)) then
        flkind(i) = trim(flkind(i))//'.'
        nwrit = nwrit+1
        writf(nwrit) = i
        endif
       enddo
       write(*,990) (writf(i),flkind(writf(i)),i=1,nwrit)
990    format(/'  The following files have been created:',
     X  /(2x,2(i3,':',1x,a40)))
       return
       end subroutine



      end module
