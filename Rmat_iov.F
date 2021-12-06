      program Rmat_IOV
      use input
      use rmat
      use mesh
      use gauss
      use precision
      implicit none
      real*8 ::  t1, t2 , tr, tg, tn
      integer :: ir

      call readin()

      ! use Lagrange mesh method to solve the inhomogeneous equation
C     do nr=15,80,5
C     write(99,*) "&nr=",nr
C     write(*,*) "&nr=",nr
C     do ir=100, 400, 1
      rmax=ir*0.1_dpreal
      if (allocated(rr)) deallocate(rr, rrw)
      allocate(rr(1:nr),rrw(1:nr))
      call gauleg(nr,0.0_dpreal,rmax,rr,rrw)

      call cpu_time(t1)
      call Inhomo_Rmatrix_method()
      call cpu_time(t2)
      tr=t2-t1


C     end do
C     end do
C
C     stop




      call cpu_time(t1)
C      call Inhomo_greenf_method()
      call cpu_time(t2)
      tg=t2-t1
      deallocate(rr,rrw)





      call fkind()

      write(*,*) "R-matrix takes time: ", tr
      write(*,*) "R-matrix_IOV takes time: ", tg


      end program
