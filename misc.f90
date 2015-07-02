module misc
      implicit none
      save
      contains

      subroutine check_mpi_error(ierr)
            use mpi
            use iso_fortran_env

            integer :: ierr
            integer :: double_err

            integer :: length
            character (len=256) :: err_str

            if (ierr .ne. MPI_SUCCESS) then
                  call MPI_ERROR_STRING(ierr, err_str, length, double_err)

                  if (double_err .ne. MPI_SUCCESS) then
                        write(error_unit,*) "Detected an error and failed to retrive error string"
                        stop
                  endif

                  write(error_unit,*) err_str(:length)
                  stop
            endif

      end subroutine check_mpi_error

!!!!!!!!!RANDOM NUMBER GENERATOR!!!!!!!!!!!!!!

      real function pad_ranf()
            implicit none
            call random_number(pad_ranf)
      end function pad_ranf



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine seed_mpi(my_rank)
            integer, intent(in):: my_rank
            integer:: time(8)
            integer, dimension(12):: seed

            call date_and_time(values=time)

            seed(:) = time(4) * ( 360000*time(5) + 6000*time(6) + 100*time(7) + time(8)) + my_rank*100
            call random_seed(PUT=seed)

      end subroutine seed_mpi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_bndry_Eflux(b1,E,bndry_Eflux)
            use dimensions
            use inputs, only: mO,q,dt,dx,dy,km_to_m,mu0
            implicit none
            real, intent(in):: b1(nx,ny,nz,3), E(nx,ny,nz,3)
            real, intent(inout):: bndry_Eflux
            integer:: i,j,k
            real:: exb_flux, mO_q

            mO_q = mO/q
            !k=2 face
            do i = 2, nx
                  do j= 2, ny
!                        m=3
                        k=2
                        exb_flux = (mO_q)**2*(1.0/mu0)*dt*dx*dy* &
                              (E(i,j,k,1)*b1(i,j,k,2) - E(i,j,k,2)*b1(i,j,k,1))* &
                              km_to_m**3
                        bndry_Eflux = bndry_Eflux + exb_flux

                  enddo
            enddo

            !k=nx face

            do i =2,nx
                  do j=2,ny
                        k = nz-1
                        exb_flux = (mO_q)**2*(1.0/mu0)*dt*dx*dy* &
                              (E(i,j,k,1)*b1(i,j,k,2) - E(i,j,k,2)*b1(i,j,k,1))* &
                              km_to_m**3
                        bndry_Eflux = bndry_Eflux - exb_flux

                  enddo
            enddo

      end subroutine get_bndry_Eflux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_beta(Ni_tot_sys,beta)
            use dimensions
            use grid, only: qx,qy,qz
            use inputs, only: np_top
            implicit none
            integer(4), intent(in):: Ni_tot_sys
            real, intent(out):: beta
            real:: vol

            vol = ((qx(nx-1)-qx(1))*(qy(ny-1)-qy(1))*(qz(nz)-qz(1)))
            beta = (Ni_tot_sys/vol)/np_top

            write(*,*) 'beta....',beta

      end subroutine get_beta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module misc
