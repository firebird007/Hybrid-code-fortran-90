module boundary
      implicit none
      contains

      subroutine periodic(b)
            use dimensions
            implicit none
            real, intent(inout):: b(nx,ny,nz,3)
            integer:: i,j,k,m

!       X direction

            do j=1,ny
                  do k=1,nz
                        do m=1,3
                              b(1,j,k,m) = b(nx-1,j,k,m)
                              b(nx,j,k,m) = b(2,j,k,m)
                        enddo
                  enddo
            enddo

!       Y direction

            do i=1,nx
                  do k=1,nz
                        do m=1,3
                              b(i,1,k,m) = b(i,ny-1,k,m)
                              b(i,ny,k,m) = b(i,2,k,m)
                        enddo
                  enddo
            enddo

!       Z direction

            do i=1,nx
                  do j=1,ny
                        do m=1,3
                              b(i,j,1,m) = b(i,j,nz-1,m)
                              b(i,j,nz,m) = b(i,j,2,m)
                        enddo
                  enddo
            enddo

      end subroutine periodic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine periodic_scalar(b)
            use dimensions
            implicit none
            real, intent(inout):: b(nx,ny,nz)
            integer:: i,j,k

!       X surfaces
            do j=1,ny
                  do k=1,nz
                        b(1,j,k) = b(nx-1,j,k)
                        b(nx,j,k) = b(2,j,k)
                  enddo
            enddo

!       Y surfaces
            do i=1,nx
                  do k=1,nz
                        b(i,1,k) = b(i,ny-1,k)
                        b(i,ny,k) = b(i,2,k)
                  enddo
            enddo

!       Z surfaces
            do i=1,nx
                  do j=1,ny
                        b(i,j,1) = b(i,j,nz-1)
                        b(i,j,nz) = b(i,j,2)
                  enddo
            enddo

      end subroutine periodic_scalar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module boundary

