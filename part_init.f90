module part_init
      implicit none
      save
      contains

      subroutine Energy_diag(Evp,Euf,EB1,EB1x,EB1y,EB1z,EE,EeP)
            use dimensions
            use mpi
            use mult_proc, only: my_rank
            use grid, only: dx_cell,dy_cell,dz_cell
            use inputs, only: mion,q,mu0,mO,km_to_m,epsilon
            use var_arrays, only: vp,b0,b1,E,nu,up,np,Ni_tot,beta,beta_p,input_E,bndry_Eflux,m_arr
            implicit none
            real, intent(out):: Euf,EB1,EB1x,EB1y,EB1z,EE,EeP,Evp
            real:: denf,m_q,recvbuf,total_E,vol
            real:: S_Evp,S_input_E
            integer:: count, ierr
            integer:: i,j,k,m,l

            count = 1
            m_q = mion/q

            Euf = 0.0
            EB1 = 0.0
            EB1x = 0.0
            EB1y = 0.0
            EB1z = 0.0
            EE = 0.0
            EeP = 0.0

            do i=1,nx-1
!                  do j = 1,ny-1
                   j=2
                        do k=1,nz-1
                              vol= dx_cell(i)*dy_cell(j)*dz_cell(k)*km_to_m**3
                              EB1x=EB1x + (vol/(2.0*mu0))*(m_q*b1(i,j,k,1))**2
                              EB1y=EB1y + (vol/(2.0*mu0))*(m_q*b1(i,j,k,2))**2
                              EB1z=EB1z + (vol/(2.0*mu0))*(m_q*b1(i,j,k,3))**2
                                    do m=1,3
                                          denf = np(i,j,k)/(km_to_m**3)
                                          Euf = Euf + 0.5*mO*denf*vol*(up(i,j,k,m)*km_to_m)**2
                                          EB1 = EB1 + (vol/(2.0*mu0))*(m_q*(b1(i,j,k,m)-b0(i,j,k,m)))**2
                                          EE = EE + (epsilon*vol/2.0)*(m_q*E(i,j,k,m)*km_to_m)**2
                                    enddo
                        enddo
!                  enddo
            enddo

            Evp = 0.0
            do l=1, Ni_tot
                  do m=1,3
                        Evp = Evp + 0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 / beta*beta_p(l)
                  enddo
            enddo


            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

            call MPI_ALLREDUCE(Evp,recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            S_Evp = recvbuf

            call MPI_ALLREDUCE(input_E,recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            S_input_E = recvbuf

            total_E = S_Evp + EE + EB1


            if (my_rank .eq. 0) then
                  write(*,*) 'Normalized energy.................',total_E/S_input_E, my_rank
                  write(*,*) 'Normalized energy (bndry).........',total_E/(S_input_E+bndry_Eflux)
            endif


      end subroutine Energy_diag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine load_Maxwellian(vth,Ni_tot_1,mass,mratio)
            use dimensions
            use misc
            use inputs, only: PI, dx, dy, km_to_m
            use grid, only: qx,qy,qz,dz_grid
            use gutsp
            use var_arrays, only: np,vp,vp1,xp,input_p,up,Ni_tot,input_E,ijkp,m_arr,mrat,beta,beta_p,wght,grav,temp_p
            implicit none
            integer(4), intent(in):: Ni_tot_1
            real, intent(in):: mratio, mass, vth

            real:: vth2, vx, vy, vz
            integer:: l,m,i,j,k

!            v1=1.0

            do l = Ni_tot_1,Ni_tot
                  xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
                  xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
                  xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
                  m_arr(l) = mass
                  mrat(l) = mratio
                  beta_p(l) = 1
                  !beta_particle + &       ! Comment for no density gradient
                   !     amp*(1-exp(-((xp(l,3)-qz(nz/2-disp))/ &         !Gaussian distribution
                    !    (grad*dz_grid(nz/2-disp)))**2))
!     x     (.5*amp*(-tanh((xp(l,3)-qz(nz/2-disp))                      !hyperbolic tangent density
!     x     /(grad*dz_grid(nz/2-disp)))
!     x     +tanh((xp(l,3)-qz(nz/2+disp))/
!     x     (grad*dz_grid(nz/2+disp))))+1.0) !Weighting function hyperbolic tangent
!          write(*,*) beta_p(l), 'beta'
!         ijkp(l,1) = floor(xp(l,1)/dx)
!!!!!!!!!!!!!Get P-index!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  i=1
!                  do
!                        if (xp(l,1) .le. qx(i)) exit
!                        i = i+1
!                  enddo
!                  i=i-1

!                  ijkp(l,1) = i
!                  ijkp(l,2) = floor(xp(l,2)/dy)

!                  k=1
!                  do
!                        if (xp(l,3) .le. qz(k)) exit
!                        k=k+1
!                  enddo
!                  k=k-1

!                  ijkp(l,3) = k
!!!!!!!!!!!!!End get P-index!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
                  call get_pindex(i,j,k,l)
!                  vth2=sqrt(vth*vth*beta_p(l)) !thermal speed dependent on np to set up pressure balance for density gradient

                  vx = vth*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf()) !remember to add in vsw to get the flow velocity
                  vy = vth*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  vz = vth*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())

!                  ii = ijkp(l,1)
!                  kk = ijkp(l,3)

!                  vp(l,1) = -0.0*(exp(-(xp(l,3)-qz(nz/2))**2/(10.*delz)**2)
!               x        *exp(-(xp(l,1)-qx(nx/2))**2/(10.*dx)**2))+vx
                  vp(l,1) = vx+57.0*exp(-(xp(l,3)-qz(nz/2))**2/(5*dz_grid(nz/2))**2) !Gaussian velocity perturbation
                  vp(l,2) = vy
                  vp(l,3) = vz

                  do m=1,3
                        vp1(l,m) = vp(l,m)
                        input_E = input_E + 0.5*m_arr(l)*(vp(l,m)*km_to_m)**2/beta * beta_p(l)
                        input_p(m) = input_p(m) + m_arr(l) * vp(l,m) / beta * beta_p(l)
                  enddo

            enddo

            call get_interp_weights()
            call update_np()
            call update_up(vp)

            ! Add a centrifugal gravity term to keep the plasma confined to the torus.  Use T * dn/dz = nmg.
            ! Depends on the density gradient.  Currently set as a gaussian.

!            write(*,*) 'vth.................', vth
!            write(*,*) 'boltzman............', kboltz
!            write(*,*) 'temperature(analytic)..', Temp
            call get_temperature()
!            write(*,*) 'temperature (2,2,100)..', temp_p(2,2,2:10)/1.6e-19
!            stop

            do i=1,nx
            do j=1,ny
            do k=1,nz
!                  grav(i,j,k) = amp*2/(grad*dz_grid(nz/2-disp))**2*(qz(nz/2-disp)-qz(k))* &
!                        exp(-((qz(k)-qz(nz/2-disp))/(grad*dz_grid(nz/2-disp)))**2) * Temp / np(2,2,k)
                  !grav(i,j,k) = -2*Tempcalc/(mion*(grad*dz_grid(nz/2-disp))**2)*(qz(k)-qz(nz/2-disp))
                  grav(i,j,k) = 0
!                  write(*,*) 'gravity.....', grav(i,j,k), i,j,k
!                  grav(i,j,k) = 8.0;
            enddo
            enddo
            enddo
      end subroutine load_Maxwellian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module part_init
