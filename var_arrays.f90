module Var_Arrays
      use dimensions
      implicit none
      save
      real::      b0(nx,ny,nz,3), &     !ambient mag field
                  b1(nx,ny,nz,3), &     !1st order mag field
                  b12(nx,ny,nz,3), &    !b1 at previous time step
                  b1p2(nx,ny,nz,3), &   !temp b1 at time level m+1
                  bt(nx,ny,nz,3), &     !total mag field, mc covarient
                  btc(nx,ny,nz,3), &    !btmf at cell center for particle move
                  np(nx,ny,nz), &       !particle ion density at level n, n+1/2
                  vp(Ni_max,3), &       !particle velocity at t level n+1/2
                  vp1(Ni_max,3), &      !particle velocity at t level n
                  vplus(Ni_max,3), &    !v+ used in velocity update
                  vminus(Ni_max,3), &   !v- used in velocity update
                  up(nx,ny,nz,3), &     !particle flow at n, n+1/2
                  xp(Ni_max,3), &       !coordinates of ion particles
                  aj(nx,ny,nz,3), &     !curlB/(alpha*n)
                  nu(nx,ny,nz), &       !collision frequency
                  Ep(Ni_max,3), &       !Ion particle electric field
                  E(nx,ny,nz,3), &        !E field from electron mom. eq.
                  temp_p(nx,ny,nz), &   !temperature
                  mnp(nx,ny,nz), &      !mass density
                  beta, beta_p(Ni_max), &       !variable for particle scaling
                  mrat(Ni_max), &       !mass array for mulit-ion species
                  m_arr(Ni_max), &
                  input_p(3), &
                  input_E, input_Eb, bndry_Eflux, &
                  grav(nx,ny,nz)            !gravity term

      integer(4):: Ni_tot

      !Location (indices) of particles in the grid

      integer:: ijkp(Ni_max,3)

      !Weight variables for trilinear interpolation

      real:: wght(Ni_max,8)



end module Var_Arrays

