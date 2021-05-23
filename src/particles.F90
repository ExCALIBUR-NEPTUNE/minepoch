MODULE particles

  USE boundary
  USE partlist
  USE deltaf_loader2, ONLY: params_local, params_local_all

  IMPLICIT NONE

!Store some pieces of this to speed up the current evaluation.  
TYPE fields_eval_tmps
    REAL(num), DIMENSION(sf_min-1:sf_max+1) :: gx, gy, gz
    INTEGER :: cell_x1, cell_y1, cell_z1
 END TYPE fields_eval_tmps

! Some numerical factors needed for various particle-fields routines.  
    REAL(num) :: i_yz, i_xz, i_xy ! can't store these as particle steps may vary
    REAL(num) :: idx, idy, idz

    
    ! For now, fluid equations for a single species
    ! density is mass density, p_fluid is momentum density
    REAL(num), DIMENSION(:), ALLOCATABLE :: dens_fluid_next, forcet, dens_fluid
    REAL(num), DIMENSION(:), ALLOCATABLE :: dens_fluid0, p_fluid0
    REAL(num), DIMENSION(:), ALLOCATABLE :: p_fluid_next, p_fluid, pressuret
    REAL(num) :: dt_fluid, dx_fluid
    INTEGER :: nx_fluid
CONTAINS

  SUBROUTINE push_particles_2ndstep
    TYPE(particle_species), POINTER :: species, next_species
    INTEGER :: ispecies

    !Revert current back to half-step values.
    jx = jx - jx_d
    jy = jy - jy_d
    jz = jz - jz_d   

    next_species => species_list
    DO ispecies = 1, n_species
       species => next_species
       next_species => species%next

       IF (species%immobile) CYCLE

       IF (species%is_driftkinetic) THEN
          CALL push_particles_dk1(species)
       END IF
       
    ENDDO

    CALL current_bcs
    CALL particle_bcs

  END SUBROUTINE push_particles_2ndstep


  SUBROUTINE push_particles

    ! 2nd order accurate particle pusher using parabolic weighting
    ! on and off the grid. The calculation of J looks rather odd
    ! Since it works by solving d(rho)/dt = div(J) and doing a 1st order
    ! Estimate of rho(t+1.5*dt) rather than calculating J directly
    ! This gives exact charge conservation on the grid

    ! Contains the integer cell position of the particle in x, y, z
    INTEGER :: cell_x1, cell_x2, cell_x3
    INTEGER :: cell_y1, cell_y2, cell_y3
    INTEGER :: cell_z1, cell_z2, cell_z3

    ! Xi (space factor see page 38 in manual)
    ! The code now uses gx and hx instead of xi0 and xi1

    ! J from a given particle, can be spread over up to 3 cells in
    ! Each direction due to parabolic weighting. We allocate 4 or 5
    ! Cells because the position of the particle at t = t+1.5dt is not
    ! known until later. This part of the algorithm could probably be
    ! Improved, but at the moment, this is just a straight copy of
    ! The core of the PSC algorithm
    INTEGER, PARAMETER :: sf0 = sf_min, sf1 = sf_max
    REAL(num) :: jxh
    REAL(num), DIMENSION(sf0-1:sf1+1) :: jyh
    REAL(num), DIMENSION(sf0-1:sf1+1,sf0-1:sf1+1) :: jzh

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x, part_y, part_z
    REAL(num) :: part_ux, part_uy, part_uz
    REAL(num) :: part_q, part_mc, ipart_mc, part_weight

    ! Contains the floating point version of the cell number (never actually
    ! used)
    REAL(num) :: cell_x_r, cell_y_r, cell_z_r

    ! The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x, cell_frac_y, cell_frac_z

    ! Weighting factors as Eqn 4.77 page 25 of manual
    ! Eqn 4.77 would be written as
    ! F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
    ! Defined at the particle position
    REAL(num), DIMENSION(sf_min-1:sf_max+1) :: gx, gy, gz

    ! Defined at the particle position - 0.5 grid cell in each direction
    ! This is to deal with the grid stagger
    REAL(num), DIMENSION(sf_min-1:sf_max+1) :: hx, hy, hz

    ! Fields at particle location
    REAL(num) :: ex_part, ey_part, ez_part, bx_part, by_part, bz_part

    ! P+, P- and Tau variables from Boris1970, page27 of manual
    REAL(num) :: uxp, uxm, uyp, uym, uzp, uzm
    REAL(num) :: tau, taux, tauy, tauz, taux2, tauy2, tauz2

    ! charge to mass ratio modified by normalisation
    REAL(num) :: cmratio, ccmratio

    ! Used by J update
    INTEGER :: xmin, xmax, ymin, ymax, zmin, zmax
    REAL(num) :: wx, wy, wz

    ! Temporary variables
    REAL(num) :: idt, dto2, dtco2
    REAL(num) :: fcx, fcy, fcz, fjx, fjy, fjz
    REAL(num) :: root, dtfac, gamma, third
    REAL(num) :: delta_x, delta_y, delta_z
    REAL(num) :: xfac1, xfac2, yfac1, yfac2, zfac1, zfac2
    REAL(num) :: gz_iz, hz_iz, hygz, hyhz, hzyfac1, hzyfac2, yzfac
    INTEGER :: ispecies, ix, iy, iz, dcellx, dcelly, dcellz, cx, cy, cz
    INTEGER(i8) :: ipart
    ! Particle weighting multiplication factor
    REAL(num) :: cf2
    REAL(num), PARAMETER :: fac = (1.0_num / 24.0_num)**c_ndims
    INTEGER :: isubstep
    TYPE(particle), POINTER :: current, next
    TYPE(particle_species), POINTER :: species, next_species

    jx = 0.0_num
    jy = 0.0_num
    jz = 0.0_num

    jx_d = 0.0_num
    jy_d = 0.0_num
    jz_d = 0.0_num

    gx = 0.0_num
    gy = 0.0_num
    gz = 0.0_num

    ! Unvarying multiplication factors

    idx = 1.0_num / dx
    idy = 1.0_num / dy
    idz = 1.0_num / dz
    idt = 1.0_num / dt
    dto2 = dt / 2.0_num
    dtco2 = c * dto2
    dtfac = 0.5_num * dt * fac
    third = 1.0_num / 3.0_num

    i_yz = idy * idz * fac
    i_xz = idx * idz * fac
    i_xy = idx * idy * fac

    next_species => species_list
    DO ispecies = 1, n_species
       species => next_species
       next_species => species%next

       IF (species%immobile) CYCLE

       IF (species%is_driftkinetic) THEN
          CALL push_particles_dk0(species)
       ELSE
          DO isubstep=1,species%nsubstep
             CALL push_particles_lorentz_split(dt/species%nsubstep)
          END DO
       END IF

    ENDDO

    CALL current_bcs
    CALL particle_bcs
    
  contains

    SUBROUTINE push_particles_lorentz
       REAL(num) :: idtyz, idtxz, idtxy
       idtyz = idt * idy * idz * fac
       idtxz = idt * idx * idz * fac
       idtxy = idt * idx * idy * fac

      current => species%attached_list%head

      IF (.NOT. particles_uniformly_distributed) THEN
         part_weight = species%weight
         fcx = idtyz * part_weight
         fcy = idtxz * part_weight
         fcz = idtxy * part_weight
      ENDIF

      !DEC$ VECTOR ALWAYS
      DO ipart = 1, species%attached_list%count
         next => current%next
         IF (particles_uniformly_distributed) THEN
            part_weight = current%weight
            fcx = idtyz * part_weight
            fcy = idtxz * part_weight
            fcz = idtxy * part_weight
         ENDIF
         part_q   = current%charge
         part_mc  = c * current%mass
         ipart_mc = 1.0_num / part_mc
         cmratio  = part_q * dtfac * ipart_mc
         ccmratio = c * cmratio

         ! Copy the particle properties out for speed
         part_x  = current%part_pos(1) - x_grid_min_local
         part_y  = current%part_pos(2) - y_grid_min_local
         part_z  = current%part_pos(3) - z_grid_min_local
         part_ux = current%part_p(1) * ipart_mc
         part_uy = current%part_p(2) * ipart_mc
         part_uz = current%part_p(3) * ipart_mc

         ! Calculate v(t) from p(t)
         ! See PSC manual page (25-27)
         root = dtco2 / SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)

         ! Move particles to half timestep position to first order
         part_x = part_x + part_ux * root
         part_y = part_y + part_uy * root
         part_z = part_z + part_uz * root

         ! Grid cell position as a fraction.
         cell_x_r = part_x * idx
         cell_y_r = part_y * idy
         cell_z_r = part_z * idz
         ! Round cell position to nearest cell
         cell_x1 = FLOOR(cell_x_r + 0.5_num)
         ! Calculate fraction of cell between nearest cell boundary and particle
         cell_frac_x = REAL(cell_x1, num) - cell_x_r
         cell_x1 = cell_x1 + 1

         cell_y1 = FLOOR(cell_y_r + 0.5_num)
         cell_frac_y = REAL(cell_y1, num) - cell_y_r
         cell_y1 = cell_y1 + 1

         cell_z1 = FLOOR(cell_z_r + 0.5_num)
         cell_frac_z = REAL(cell_z1, num) - cell_z_r
         cell_z1 = cell_z1 + 1

         ! Particle weight factors as described in the manual, page25
         ! These weight grid properties onto particles
         ! Also used to weight particle properties onto grid, used later
         ! to calculate J
         ! NOTE: These weights require an additional multiplication factor!
#include "bspline3/gx.inc"

         ! Now redo shifted by half a cell due to grid stagger.
         ! Use shifted version for ex in X, ey in Y, ez in Z
         ! And in Y&Z for bx, X&Z for by, X&Y for bz
         cell_x2 = FLOOR(cell_x_r)
         cell_frac_x = REAL(cell_x2, num) - cell_x_r + 0.5_num
         cell_x2 = cell_x2 + 1

         cell_y2 = FLOOR(cell_y_r)
         cell_frac_y = REAL(cell_y2, num) - cell_y_r + 0.5_num
         cell_y2 = cell_y2 + 1

         cell_z2 = FLOOR(cell_z_r)
         cell_frac_z = REAL(cell_z2, num) - cell_z_r + 0.5_num
         cell_z2 = cell_z2 + 1

         dcellx = 0
         dcelly = 0
         dcellz = 0
         ! NOTE: These weights require an additional multiplication factor!
#include "bspline3/hx_dcell.inc"

         ! These are the electric and magnetic fields interpolated to the
         ! particle position. They have been checked and are correct.
         ! Actually checking this is messy.
#include "bspline3/e_part.inc"
#include "bspline3/b_part.inc"

         ! update particle momenta using weighted fields
         uxm = part_ux + cmratio * ex_part
         uym = part_uy + cmratio * ey_part
         uzm = part_uz + cmratio * ez_part

         ! Half timestep, then use Boris1970 rotation, see Birdsall and Langdon
         root = ccmratio / SQRT(uxm**2 + uym**2 + uzm**2 + 1.0_num)

         taux = bx_part * root
         tauy = by_part * root
         tauz = bz_part * root

         taux2 = taux**2
         tauy2 = tauy**2
         tauz2 = tauz**2

         tau = 1.0_num / (1.0_num + taux2 + tauy2 + tauz2)

         uxp = ((1.0_num + taux2 - tauy2 - tauz2) * uxm &
              + 2.0_num * ((taux * tauy + tauz) * uym &
              + (taux * tauz - tauy) * uzm)) * tau
         uyp = ((1.0_num - taux2 + tauy2 - tauz2) * uym &
              + 2.0_num * ((tauy * tauz + taux) * uzm &
              + (tauy * taux - tauz) * uxm)) * tau
         uzp = ((1.0_num - taux2 - tauy2 + tauz2) * uzm &
              + 2.0_num * ((tauz * taux + tauy) * uxm &
              + (tauz * tauy - taux) * uym)) * tau

         ! Rotation over, go to full timestep
         part_ux = uxp + cmratio * ex_part
         part_uy = uyp + cmratio * ey_part
         part_uz = uzp + cmratio * ez_part

         ! Calculate particle velocity from particle momentum
         gamma = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)
         root = dtco2 / gamma

         delta_x = part_ux * root
         delta_y = part_uy * root
         delta_z = part_uz * root

         ! Move particles to end of time step at 2nd order accuracy
         part_x = part_x + delta_x
         part_y = part_y + delta_y
         part_z = part_z + delta_z

         ! particle has now finished move to end of timestep, so copy back
         ! into particle array
         current%part_pos = (/ part_x + x_grid_min_local, &
              part_y + y_grid_min_local, part_z + z_grid_min_local /)
         current%part_p   = part_mc * (/ part_ux, part_uy, part_uz /)

         ! Original code calculates densities of electrons, ions and neutrals
         ! here. This has been removed to reduce memory footprint

         ! Now advance to t+1.5dt to calculate current. This is detailed in
         ! the manual between pages 37 and 41. The version coded up looks
         ! completely different to that in the manual, but is equivalent.
         ! Use t+1.5 dt so that can update J to t+dt at 2nd order
         part_x = part_x + delta_x
         part_y = part_y + delta_y
         part_z = part_z + delta_z

         cell_x_r = part_x * idx
         cell_y_r = part_y * idy
         cell_z_r = part_z * idz

         cell_x3 = FLOOR(cell_x_r + 0.5_num)
         cell_frac_x = REAL(cell_x3, num) - cell_x_r
         cell_x3 = cell_x3 + 1

         cell_y3 = FLOOR(cell_y_r + 0.5_num)
         cell_frac_y = REAL(cell_y3, num) - cell_y_r
         cell_y3 = cell_y3 + 1

         cell_z3 = FLOOR(cell_z_r + 0.5_num)
         cell_frac_z = REAL(cell_z3, num) - cell_z_r
         cell_z3 = cell_z3 + 1

         hx = 0.0_num
         hy = 0.0_num
         hz = 0.0_num

         dcellx = cell_x3 - cell_x1
         dcelly = cell_y3 - cell_y1
         dcellz = cell_z3 - cell_z1
         ! NOTE: These weights require an additional multiplication factor!
#include "bspline3/hx_dcell.inc"

         ! Now change Xi1* to be Xi1*-Xi0*. This makes the representation of
         ! the current update much simpler
         hx = hx - gx
         hy = hy - gy
         hz = hz - gz

         ! Remember that due to CFL condition particle can never cross more
         ! than one gridcell in one timestep

         xmin = sf_min + (dcellx - 1) / 2
         xmax = sf_max + (dcellx + 1) / 2

         ymin = sf_min + (dcelly - 1) / 2
         ymax = sf_max + (dcelly + 1) / 2

         zmin = sf_min + (dcellz - 1) / 2
         zmax = sf_max + (dcellz + 1) / 2

         fjx = fcx * part_q
         fjy = fcy * part_q
         fjz = fcz * part_q

         jzh = 0.0_num
         DO iz = zmin, zmax
            cz = cell_z1 + iz
            zfac1 =         gz(iz) + 0.5_num * hz(iz)
            zfac2 = third * hz(iz) + 0.5_num * gz(iz)

            gz_iz = gz(iz)
            hz_iz = hz(iz)

            jyh = 0.0_num
            DO iy = ymin, ymax
               cy = cell_y1 + iy
               yfac1 =         gy(iy) + 0.5_num * hy(iy)
               yfac2 = third * hy(iy) + 0.5_num * gy(iy)

               hygz = hy(iy) * gz_iz
               hyhz = hy(iy) * hz_iz
               yzfac = gy(iy) * zfac1 + hy(iy) * zfac2
               hzyfac1 = hz_iz * yfac1
               hzyfac2 = hz_iz * yfac2

               jxh = 0.0_num
               DO ix = xmin, xmax
                  cx = cell_x1 + ix
                  xfac1 =         gx(ix) + 0.5_num * hx(ix)
                  xfac2 = third * hx(ix) + 0.5_num * gx(ix)

                  wx = hx(ix) * yzfac
                  wy = xfac1 * hygz + xfac2 * hyhz
                  wz = gx(ix) * hzyfac1 + hx(ix) * hzyfac2

                  ! This is the bit that actually solves d(rho)/dt = -div(J)
                  jxh = jxh - fjx * wx
                  jyh(ix) = jyh(ix) - fjy * wy
                  jzh(ix, iy) = jzh(ix, iy) - fjz * wz

                  jx(cx, cy, cz) = jx(cx, cy, cz) + jxh
                  jy(cx, cy, cz) = jy(cx, cy, cz) + jyh(ix)
                  jz(cx, cy, cz) = jz(cx, cy, cz) + jzh(ix, iy)
               ENDDO
            ENDDO
         ENDDO

         current => next
      ENDDO

    END SUBROUTINE push_particles_lorentz

    SUBROUTINE push_particles_lorentz_split(dt_sub)
      REAL(num), intent(IN) :: dt_sub 
      TYPE(fields_eval_tmps) :: st_half
      REAL(num), DIMENSION(3) :: part_pos_t1p5, pos_half, Bvec, Evec
      REAL(num) :: weight_back, part_qfac
      REAL(num), DIMENSION(3) :: force, part_v
      REAL(num) :: idt, dto2, dtco2, idt0
      REAL(num) :: dtfac
      
      idt = 1.0_num / dt_sub
      idt0= 1.0_num / dt
      dto2 = dt_sub / 2.0_num
      dtco2 = c * dto2
      dtfac = 0.5_num * dt_sub * fac

      IF (species%solve_fluid) CALL initstep_fluid
      
      current => species%attached_list%head

      IF (.NOT. particles_uniformly_distributed) THEN
         part_weight = species%weight 
      ENDIF

      !DEC$ VECTOR ALWAYS
      DO ipart = 1, species%attached_list%count
         next => current%next
         IF (particles_uniformly_distributed) THEN
            part_weight = current%weight
         ENDIF

         part_q    = current%charge
         part_qfac = part_q * idt0
         part_mc  = c * current%mass
         ipart_mc = 1.0_num / part_mc
         cmratio  = part_q * dtfac * ipart_mc
         ccmratio = c * cmratio

         ! Copy the particle properties out for speed
         part_ux = current%part_p(1) * ipart_mc
         part_uy = current%part_p(2) * ipart_mc
         part_uz = current%part_p(3) * ipart_mc

         ! Calculate v(t) from p(t)
         ! See PSC manual page (25-27)
         root = dtco2 / SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)

         ! Move particles to half timestep position to first order
         pos_half = current%part_pos + root * (/ part_ux, part_uy, part_uz /)
         CALL get_fields_at_point_store(pos_half,Bvec,Evec,st_half)

         ! update particle momenta using weighted fields
         uxm = part_ux + cmratio * Evec(1)
         uym = part_uy + cmratio * Evec(2)
         uzm = part_uz + cmratio * Evec(3)

         ! Half timestep, then use Boris1970 rotation, see Birdsall and Langdon
         root = ccmratio / SQRT(uxm**2 + uym**2 + uzm**2 + 1.0_num)

         taux = Bvec(1) * root
         tauy = Bvec(2) * root
         tauz = Bvec(3) * root

         taux2 = taux**2
         tauy2 = tauy**2
         tauz2 = tauz**2

         tau = 1.0_num / (1.0_num + taux2 + tauy2 + tauz2)

         uxp = ((1.0_num + taux2 - tauy2 - tauz2) * uxm &
              + 2.0_num * ((taux * tauy + tauz) * uym &
              + (taux * tauz - tauy) * uzm)) * tau
         uyp = ((1.0_num - taux2 + tauy2 - tauz2) * uym &
              + 2.0_num * ((tauy * tauz + taux) * uzm &
              + (tauy * taux - tauz) * uxm)) * tau
         uzp = ((1.0_num - taux2 - tauy2 + tauz2) * uzm &
              + 2.0_num * ((tauz * taux + tauy) * uxm &
              + (tauz * tauy - taux) * uym)) * tau

         ! Rotation over, go to full timestep
         part_ux = uxp + cmratio * Evec(1)
         part_uy = uyp + cmratio * Evec(2)
         part_uz = uzp + cmratio * Evec(3)

         ! Calculate particle velocity from particle momentum
         gamma = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)
         root = dtco2 / gamma

         delta_x = part_ux * root
         delta_y = part_uy * root
         delta_z = part_uz * root

         ! Move particles to end of time step at 2nd order accuracy
         ! particle has now finished move to end of timestep, place data
         ! in particle array
         current%part_pos = pos_half + (/delta_x, delta_y, delta_z/)
         current%part_p   = part_mc * (/ part_ux, part_uy, part_uz /)

         ! Delta-f calcuation: subtract background from
         ! calculated current.
         IF(species%use_deltaf) THEN
            weight_back = current%pvol * f0(species, current, current%mass)
         END IF
            
         ! Now advance to t+1.5dt to calculate current. This is detailed in
         ! the manual between pages 37 and 41. The version coded up looks
         ! completely different to that in the manual, but is equivalent.
         ! Use t+1.5 dt so that can update J to t+dt at 2nd order
         part_pos_t1p5 = current%part_pos + (/ delta_x, delta_y, delta_z /)
         !Current deposition uses position at t+0.5dt and t+1.5dt, particle
         !assumed to travel in direct line between these locations. Second order 
         !in time for evaluation of current at t+dt 
         CALL current_deposition_store(st_half,part_pos_t1p5,(part_weight*part_qfac),.false.)

         if (species%solve_fluid) then
            part_v = (/part_ux, part_uy,part_uz/)
            force = part_q*(Evec+cross(part_v,Bvec))
            CALL assignweight_fluid(current%part_pos, &
                 & (part_weight-weight_back)*current%part_p(1), &
                 & (part_weight-weight_back)*current%mass,force )
         end if
         
         current => next
      ENDDO

      IF (species%solve_fluid) then
         CALL fluideq_diag
         CALL update_fluideq
      end if
      
    END SUBROUTINE push_particles_lorentz_split

  END SUBROUTINE push_particles

    !Drift-kinetic push. Because we can't do leapfrog, do an explicit two-step scheme.
    ! -first substep is just Euler.
    ! -particles stored at half-timestep to make current update easier.
  SUBROUTINE push_particles_dk0(species)
      TYPE(fields_eval_tmps) :: st_0
      TYPE(particle_species), INTENT(INOUT) :: species
      TYPE(particle), POINTER :: current, next
      REAL(num), DIMENSION(3) :: dRdt, pos_0, pos_h
      REAL(num), DIMENSION(3) :: bdir, drifts_mu, drifts_vpll
      REAL(num), DIMENSION(3) :: drifts_ExB
      REAL(num), DIMENSION(3) :: Evec, Bvec
      REAL(num), DIMENSION(3,3):: Btens

      REAL(num) :: part_mu, part_u, bdotBmag
      REAL(num) :: part_u_0,part_u_h, dudt
      REAL(num) :: part_q, part_weight, part_qfac
      INTEGER(i8) :: ipart
      REAL(num) :: idt
      idt = 1.0_num/dt


      current => species%attached_list%head

      IF (.NOT. particles_uniformly_distributed) THEN
         part_weight = species%weight
      ENDIF

      !DEC$ VECTOR ALWAYS
      DO ipart = 1, species%attached_list%count
         next => current%next
         IF (particles_uniformly_distributed) THEN
            part_weight = current%weight
         ENDIF

         part_q   = current%charge
         part_qfac= part_q * idt
         ! Do nonrelativistic drift-kinetics for the moment.
         part_u   = current%part_p(1)
         part_mu  = current%part_p(2)

         IF (current%part_p(3)<0.5_num) THEN
            write (18,*) ipart,current%part_pos(1),current%part_pos(2),current%part_pos(3),part_u,part_mu,current%part_p(3)
         END IF

         CALL get_fields_at_point_store(current%part_pos,Bvec,Evec,st_0,Btens)
         CALL get_drifts(current%part_pos,Evec,Bvec,Btens,drifts_ExB,bdir, &
              &  drifts_mu,drifts_vpll, bdotBmag)

         dRdt = bdir * part_u
         dudt = - (part_mu * bdotBmag / current%mass) 
         ! Move particles to half timestep position to first order
         pos_0 = current%part_pos + dRdt * dt
         part_u_0 = part_u + dudt * dt

         !Half step position.
         pos_h = 0.5*(current%part_pos + pos_0)
         part_u_h = 0.5*(part_u + part_u_0)
         current%work(1:3)  = pos_h
         current%work(4)    = part_u_h

         !Do current deposition using lowest order current. (current at t_{N+1})
         !Before we apply this current, E+B need to be stored in a temporary;
         !this is the lowest order current.
         CALL current_deposition_store(st_0,pos_0,(part_weight*part_qfac),.true.)

         
         current => next
      ENDDO
      
    END SUBROUTINE push_particles_dk0

    !Drift-kinetic push. Because we can't do leapfrog, do an explicit two-step scheme.
    !second substep: evaluate derivatives at t+dt using next-step fields.
    SUBROUTINE push_particles_dk1(species)
      TYPE(fields_eval_tmps) :: st_half
      TYPE(particle_species), INTENT(INOUT) :: species
      TYPE(particle), POINTER :: current, next
      REAL(num), DIMENSION(3) :: pos_0, pos_h, dRdt_1
      REAL(num), DIMENSION(3) :: bdir, drifts_mu, drifts_vpll
      REAL(num), DIMENSION(3) :: drifts_ExB
      REAL(num), DIMENSION(3) :: Evec, Bvec
      REAL(num), DIMENSION(3,3):: Btens

      REAL(num) :: part_mu, part_u, bdotBmag
      REAL(num) :: part_u_h, dudt_1
      REAL(num) :: part_q, part_weight, part_qfac
      INTEGER(i8) :: ipart
      REAL(num) :: idt
      idt = 1.0_num/dt
      
      current => species%attached_list%head

      IF (.NOT. particles_uniformly_distributed) THEN
         part_weight = species%weight
      ENDIF

      !DEC$ VECTOR ALWAYS
      DO ipart = 1, species%attached_list%count
         next => current%next
         IF (particles_uniformly_distributed) THEN
            part_weight = current%weight
         ENDIF

         part_q   = current%charge
         part_qfac= part_q * idt
         ! Do nonrelativistic drift-kinetics for the moment.
         part_u   = current%part_p(1)
         part_mu  = current%part_p(2)
 
         !Half step position.
         pos_h = current%work(1:3)
         part_u_h = current%work(4)

         pos_0 = current%part_pos
         ! From this point on: need to advance fields half a step in time: we could do this 
         ! using the PIC particle currents + estimated drift currents:

         CALL get_fields_at_point_store(pos_h,Bvec,Evec,st_half,Btens)
         CALL get_drifts(pos_h,Evec,Bvec,Btens,drifts_ExB,bdir, &
              &  drifts_mu,drifts_vpll, bdotBmag)
         dRdt_1 = bdir * part_u_h
         dudt_1 = - (part_mu * bdotBmag / current%mass) 
         
         ! particle has now finished move to end of timestep, so copy back
         ! into particle array
         part_u = part_u + dudt_1 * dt
         current%part_pos = current%part_pos + dRdt_1 * dt 
         current%part_p(1:2)   = (/ part_u, part_mu /)

         !This is the current between step N+1/2 and N+3/2 so 2nd order at N+1
         CALL current_deposition(pos_0,current%part_pos,(part_weight*part_q),.true.)

         current => next
      ENDDO

    END SUBROUTINE push_particles_dk1


  SUBROUTINE get_drifts(pos,Evec,Bvec,Btens,ExB,bdir,drifts_mu,drifts_vpll,bdir_dotgradBmag)    
    REAL(num), DIMENSION(3), INTENT(INOUT) :: pos, ExB, Evec, Bvec
    REAL(num), DIMENSION(3,3), INTENT(INOUT) :: Btens
    REAL(num), DIMENSION(3), INTENT(OUT)   :: & 
         & drifts_mu,drifts_vpll
    REAL(num), INTENT(OUT)   :: bdir_dotgradBmag 

    REAL(num), DIMENSION(3) :: Bdir,BdotgradB,gradBmag
    REAL(num) :: Bsq,Bnorm

    Bsq = dot_product(Bvec,Bvec)
    Bnorm = sqrt(Bsq)
    bdir = Bvec/Bnorm

    !Maybe can do these with matmul?
    gradBmag = (Btens(:,1)*Bvec(1) + Btens(:,2)*Bvec(2) + Btens(:,3)*Bvec(3))/Bnorm 
    BdotgradB = Bvec(1)*Btens(1,:) + Bvec(2)*Btens(2,:) + Bvec(3)*Btens(3,:)

    bdir_dotgradBmag = dot_product(bdir,gradBmag)

    ExB = cross(Evec,Bvec)/Bsq

    drifts_vpll = cross(BdotgradB,Bvec)/(Bnorm*Bsq)
    drifts_mu   = cross(gradBmag,Bvec)/Bsq
  END SUBROUTINE get_drifts


  FUNCTION cross(a,b)
    REAL(num), DIMENSION(3) :: cross
    REAL(num), DIMENSION(3) :: a,b
    cross(1) = a(2)*b(3)-a(3)*b(2)
    cross(2) = a(3)*b(1)-a(1)*b(3)
    cross(3) = a(1)*b(2)-a(2)*b(1)
  END FUNCTION cross

  !Find derivative of weight array with respect to position.
  !This is stored in array(:,2)
  !array(:,1) is the weight array itself.
  !cell_frac is negative position in cell relative to centre?
  SUBROUTINE h_derivs(array,array1,offset,cell_frac,id_cell)
    REAL(num), DIMENSION(sf_min-1:sf_max+1,2), INTENT(INOUT) :: array
    REAL(num), DIMENSION(sf_min-1:sf_max+1), INTENT(INOUT) :: array1
    REAL(num), INTENT(IN) :: cell_frac,id_cell
    INTEGER, INTENT(IN)   :: offset
    REAL(num) :: cf2

    cf2 = cell_frac*cell_frac
    array = 0.0_num
    array(offset-2,2) =  4.0_num*(0.5_num + cell_frac)**3
    array(offset-1,2) =  11.0_num &
         + 4.0_num * cell_frac * (3.0_num - 3.0_num*cell_frac - 4.0_num*cf2) 
    array(offset  ,2) =     6.0_num * (4.0_num*cell_frac) * (cf2 - 1.25_num) 
    array(offset+1,2) = -11.0_num &
         + 4.0_num * cell_frac * (3.0_num + 3.0_num*cell_frac - 4.0_num*cf2) 
    array(offset+2,2) = -4.0_num*(0.5_num - cell_frac)**3

    !Because cell_frac is negatively proportional to position.
    array(:,2) =-array(:,2)*id_cell
    array(:,1) = array1
  END SUBROUTINE h_derivs

  !Rewrite finite-difference evaluation include file as a function for
  !testing purposes. 
  SUBROUTINE hfun(array,offset,cell_frac)
    REAL(num), DIMENSION(sf_min-1:sf_max+1), INTENT(INOUT) :: array
    REAL(num), INTENT(IN) :: cell_frac
    INTEGER, INTENT(IN)   :: offset
    REAL(num) :: cf2

    cf2 = cell_frac**2
    array(offset-2) = (0.5_num + cell_frac)**4
    array(offset-1) = 4.75_num + 11.0_num * cell_frac &
         + 4.0_num * cf2 * (1.5_num - cell_frac - cf2)
    array(offset+0) = 14.375_num + 6.0_num * cf2 * (cf2 - 2.5_num)
    array(offset+1) = 4.75_num - 11.0_num * cell_frac &
         + 4.0_num * cf2 * (1.5_num + cell_frac - cf2)
    array(offset+2) = (0.5_num - cell_frac)**4
  END SUBROUTINE hfun


  !Need this for general grad B drift, calculating B*, etc.
  !For the moment, aiming for code that is easy to modify/test rather
  !than max performance.
  SUBROUTINE calc_Btens(Btens,hdx,hdy,hdz,gdx,gdy,gdz,idx,idy,idz, &
       & cell_x1,cell_x2,cell_y1,cell_y2,cell_z1,cell_z2)
    INTEGER :: cell_x1,cell_x2,cell_y1,cell_y2,cell_z1,cell_z2
    REAL(num), DIMENSION(3,3) :: Btens
    REAL(num), DIMENSION(sf_min-1:sf_max+1,2) :: gdx, gdy, gdz
    REAL(num), DIMENSION(sf_min-1:sf_max+1,2) :: hdx, hdy, hdz
    REAL(num) :: idx,idy,idz
    REAL(num), PARAMETER :: fac = (1.0_num / 24.0_num)**c_ndims

    INTEGER :: cello_x,cello_y,cello_z
    INTEGER :: xp,yp,zp !Keep track of which derivative is being taken.
    INTEGER :: ii

    Btens = 0.0_num
    ! Calculate grad-B tensor
    do cello_x = -2,2
       do cello_y = -2,2
          do cello_z = -2,2
             do ii=1,3
                xp=1+kronecker_delta(ii,1)
                yp=1+kronecker_delta(ii,2)
                zp=1+kronecker_delta(ii,3)
                Btens(ii,1) = Btens(ii,1) + &
                     &  gdx(cello_x,xp)*hdy(cello_y,yp)*hdz(cello_z,zp)  &
                     & *bx(cell_x1+cello_x,cell_y2+cello_y,cell_z2+cello_z)
                Btens(ii,2) = Btens(ii,2) + &
                     &  hdx(cello_x,xp)*gdy(cello_y,yp)*hdz(cello_z,zp)  &
                     & *by(cell_x2+cello_x,cell_y1+cello_y,cell_z2+cello_z)
                Btens(ii,3) = Btens(ii,3) + &
                     &  hdx(cello_x,xp)*hdy(cello_y,yp)*gdz(cello_z,zp)  &
                     & *bz(cell_x2+cello_x,cell_y2+cello_y,cell_z1+cello_z)
             end do
          end do
       end do
    end do
  END SUBROUTINE calc_Btens

  INTEGER FUNCTION kronecker_delta(a,b)
    INTEGER, INTENT(IN) :: a,b
    IF (a==b) THEN
       kronecker_delta=1
    ELSE
       kronecker_delta=0
    END IF
  END FUNCTION KRONECKER_DELTA

  !Do current deposition of particle moving along straight line
  !from pos0 to pos1, with chargeweight = weight*charge
  SUBROUTINE current_deposition(pos0,pos1,chargeweight,drift_switch)
    REAL(num), DIMENSION(3), INTENT(INOUT) :: pos0,pos1
    REAL(num), INTENT(IN) :: chargeweight
    LOGICAL, INTENT(IN) :: drift_switch
    
    TYPE(fields_eval_tmps)  :: st0

    CALL calc_stdata(pos0,st0)
    CALL current_deposition_store(st0,pos1,chargeweight,drift_switch)
  END SUBROUTINE current_deposition

  !Do current deposition, reusing some precalculated data (in st)
  !Particle has moved from pos0 (data stored in st) to pos
  SUBROUTINE current_deposition_store(st,pos,chargeweight,drift_switch)
    REAL(num), DIMENSION(3), INTENT(INOUT) :: pos
    TYPE(fields_eval_tmps), INTENT(INOUT)  :: st
    REAL(num), INTENT(IN) :: chargeweight
    LOGICAL, INTENT(IN) :: drift_switch

    REAL(num) :: cell_x_r, cell_y_r, cell_z_r
    REAL(num) :: part_x, part_y, part_z
    INTEGER :: cell_x3
    INTEGER :: cell_y3
    INTEGER :: cell_z3
    INTEGER :: dcellx, dcelly, dcellz
    REAL(num) :: cf2
    ! Defined at the particle position - 0.5 grid cell in each direction
    ! This is to deal with the grid stagger
    REAL(num), DIMENSION(sf_min-1:sf_max+1) :: hx, hy, hz
    ! J from a given particle, can be spread over up to 3 cells in
    ! Each direction due to parabolic weighting. We allocate 4 or 5
    ! Cells because the position of the particle at t = t+1.5dt is not
    ! known until later. This part of the algorithm could probably be
    ! Improved, but at the moment, this is just a straight copy of
    ! The core of the PSC algorithm
    INTEGER, PARAMETER :: sf0 = sf_min, sf1 = sf_max
    REAL(num) :: jxh
    REAL(num), DIMENSION(sf0-1:sf1+1) :: jyh
    REAL(num), DIMENSION(sf0-1:sf1+1,sf0-1:sf1+1) :: jzh
    REAL(num) :: fjx, fjy, fjz
    ! The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x, cell_frac_y, cell_frac_z
    INTEGER :: ix, iy, iz, cx, cy, cz
    REAL(num) :: gz_iz, hz_iz, hygz, hyhz, hzyfac1, hzyfac2, yzfac
    ! Used by J update
    INTEGER :: xmin, xmax, ymin, ymax, zmin, zmax
    REAL(num) :: wx, wy, wz
    REAL(num), PARAMETER :: third=(1.0_num / 3.0_num)
    REAL(num) :: xfac1, xfac2, yfac1, yfac2, zfac1, zfac2
  

    part_x = pos(1) - x_grid_min_local
    part_y = pos(2) - y_grid_min_local
    part_z = pos(3) - z_grid_min_local

    cell_x_r = part_x * idx
    cell_y_r = part_y * idy
    cell_z_r = part_z * idz

    cell_x3 = FLOOR(cell_x_r + 0.5_num)
    cell_frac_x = REAL(cell_x3, num) - cell_x_r
    cell_x3 = cell_x3 + 1

    cell_y3 = FLOOR(cell_y_r + 0.5_num)
    cell_frac_y = REAL(cell_y3, num) - cell_y_r
    cell_y3 = cell_y3 + 1
    
    cell_z3 = FLOOR(cell_z_r + 0.5_num)
    cell_frac_z = REAL(cell_z3, num) - cell_z_r
    cell_z3 = cell_z3 + 1

    hx = 0.0_num
    hy = 0.0_num
    hz = 0.0_num

    dcellx = cell_x3 - st%cell_x1
    dcelly = cell_y3 - st%cell_y1
    dcellz = cell_z3 - st%cell_z1
    ! NOTE: These weights require an additional multiplication factor!
#include "bspline3/hx_dcell.inc"

    ! Now change Xi1* to be Xi1*-Xi0*. This makes the representation of
    ! the current update much simpler
    hx = hx - st%gx
    hy = hy - st%gy
    hz = hz - st%gz

    ! Remember that due to CFL condition particle can never cross more
    ! than one gridcell in one timestep

    xmin = sf_min + (dcellx - 1) / 2
    xmax = sf_max + (dcellx + 1) / 2
    
    ymin = sf_min + (dcelly - 1) / 2
    ymax = sf_max + (dcelly + 1) / 2
    
    zmin = sf_min + (dcellz - 1) / 2
    zmax = sf_max + (dcellz + 1) / 2

    fjx = i_yz * chargeweight
    fjy = i_xz * chargeweight
    fjz = i_xy * chargeweight
    
    jzh = 0.0_num
    DO iz = zmin, zmax
       cz = st%cell_z1 + iz
       zfac1 =         st%gz(iz) + 0.5_num * hz(iz)
       zfac2 = third * hz(iz) + 0.5_num * st%gz(iz)
       
       gz_iz = st%gz(iz)
       hz_iz = hz(iz)
       
       jyh = 0.0_num
       DO iy = ymin, ymax
          cy = st%cell_y1 + iy
          yfac1 =         st%gy(iy) + 0.5_num * hy(iy)
          yfac2 = third * hy(iy) + 0.5_num * st%gy(iy)
          
          hygz = hy(iy) * gz_iz
          hyhz = hy(iy) * hz_iz
          yzfac = st%gy(iy) * zfac1 + hy(iy) * zfac2
          hzyfac1 = hz_iz * yfac1
          hzyfac2 = hz_iz * yfac2
          
          jxh = 0.0_num
          DO ix = xmin, xmax
             cx = st%cell_x1 + ix
             xfac1 =         st%gx(ix) + 0.5_num * hx(ix)
             xfac2 = third * hx(ix) + 0.5_num * st%gx(ix)
             
             wx = hx(ix) * yzfac
             wy = xfac1 * hygz + xfac2 * hyhz
             wz = st%gx(ix) * hzyfac1 + hx(ix) * hzyfac2
             
             ! This is the bit that actually solves d(rho)/dt = -div(J)
             jxh = jxh - fjx * wx
             jyh(ix) = jyh(ix) - fjy * wy
             jzh(ix, iy) = jzh(ix, iy) - fjz * wz
             if (.NOT.drift_switch) then
                jx(cx, cy, cz) = jx(cx, cy, cz) + jxh
                jy(cx, cy, cz) = jy(cx, cy, cz) + jyh(ix)
                jz(cx, cy, cz) = jz(cx, cy, cz) + jzh(ix, iy)
             else 
                jx_d(cx, cy, cz) = jx_d(cx, cy, cz) + jxh
                jy_d(cx, cy, cz) = jy_d(cx, cy, cz) + jyh(ix)
                jz_d(cx, cy, cz) = jz_d(cx, cy, cz) + jzh(ix, iy)
             end if
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE current_deposition_store

  SUBROUTINE get_fields_at_point(pos,bvec,evec,btens)
    REAL(num), DIMENSION(3),   INTENT(INOUT) :: pos,bvec,evec
    REAL(num), DIMENSION(3,3), INTENT(INOUT) :: btens
    TYPE(fields_eval_tmps) :: st
    CALL get_fields_at_point_store(pos,bvec,evec,st,btens)
  END SUBROUTINE get_fields_at_point
 
  ! Utility routine for calculating cell offsets and weights
  ! at a position.
  SUBROUTINE calc_stdata(pos,st)
    TYPE(fields_eval_tmps), INTENT(INOUT) :: st
    REAL(num), DIMENSION(3),   INTENT(INOUT) :: pos
    INTEGER :: cell_x1, cell_y1, cell_z1
    REAL(num) :: cell_x_r, cell_y_r, cell_z_r
    REAL(num) :: part_x, part_y, part_z
    REAL(num), DIMENSION(sf_min-1:sf_max+1) :: gx, gy, gz
    REAL(num) :: cell_frac_x, cell_frac_y, cell_frac_z
    REAL(num) :: cf2

    part_x = pos(1) - x_grid_min_local
    part_y = pos(2) - y_grid_min_local
    part_z = pos(3) - z_grid_min_local

    ! Grid cell position as a fraction.
    cell_x_r = part_x * idx
    cell_y_r = part_y * idy
    cell_z_r = part_z * idz
    ! Round cell position to nearest cell
    cell_x1 = FLOOR(cell_x_r + 0.5_num)
    ! Calculate fraction of cell between nearest cell boundary and particle
    cell_frac_x = REAL(cell_x1, num) - cell_x_r
    cell_x1 = cell_x1 + 1

    cell_y1 = FLOOR(cell_y_r + 0.5_num)
    cell_frac_y = REAL(cell_y1, num) - cell_y_r
    cell_y1 = cell_y1 + 1

    cell_z1 = FLOOR(cell_z_r + 0.5_num)
    cell_frac_z = REAL(cell_z1, num) - cell_z_r
    cell_z1 = cell_z1 + 1

    ! Particle weight factors as described in the manual, page25
    ! These weight grid properties onto particles
    ! Also used to weight particle properties onto grid, used later
    ! to calculate J
    ! NOTE: These weights require an additional multiplication factor!
    gx=0
    gy=0
    gz=0
#include "bspline3/gx.inc"

    !Temporary storage for current deposition.
    st%gx = gx
    st%gy = gy
    st%gz = gz

    st%cell_x1 = cell_x1
    st%cell_y1 = cell_y1
    st%cell_z1 = cell_z1    
  END SUBROUTINE calc_stdata


  ! Evaluate fields at a point.
  SUBROUTINE get_fields_at_point_store(pos,bvec,evec,st,btens)
    TYPE(fields_eval_tmps), INTENT(INOUT) :: st
    REAL(num), DIMENSION(3),   INTENT(INOUT) :: pos,bvec,evec
    REAL(num), DIMENSION(3,3), INTENT(INOUT), OPTIONAL :: btens
    REAL(num), PARAMETER :: fac = (1.0_num / 24.0_num)**c_ndims

    ! Fields at particle location
    REAL(num) :: ex_part, ey_part, ez_part, bx_part, by_part, bz_part
    REAL(num) :: cell_x_r, cell_y_r, cell_z_r
    INTEGER :: cell_x1, cell_x2
    INTEGER :: cell_y1, cell_y2
    INTEGER :: cell_z1, cell_z2
    INTEGER :: dcellx, dcelly, dcellz
    REAL(num) :: cf2
    REAL(num) :: part_x, part_y, part_z
    ! The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x, cell_frac_y, cell_frac_z

    ! Weighting factors as Eqn 4.77 page 25 of manual
    ! Eqn 4.77 would be written as
    ! F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
    ! Defined at the particle position
    REAL(num), DIMENSION(sf_min-1:sf_max+1) :: gx, gy, gz
    ! Spatial derivative of same
    REAL(num), DIMENSION(sf_min-1:sf_max+1,2) :: gdx, gdy, gdz
    ! Defined at the particle position - 0.5 grid cell in each direction
    ! This is to deal with the grid stagger
    REAL(num), DIMENSION(sf_min-1:sf_max+1) :: hx, hy, hz
    ! Spatial derivative of same
    REAL(num), DIMENSION(sf_min-1:sf_max+1,2) :: hdx, hdy, hdz

    gx=0
    gy=0
    gz=0

    part_x = pos(1) - x_grid_min_local
    part_y = pos(2) - y_grid_min_local
    part_z = pos(3) - z_grid_min_local

    ! Grid cell position as a fraction.
    cell_x_r = part_x * idx
    cell_y_r = part_y * idy
    cell_z_r = part_z * idz
    ! Round cell position to nearest cell
    cell_x1 = FLOOR(cell_x_r + 0.5_num)
    ! Calculate fraction of cell between nearest cell boundary and particle
    cell_frac_x = REAL(cell_x1, num) - cell_x_r
    cell_x1 = cell_x1 + 1

    cell_y1 = FLOOR(cell_y_r + 0.5_num)
    cell_frac_y = REAL(cell_y1, num) - cell_y_r
    cell_y1 = cell_y1 + 1

    cell_z1 = FLOOR(cell_z_r + 0.5_num)
    cell_frac_z = REAL(cell_z1, num) - cell_z_r
    cell_z1 = cell_z1 + 1

    ! Particle weight factors as described in the manual, page25
    ! These weight grid properties onto particles
    ! Also used to weight particle properties onto grid, used later
    ! to calculate J
    ! NOTE: These weights require an additional multiplication factor!
#include "bspline3/gx.inc"
    IF(present(Btens)) THEN
       CALL h_derivs(gdx,gx,0,cell_frac_x,idx)
       CALL h_derivs(gdy,gy,0,cell_frac_y,idy)
       CALL h_derivs(gdz,gz,0,cell_frac_z,idz)
    END IF
    ! Now redo shifted by half a cell due to grid stagger.
    ! Use shifted version for ex in X, ey in Y, ez in Z
    ! And in Y&Z for bx, X&Z for by, X&Y for bz
    cell_x2 = FLOOR(cell_x_r)
    cell_frac_x = REAL(cell_x2, num) - cell_x_r + 0.5_num
    cell_x2 = cell_x2 + 1

    cell_y2 = FLOOR(cell_y_r)
    cell_frac_y = REAL(cell_y2, num) - cell_y_r + 0.5_num
    cell_y2 = cell_y2 + 1

    cell_z2 = FLOOR(cell_z_r)
    cell_frac_z = REAL(cell_z2, num) - cell_z_r + 0.5_num
    cell_z2 = cell_z2 + 1

    dcellx = 0
    dcelly = 0
    dcellz = 0
    ! NOTE: These weights require an additional multiplication factor!
#include "bspline3/hx_dcell.inc"
    ! These are the electric and magnetic fields interpolated to the
    ! particle position. They have been checked and are correct.
    ! Actually checking this is messy.
#include "bspline3/e_part.inc"
#include "bspline3/b_part.inc"
    Evec(1) = ex_part
    Evec(2) = ey_part
    Evec(3) = ez_part
    Bvec(1) = Bx_part
    Bvec(2) = By_part
    Bvec(3) = Bz_part

     !Temporary storage for current deposition.
     st%gx = gx
     st%gy = gy
     st%gz = gz

     st%cell_x1 = cell_x1
     st%cell_y1 = cell_y1
     st%cell_z1 = cell_z1    
     
     IF(present(Btens)) THEN
       CALL h_derivs(hdx,hx,dcellx,cell_frac_x,idx)
       CALL h_derivs(hdy,hy,dcelly,cell_frac_y,idy)
       CALL h_derivs(hdz,hz,dcellz,cell_frac_z,idz)
       CALL calc_Btens(Btens,hdx,hdy,hdz,gdx,gdy,gdz,idx,idy,idz, &
            & cell_x1,cell_x2,cell_y1,cell_y2,cell_z1,cell_z2)  
     END IF
   END SUBROUTINE get_fields_at_point_store

  SUBROUTINE postsetup_testing
    REAL(num), DIMENSION(3)   :: pos,bvec,evec
    REAL(num), DIMENSION(3,3) :: btens
    REAL(num), DIMENSION(3)   ::  drifts_mu,drifts_vpll,ExB,bdir
    REAL(num)                 :: bdir_dotgradBmag 

    INTEGER :: i,j 

    WRITE (*,*) 'xminmax',x_grid_min_local, x_grid_max_local
    pos(1) =  x_grid_min_local
    pos(2) =  y_grid_min_local
    pos(3) =  z_grid_min_local

    CALL  get_fields_at_point(pos,bvec,evec,btens)
    WRITE (*,*) 'pos',pos(1),pos(2),pos(3)
    WRITE (*,*) 'bvec',bvec(1),bvec(2),bvec(3)
    WRITE (*,*) 'evec',evec(1),evec(2),evec(3)
    DO i=1,3
       DO j=1,3
          WRITE (*,*) btens(i,j)
       END DO
    END DO

    OPEN(unit=24,file='tfields')
    OPEN(unit=25,file='tdrifts')

    Do i=1,nx
       pos(1) =  x_grid_min_local + i*dx*0.2
       CALL  get_fields_at_point(pos,bvec,evec,btens)
       WRITE (24,*) pos(1),bvec(1),bvec(2),bvec(3),evec(1),evec(2),evec(3)
       CALL  get_drifts(pos,Evec,Bvec,Btens,ExB,bdir,drifts_mu,drifts_vpll,bdir_dotgradBmag)    
       !WRITE (24,*) pos(1),bvec(1),bvec(2),bvec(3),btens(1,1),btens(1,2),btens(1,3)
       WRITE (25,*) pos(1),ExB(1),ExB(2),ExB(3),drifts_mu(1),drifts_mu(2),drifts_mu(3)
    END DO
    !CALL MPI_EXIT()
    STOP  
  END SUBROUTINE postsetup_testing

  ! Background distribution function used for delta-f calculations.
  ! Specialise to a drifting (tri)-Maxwellian to simplify and ensure
  ! zero density/current divergence.
  ! Can effectively switch off deltaf method by setting zero background density.

  FUNCTION f0(species, current, mass)
    TYPE(particle), POINTER, INTENT(IN) :: current
    REAL(num), INTENT(IN) :: mass
    REAL(num) :: f0
    REAL(num) :: density
    REAL(num) :: f0_exponent, norm, two_kb_mass, two_pi_kb_mass3
    REAL(num), DIMENSION(3) :: tl, dl
    TYPE(particle_species),  INTENT(INOUT) :: species
    
    IF (species%use_deltaf) THEN
       two_kb_mass = 2.0_num * kb * mass
       two_pi_kb_mass3 = (pi * two_kb_mass)**3
       !DO i=1,3
       !   CALL params_local(current, species%temp(:,:,:,i), &
       !        & species%drift(:,:,:,i), tl(i), dl(i))
       !END DO
       CALL params_local_all(current, species%temp, &
               & species%drift,species%density, tl, dl, density)

       CALL fluid_quants(current%part_pos,density,dl)

       ! Per-particle momentum
       dl = current%mass*dl/density
       density = density/current%mass
       
       !WRITE (*,*) current%part_p(1), dl(1)
  !     f0_exponent = ((current%part_p(1) - dl(1))**2 / tl(1) &
  !                  + (current%part_p(2) - dl(2))**2 / tl(2) &
  !                  + (current%part_p(3) - dl(3))**2 / tl(3)) / two_kb_mass
       !     norm = density / SQRT(two_pi_kb_mass3 * tl(1) * tl(2) * tl(3))
       ! Single-temperature distribution.
       f0_exponent = ((current%part_p(1) - dl(1))**2 / tl(1)) / two_kb_mass
       norm = density / SQRT(two_pi_kb_mass3 * tl(1))
       f0 = norm * EXP(-f0_exponent)
    ELSE
       f0 = 0.0_num
    END IF

  END FUNCTION f0

  SUBROUTINE fluid_quants(r,density,mom)
    REAL(num), DIMENSION(3), INTENT(IN)    :: r
    REAL(num), DIMENSION(3), INTENT(INOUT) :: mom
    REAL(num) :: density
    INTEGER :: ix,ixp,ix0
    REAL(num) :: x,cell_x,xsc
    
    x = r(1)
    xsc = x/dx_fluid
    ix0=floor(xsc)
    ix =modulo(ix0-1,nx_fluid)+1
    ixp = ix+1
    ixp = modulo(ixp-1,nx_fluid)+1
    cell_x = xsc - ix0

    mom=0
    density = (1-cell_x)*dens_fluid(ix) + cell_x*dens_fluid(ixp)
    mom(1)  = (1-cell_x)*p_fluid(ix)    + cell_x*p_fluid(ixp)

  END SUBROUTINE fluid_quants
  
  !Piecewise linear evaluation of density.
  REAL(num) FUNCTION density_fluid(r)
    REAL(num), DIMENSION(3) :: r

    INTEGER :: ix,ixp
    REAL(num) :: x,cell_x
    x = r(1)
    ix =floor(r(1)/dx_fluid)
    ix =modulo(ix-1,nx_fluid)+1
    ixp = ix+1
    ixp = modulo(ixp-1,nx_fluid)+1
    
    cell_x = x - ix*dx_fluid

    density_fluid = (1-cell_x)*dens_fluid(ix) + cell_x*dens_fluid(ixp)
  END FUNCTION density_fluid

  !Piecewise linear evaluation of density.
  !px is the macroparticle momentum
  !weight is the macroparticle mass
  SUBROUTINE assignweight_fluid(r,px,weight,mforce)
    REAL(num), DIMENSION(:), INTENT(INOUT) :: r,mforce
    REAL(num), INTENT(IN)                  :: px,weight
    
    INTEGER :: ix,ixp,ix0
    REAL(num) :: x,cell_x,weight_dV,xsc,dydz

    dydz = (z_max-z_min)*(y_max-y_min)
    
    x = r(1)
    xsc = x/dx_fluid
    ix0=floor(xsc)
    ix =modulo(ix0-1,nx_fluid)+1
    ixp = ix+1
    ixp = modulo(ixp-1,nx_fluid)+1
    
    cell_x = xsc - ix0

    weight_dV = 1.0/(dx_fluid*dydz)
    
    dens_fluid0(ixp) = dens_fluid0(ixp)  + cell_x*weight*weight_dV
    dens_fluid0(ix)  = dens_fluid0(ix) + (1.0_num-cell_x)*weight*weight_dV
    !dens_fluid0(ix)  = dens_fluid0(ix) + weight_dx
    
    p_fluid0(ixp)= p_fluid0(ixp) + px*cell_x*weight_dV
    p_fluid0(ix) = p_fluid0(ix)  + px*(1.0_num-cell_x)*weight_dV
    !p_fluid0(ix) = p_fluid0(ix)  + px*weight_dx
    
    pressuret(ix) = pressuret(ix) + px*px*weight_dV
    forcet(ix) = forcet(ix) + weight_dV*mforce(1)
    
  END SUBROUTINE assignweight_fluid

  SUBROUTINE fluideq_diag
    INTEGER :: ix
    REAL(num) :: xg
    DO ix=1,nx_fluid
       xg = ix*dx_fluid
       WRITE (25,*) time,xg,dens_fluid0(ix),p_fluid0(ix),dens_fluid(ix),p_fluid(ix)
    END DO
  END SUBROUTINE fluideq_diag

  ! Solve a simple fluid equation for moment update.
  ! -lowest order Conservative finite volume 
  SUBROUTINE update_fluideq
    INTEGER :: ix, lr, ip
    REAL(num) :: flux(2), pflux(2), kterm

    !Just solve free-streaming cold gas for now.
    pressuret = 0.0
    forcet = 0.0

    DO ix=1,nx_fluid
       DO lr=1,2
          IF(p_fluid(ix)>0) THEN
             ip = ix+lr-2
          ELSE
             ip = ix+lr-1
          END IF
          ip = MODULO(ip-1,nx_fluid)+1
          flux(lr)  = p_fluid(ip)
          pflux(lr) = p_fluid(ip)**2/dens_fluid(ip) + pressuret(ip)
       END DO
       dens_fluid_next(ix) = dens_fluid(ix) + (flux(1) - flux(2))*dt_fluid/dx_fluid
       kterm = forcet(ix)
       p_fluid_next(ix)    = p_fluid(ix) &
            & + (pflux(1) - pflux(2))*dt_fluid/dx_fluid + kterm
       !p_fluid_next(ix)    = p_fluid(ix)        
    END DO
    p_fluid = p_fluid_next
    dens_fluid = dens_fluid_next
  END SUBROUTINE update_fluideq

  SUBROUTINE setup_fluid
    
    nx_fluid = nx
    dx_fluid = dx
    ALLOCATE(p_fluid(1:nx_fluid))
    ALLOCATE(p_fluid0(1:nx_fluid))
    ALLOCATE(dens_fluid(1:nx_fluid))
    ALLOCATE(dens_fluid0(1:nx_fluid))
    ALLOCATE(p_fluid_next(1:nx_fluid))
    ALLOCATE(dens_fluid_next(1:nx_fluid))
    ALLOCATE(forcet(1:nx_fluid))
    ALLOCATE(pressuret(1:nx_fluid))
    forcet = 0.0
    pressuret = 0.0
    dens_fluid = 0.0
    dens_fluid0= 0.0
    p_fluid_next = 0.0
    p_fluid = 0.0
    p_fluid0= 0.0
    dens_fluid_next = 0.0
    dt_fluid = dt

    CALL load_fluid
   
  END SUBROUTINE setup_fluid
    
  SUBROUTINE load_fluid
    REAL(num), DIMENSION(3) :: tl,dl,mforce
    REAL(num) :: p,part_weight,dens, dydz, volfac, x_pos
    REAL(num), DIMENSION(4) :: gaussp, gaussw
    INTEGER :: ig,ix
    TYPE(particle), POINTER :: first_part
    TYPE(particle), POINTER :: dummy_part
    ALLOCATE(dummy_part)
    
    gaussp = dx_fluid * (0.5 + 0.5 * (/-0.861136, -0.339981, 0.339981, 0.861136/) )
    gaussw = dx_fluid * 0.5 * (/ 0.347855,  0.652145, 0.652145, 0.347855/)

    dydz = (z_max-z_min)*(y_max-y_min)
    volfac = dydz

    ! Pick first species, and first particle from that: assumed 1-species case
    ! and uniform mass, charge.
    
    first_part => species_list%attached_list%head
  
    do ix=1,nx
       do ig = 1,4
          x_pos = ix*dx_fluid + gaussp(ig)
          dummy_part%part_pos = (/ x_pos, 0.0_num, 0.0_num /)
          CALL params_local_all(dummy_part, species_list%temp, &
               & species_list%drift, species_list%density, tl, dl, dens)
          part_weight = first_part%mass*volfac*dens*gaussw(ig)
          !Drift is momentum/c units per particle.
          p = (dl(1))*volfac*dens*gaussw(ig)
          mforce = (/ 0.0_num, 0.0_num, 0.0_num/)
          CALL assignweight_fluid(dummy_part%part_pos,p,part_weight,mforce)
          !WRITE (*,*) x_pos,ig,part_weight,dens,dl(1)
       end do
       !WRITE (*,*) ix,dens_fluid0(ix)
    end do
    dens_fluid = dens_fluid0
    p_fluid = p_fluid0
    
    dens_fluid0= 0.0
    p_fluid0= 0.0
   
  END SUBROUTINE load_fluid
    
    
  SUBROUTINE solve_fluid
    REAL(num) :: t, t_fin, vmax, xg
    INTEGER   :: it, ix, nt_fluid

    CALL setup_fluid
    
    vmax=1.0
    DO ix=1,nx_fluid
       xg = ix*dx_fluid
       dens_fluid(ix) = 1.0_num !+ 0.3_num * cos(xg)
       p_fluid(ix)    = 1.0_num + 0.3_num * cos(xg)
       vmax = max(vmax,abs(p_fluid(ix)/dens_fluid(ix)))
    END DO

    !Courant?
    dt_fluid = 0.5*dx_fluid/vmax

    t_fin = 10.0_num
    nt_fluid = FLOOR(t_fin/dt_fluid) + 1
    nt_fluid = 400
    
    DO it = 1,nt_fluid
       t = it*dt_fluid
       DO ix=1,nx_fluid
          xg = ix*dx_fluid
          WRITE (25,*) t,xg,dens_fluid(ix),p_fluid(ix)
       END DO
       CALL update_fluideq
    end do
    
    STOP
  END SUBROUTINE SOLVE_FLUID

  subroutine initstep_fluid
      p_fluid0    = 0.0
      dens_fluid0 = 0.0
      pressuret   = 0.0
      forcet      = 0.0
  end subroutine initstep_fluid
  
END MODULE particles
