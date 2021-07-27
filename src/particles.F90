MODULE particles

  USE boundary
  USE current_deposition
  USE partlist
  USE deltaf_loader2, ONLY: params_local, params_local_all
  USE utilities

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: push_particles, push_particles_2ndstep
  PUBLIC :: setup_fluid
  PUBLIC :: setup_particle_push

  ! Some numerical factors needed for various particle-fields routines.
  REAL(num), PARAMETER, PRIVATE :: fac = (1.0_num / 24.0_num)**c_ndims
  REAL(num), PRIVATE :: idx, idy, idz

  ! For now, fluid equations for a single species
  ! density is mass density, p_fluid is momentum density
  REAL(num), DIMENSION(:), ALLOCATABLE :: dens_fluid_next, forcet, dens_fluid
  REAL(num), DIMENSION(:), ALLOCATABLE :: dens_fluid0, p_fluid0
  REAL(num), DIMENSION(:), ALLOCATABLE :: p_fluid_next, p_fluid, pressuret
  REAL(num) :: dt_fluid, dx_fluid
  INTEGER :: nx_fluid

CONTAINS

  ! Initialise unvarying module variables
  SUBROUTINE setup_particle_push

    ! Unvarying multiplication factors
    idx = 1.0_num / dx
    idy = 1.0_num / dy
    idz = 1.0_num / dz

    ! Now initialise current depositon module
    CALL setup_current_deposition

  END SUBROUTINE setup_particle_push



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

    INTEGER :: ispecies, isubstep
    TYPE(particle_species), POINTER :: species, next_species

    ! Reset current arrays
    jx = 0.0_num
    jy = 0.0_num
    jz = 0.0_num

    jx_d = 0.0_num
    jy_d = 0.0_num
    jz_d = 0.0_num

    next_species => species_list
    DO ispecies = 1, n_species
       species => next_species
       next_species => species%next

       IF (species%immobile) CYCLE

       IF (species%is_driftkinetic) THEN
          CALL push_particles_dk0(species)
       ELSE
          DO isubstep=1,species%nsubstep
             CALL push_particles_lorentz_split(species)
          END DO
       END IF

    ENDDO

    CALL current_bcs
    CALL particle_bcs

  END SUBROUTINE push_particles



  SUBROUTINE push_particles_lorentz_split(species)
    TYPE(particle_species), INTENT(INOUT) :: species
    ! 2nd order accurate particle pusher using parabolic weighting
    ! on and off the grid. The calculation of J looks rather odd
    ! Since it works by solving d(rho)/dt = div(J) and doing a 1st order
    ! Estimate of rho(t+1.5*dt) rather than calculating J directly
    ! This gives exact charge conservation on the grid

    ! J from a given particle, can be spread over up to 3 cells in
    ! Each direction due to parabolic weighting. We allocate 4 or 5
    ! Cells because the position of the particle at t = t+1.5dt is not
    ! known until later. This part of the algorithm could probably be
    ! Improved, but at the moment, this is just a straight copy of
    ! The core of the PSC algorithm
    INTEGER, PARAMETER :: sf0 = sf_min, sf1 = sf_max

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_ux, part_uy, part_uz
    REAL(num) :: part_q, part_mc, ipart_mc, part_weight

    ! Weighting factors as Eqn 4.77 page 25 of manual
    ! Eqn 4.77 would be written as
    ! F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
    ! Defined at the particle position
    REAL(num), DIMENSION(sf_min-1:sf_max+1) :: gx, gy, gz
    ! P+, P- and Tau variables from Boris1970, page27 of manual
    REAL(num) :: uxp, uxm, uyp, uym, uzp, uzm
    REAL(num) :: tau, taux, tauy, tauz, taux2, tauy2, tauz2

    ! charge to mass ratio modified by normalisation
    REAL(num) :: cmratio, ccmratio

    ! Temporary variables
    REAL(num) :: root, dtfac, gamma
    REAL(num) :: delta_x, delta_y, delta_z
    INTEGER(i8) :: ipart
    TYPE(particle), POINTER :: current, next
    TYPE(fields_eval_tmps) :: st_half
    REAL(num), DIMENSION(3) :: part_pos_t1p5, pos_half, Bvec, Evec
    REAL(num) :: weight_back, part_qfac
    REAL(num), DIMENSION(3) :: force, part_v
    REAL(num) :: idt, dto2, dtco2, idt0, dt_sub

    dt_sub = dt / species%nsubstep

    gx = 0.0_num
    gy = 0.0_num
    gz = 0.0_num

    idt = 1.0_num / dt
    dto2 = dt / 2.0_num
    dtco2 = c * dto2
    dtfac = 0.5_num * dt * fac

    idt = 1.0_num / dt_sub
    idt0 = 1.0_num / dt
    dto2 = dt_sub / 2.0_num
    dtco2 = c * dto2
    dtfac = 0.5_num * dt_sub

    IF (species%solve_fluid) CALL initstep_fluid

    current => species%attached_list%head

    IF (.NOT. particles_uniformly_distributed) THEN
      part_weight = species%weight
    END IF

    !DEC$ VECTOR ALWAYS
    DO ipart = 1, species%attached_list%count
      next => current%next
      IF (particles_uniformly_distributed) THEN
        part_weight = current%weight
      END IF

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
      CALL get_fields_at_point(pos_half,Bvec,Evec,st=st_half)

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

      IF (use_esirkepov) THEN
        ! Now advance to t+1.5dt to calculate current. This is detailed in
        ! the manual between pages 37 and 41. The version coded up looks
        ! completely different to that in the manual, but is equivalent.
        ! Use t+1.5 dt so that can update J to t+dt at 2nd order
        part_pos_t1p5 = current%part_pos + (/ delta_x, delta_y, delta_z /)
        ! Current deposition uses position at t+0.5dt and t+1.5dt, particle
        ! assumed to travel in direct line between these locations. Second order
        ! in time for evaluation of current at t+dt
        CALL current_deposition_esirkepov(st_half, part_pos_t1p5, (part_weight*part_qfac), jx, jy, jz)
      ELSE
        part_v = (/part_ux, part_uy, part_uz/) * c / gamma
        CALL current_deposition_simple(current%part_pos, part_v, part_weight * part_q, jx, jy, jz)
      END IF

      IF (species%solve_fluid) THEN
        part_v = (/part_ux, part_uy,part_uz/)
        force = part_q*(Evec+cross(part_v,Bvec))
        CALL assignweight_fluid(current%part_pos, &
            & (part_weight-weight_back)*current%part_p(1), &
            & (part_weight-weight_back)*current%mass,force )
      END IF

      current => next
    END DO

    IF (species%solve_fluid) then
      CALL fluideq_diag
      CALL update_fluideq
   END IF

 END SUBROUTINE push_particles_lorentz_split



 ! Drift-kinetic push. Because we can't do leapfrog, do an explicit two-step scheme.
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

         CALL get_fields_at_point(current%part_pos,Bvec,Evec,st=st_0,btens=Btens)
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
         IF (use_esirkepov) THEN
           CALL current_deposition_esirkepov(st_0, pos_0, (part_weight*part_qfac), jx_d, jy_d, jz_d)
         ELSE
           CALL current_deposition_simple(pos_h, dRdt, part_weight * part_q, jx_d, jy_d, jz_d)
         END IF

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

         CALL get_fields_at_point(pos_h,Bvec,Evec,st=st_half,btens=Btens)
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
         IF (use_esirkepov) THEN
           CALL current_deposition_esirkepov(pos_0, current%part_pos, (part_weight*part_q), jx_d, jy_d, jz_d)
         ELSE
           pos_h = current%part_pos - dRdt_1 * 0.5_num * dt
           CALL current_deposition_simple(pos_h, dRdt_1, part_weight * part_q, jx_d, jy_d, jz_d)
         END IF

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

    CALL  get_fields_at_point(pos,bvec,evec,btens=btens)
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
       CALL  get_fields_at_point(pos,bvec,evec,btens=btens)
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

    IF (n_species <= 0) RETURN
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
