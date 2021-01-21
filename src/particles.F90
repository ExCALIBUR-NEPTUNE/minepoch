MODULE particles

  USE boundary
  USE partlist

  IMPLICIT NONE

CONTAINS


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
    REAL(num) :: idx, idy, idz
    REAL(num) :: idtyz, idtxz, idtxy
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
    REAL(num), DIMENSION(3,3) :: Btens ! Btens(i,j) derivative of j_th component of B field along i direction.

    TYPE(particle), POINTER :: current, next
    TYPE(particle_species), POINTER :: species, next_species

    jx = 0.0_num
    jy = 0.0_num
    jz = 0.0_num

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

    idtyz = idt * idy * idz * fac
    idtxz = idt * idx * idz * fac
    idtxy = idt * idx * idy * fac

    next_species => species_list
    DO ispecies = 1, n_species
       species => next_species
       next_species => species%next

       IF (species%immobile) CYCLE

       IF (species%is_driftkinetic) THEN
          CALL push_particles_dk
       ELSE
          CALL push_particles_lorentz
       END IF

    ENDDO

    CALL current_bcs
    CALL particle_bcs

  contains

    SUBROUTINE push_particles_lorentz

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

    SUBROUTINE push_particles_dk
      REAL(num), DIMENSION(3) :: dRdt, pos_0, pos_h, dRdt_1
      REAL(num), DIMENSION(3) :: bdir, drifts_mu, drifts_vpll
      REAL(num), DIMENSION(3) :: drifts_ExB
      REAL(num), DIMENSION(3) :: Evec, Bvec
      REAL(num) :: part_mu, part_u, bdotBmag
      REAL(num) :: part_u_0,part_u_h, dudt, dudt_1
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

         part_q   = current%charge
         ! Do nonrelativistic drift-kinetics for the moment.
         part_u   = current%part_p(1)
         part_mu  = current%part_p(2)
 
         IF (current%part_p(3)<0.5_num) THEN
            write (18,*) ipart,current%part_pos(1),current%part_pos(2),current%part_pos(3),part_u,part_mu,current%part_p(3)
         END IF
         !stop

         CALL get_fields_at_point(current%part_pos,Bvec,Evec,Btens)
         CALL get_drifts(current%part_pos,Evec,Bvec,Btens,drifts_ExB,bdir, &
              &  drifts_mu,drifts_vpll, bdotBmag)
         ! Ignore B_perp would also make things easier. Can we do Esirkepov based on guiding centre? (dont see why not)
         ! Probably want to do electrostatic somehow: couple DK particles to ES field only?
         ! Otherwise, need full current description. Curl of shape function for magnetisation current?
         !IF (current%part_p(3)<0.5_num) THEN
         !   write (19,*) bdir(1),bdir(2),bdir(3),(part_mu * bdotBmag / current%mass)
         !END IF
         !dRdt = bdir * part_u &
         !     & + drifts_ExB + drifts_mu*part_mu + drifts_vpll*part_u
         dRdt = bdir * part_u
         dudt = - (part_mu * bdotBmag / current%mass) 
         ! Move particles to half timestep position to first order
         pos_0 = current%part_pos + dRdt * dt
         part_u_0 = part_u + dudt * dt

         !Half step position.
         pos_h = 0.5*(current%part_pos + pos_0)
         part_u_h = 0.5*(part_u + part_u_0)

         CALL get_fields_at_point(pos_h,Bvec,Evec,Btens)
         CALL get_drifts(pos_h,Evec,Bvec,Btens,drifts_ExB,bdir, &
              &  drifts_mu,drifts_vpll, bdotBmag)
         dRdt_1 = bdir * part_u_h
         dudt_1 = - (part_mu * bdotBmag / current%mass) 
         
         ! particle has now finished move to end of timestep, so copy back
         ! into particle array
         part_u = part_u + dudt_1 * dt
         current%part_pos = current%part_pos + dRdt_1 * dt 
         current%part_p(1:2)   = (/ part_u, part_mu /)

         current => next
      ENDDO

    END SUBROUTINE push_particles_dk

  END SUBROUTINE push_particles

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
    Btens = Btens*fac
  END SUBROUTINE calc_Btens

  INTEGER FUNCTION kronecker_delta(a,b)
    INTEGER, INTENT(IN) :: a,b
    IF (a==b) THEN
       kronecker_delta=1
    ELSE
       kronecker_delta=0
    END IF
  END FUNCTION KRONECKER_DELTA

  ! Evaluate fields at a point.
  SUBROUTINE get_fields_at_point(pos,bvec,evec,btens)
    REAL(num), DIMENSION(3),   INTENT(INOUT) :: pos,bvec,evec
    REAL(num), DIMENSION(3,3), INTENT(INOUT) :: btens
    REAL(num), PARAMETER :: fac = (1.0_num / 24.0_num)**c_ndims

    ! Fields at particle location
    REAL(num) :: ex_part, ey_part, ez_part, bx_part, by_part, bz_part
    REAL(num) :: cell_x_r, cell_y_r, cell_z_r
    INTEGER :: cell_x1, cell_x2
    INTEGER :: cell_y1, cell_y2
    INTEGER :: cell_z1, cell_z2
    INTEGER :: dcellx, dcelly, dcellz
    REAL(num) :: idx, idy, idz
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

    idx = 1.0_num / dx
    idy = 1.0_num / dy
    idz = 1.0_num / dz

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
    CALL h_derivs(gdx,gx,0,cell_frac_x,idx)
    CALL h_derivs(gdy,gy,0,cell_frac_y,idy)
    CALL h_derivs(gdz,gz,0,cell_frac_z,idz)

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
    CALL h_derivs(hdx,hx,dcellx,cell_frac_x,idx)
    CALL h_derivs(hdy,hy,dcelly,cell_frac_y,idy)
    CALL h_derivs(hdz,hz,dcellz,cell_frac_z,idz)
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
    ! This is the shape function partition-of-unity factor.
    Evec  = Evec*fac
    Bvec  = Bvec*fac

    CALL calc_Btens(Btens,hdx,hdy,hdz,gdx,gdy,gdz,idx,idy,idz, &
         & cell_x1,cell_x2,cell_y1,cell_y2,cell_z1,cell_z2)  

  END SUBROUTINE get_fields_at_point

  SUBROUTINE postsetup_testing
    REAL(num), DIMENSION(3)   :: pos,bvec,evec
    REAL(num), DIMENSION(3,3) :: btens
    REAL(num), DIMENSION(3)   ::  drifts_mu,drifts_vpll,ExB,bdir
    REAL(num)                 :: bdir_dotgradBmag 

    INTEGER :: i,j !,k
    !    do i=1-ng,nx+ng
    !    do j=1-ng,ny+ng
    !    do k=1-ng,nz+ng
    !       WRITE(*,*) i,j,k,ex(i,j,k),bx(i,j,k)
    !    end do
    ! end do
    !end do
    !    STOP

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

END MODULE particles
