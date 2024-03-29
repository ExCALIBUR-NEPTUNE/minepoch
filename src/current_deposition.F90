MODULE current_deposition

  USE shared_data
  USE utilities

  IMPLICIT NONE

  INTERFACE current_deposition_esirkepov
    MODULE PROCEDURE current_deposition_esirkepov_pos, current_deposition_esirkepov_store
  END INTERFACE current_deposition_esirkepov

  PRIVATE :: calc_stdata
  PRIVATE :: current_deposition_esirkepov_pos, current_deposition_esirkepov_store

  REAL(num), PRIVATE :: idx, idy, idz, iv
  REAL(num), PRIVATE :: i_yz, i_xz, i_xy
  REAL(num), PRIVATE, PARAMETER :: fac = (1.0_num / 24.0_num)**c_ndims

CONTAINS

  !****************************************************************************
  !> \brief Set run-time constants needed for current deposition module
  !****************************************************************************
  SUBROUTINE setup_current_deposition


    idx = 1.0_num / dx
    idy = 1.0_num / dy
    idz = 1.0_num / dz

    iv = idx * idy * idz

    i_yz = idy * idz * fac
    i_xz = idx * idz * fac
    i_xy = idx * idy * fac

  END SUBROUTINE setup_current_deposition



  !****************************************************************************
  !> \brief Calculate current deposition, given position, velocity and charge.
  !>
  !> Current of a macro-particle is calculated as the product of the
  !> velocity, charge and shape function.
  !> Note that this scheme -is not- charge conserving.
  !>
  !> @param[in] pos 3D (x, y, z) position of particle.
  !> @param[in] vel 3D (vx, vy, vz) velocity of particle.
  !> @param[in] chargeweight Total charge of macro-particle.
  !> @param[in,out] jx Array to deposit x-component of current density
  !> @param[in,out] jy Array to deposit y-component of current density
  !> @param[in,out] jz Array to deposit z-component of current density
  !****************************************************************************
  SUBROUTINE current_deposition_simple(pos, vel, chargeweight, jx, jy, jz)

    REAL(num), DIMENSION(3), INTENT(IN) :: pos, vel
    REAL(num), INTENT(IN) :: chargeweight
    REAL(num), INTENT(INOUT), DIMENSION(1-jng:,1-jng:,1-jng:) :: jx, jy, jz
    REAL(num), DIMENSION(sf_min-1:sf_max+1) :: gx, gy, gz
    REAL(num), DIMENSION(sf_min-1:sf_max+1) :: hx, hy, hz
    REAL(num) :: part_x, part_y, part_z
    REAL(num) :: cell_x_r, cell_y_r, cell_z_r
    INTEGER :: cell_x1, cell_x2
    INTEGER :: cell_y1, cell_y2
    INTEGER :: cell_z1, cell_z2
    INTEGER :: dcellx, dcelly, dcellz
    ! The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x, cell_frac_y, cell_frac_z
    REAL(num) :: jxp, jyp, jzp, jxw, jyw, jzw
    REAL(num) :: cf2
    INTEGER :: ix, iy, iz

    part_x = pos(1) - x_grid_min_local
    part_y = pos(2) - y_grid_min_local
    part_z = pos(3) - z_grid_min_local

    gx = 0.0_num
    gy = 0.0_num
    gz = 0.0_num

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
#include "bspline3/gx.inc"

    ! Now redo shifted by half a cell due to grid stagger.
    cell_x2 = FLOOR(cell_x_r)
    cell_frac_x = REAL(cell_x2, num) - cell_x_r + 0.5_num
    cell_x2 = cell_x2 + 1

    cell_y2 = FLOOR(cell_y_r)
    cell_frac_y = REAL(cell_y2, num) - cell_y_r + 0.5_num
    cell_y2 = cell_y2 + 1

    cell_z2 = FLOOR(cell_z_r)
    cell_frac_z = REAL(cell_z2, num) - cell_z_r + 0.5_num
    cell_z2 = cell_z2 + 1

    dcellx = 0.0_num
    dcelly = 0.0_num
    dcellz = 0.0_num
#include "bspline3/hx_dcell.inc"

    jxp = chargeweight * vel(1) * iv * fac
    jyp = chargeweight * vel(2) * iv * fac
    jzp = chargeweight * vel(3) * iv * fac

    DO iz = sf_min, sf_max
      DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          jxw = gz(iz) * gy(iy) * hx(ix)
          jyw = gz(iz) * hy(iy) * gx(ix)
          jzw = hz(iz) * gy(iy) * gx(ix)
          jx(cell_x2+ix, cell_y1+iy, cell_z1+iz) = &
              jx(cell_x2+ix, cell_y1+iy, cell_z1+iz) + jxw * jxp
          jy(cell_x1+ix, cell_y2+iy, cell_z1+iz) = &
              jy(cell_x1+ix, cell_y2+iy, cell_z1+iz) + jyw * jyp
          jz(cell_x1+ix, cell_y1+iy, cell_z2+iz) = &
              jz(cell_x1+ix, cell_y1+iy, cell_z2+iz) + jzw * jzp
        END DO
      END DO
    END DO

  END SUBROUTINE current_deposition_simple



  !****************************************************************************
  !> \brief Calculate current deposition as particle moves from pos0 to pos1
  !>
  !> Uses the charge-conserving scheme of Esirkepov (see reference) to
  !> calculate current due to macro-particle of total charge `chargeweight`.
  !>
  !> @param[in] pos0 3D (x, y, z) initial position of particle
  !> @param[in] pos1 3D (x, y, z) final position of particle
  !> @param[in] chargeweight Total charge of macro-particle
  !> @param[in,out] jx Array to deposit x-component of current density
  !> @param[in,out] jy Array to deposit y-component of current density
  !> @param[in,out] jz Array to deposit z-component of current density
  !>
  !> Reference: (https://doi.org/10.1016/S0010-4655(00)00228-9)
  !****************************************************************************
  SUBROUTINE current_deposition_esirkepov_pos(pos0, pos1, chargeweight, jx, jy, jz)

    REAL(num), DIMENSION(3), INTENT(IN) :: pos0,pos1
    REAL(num), INTENT(IN) :: chargeweight
    REAL(num), INTENT(INOUT), DIMENSION(1-jng:,1-jng:,1-jng:) :: jx, jy, jz
    TYPE(fields_eval_tmps)  :: st0

    CALL calc_stdata(pos0,st0)
    CALL current_deposition_esirkepov_store(st0, pos1, chargeweight, jx, jy, jz)

  END SUBROUTINE current_deposition_esirkepov_pos



  !****************************************************************************
  !> \breif Calculate current deposition as particle moves from st to pos
  !>
  !> Uses the charge-conserving scheme of Esirkepov (see reference) to
  !> calculate current due to macro-particle of total charge `chargeweight`.
  !> Uses pre-computed weight-function information in place for initial
  !> particle position
  !>
  !> @param[in] st TYPE(fields_eval_tmps) containing particle weight-function
  !>               for initial position
  !> @param[in] pos 3D (x, y, z) final position of particle
  !> @param[in] chargeweight Total charge of macro-particle
  !> @param[in,out] jx Array to deposit x-component of current density
  !> @param[in,out] jy Array to deposit y-component of current density
  !> @param[in,out] jz Array to deposit z-component of current density
  !>
  !> Reference: (https://doi.org/10.1016/S0010-4655(00)00228-9)
  !****************************************************************************
  SUBROUTINE current_deposition_esirkepov_store(st, pos, chargeweight, jx, jy, jz)

    REAL(num), DIMENSION(3), INTENT(IN) :: pos
    TYPE(fields_eval_tmps), INTENT(IN)  :: st
    REAL(num), INTENT(IN) :: chargeweight
    REAL(num), INTENT(INOUT), DIMENSION(1-jng:,1-jng:,1-jng:) :: jx, jy, jz
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
    REAL(num), PARAMETER :: third = (1.0_num / 3.0_num)
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
      zfac1 = st%gz(iz) + 0.5_num * hz(iz)
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
          xfac1 = st%gx(ix) + 0.5_num * hx(ix)
          xfac2 = third * hx(ix) + 0.5_num * st%gx(ix)

          wx = hx(ix) * yzfac
          wy = xfac1 * hygz + xfac2 * hyhz
          wz = st%gx(ix) * hzyfac1 + hx(ix) * hzyfac2

          ! This is the bit that actually solves d(rho)/dt = -div(J)
          jxh = jxh - fjx * wx
          jyh(ix) = jyh(ix) - fjy * wy
          jzh(ix, iy) = jzh(ix, iy) - fjz * wz

          jx(cx, cy, cz) = jx(cx, cy, cz) + jxh
          jy(cx, cy, cz) = jy(cx, cy, cz) + jyh(ix)
          jz(cx, cy, cz) = jz(cx, cy, cz) + jzh(ix, iy)
        END DO
      END DO
    END DO

  END SUBROUTINE current_deposition_esirkepov_store



  !****************************************************************************
  !> \brief Utility routine for calculating cell offsets and weights
  !>
  !> @param[in,out] st Structure containing particle-to-grid data
  !> @param[in] pos Position of particle
  !****************************************************************************
  SUBROUTINE calc_stdata(pos,st)
    TYPE(fields_eval_tmps), INTENT(INOUT) :: st
    REAL(num), DIMENSION(3), INTENT(IN) :: pos
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
    gx = 0.0_num
    gy = 0.0_num
    gz = 0.0_num
#include "bspline3/gx.inc"

    !Temporary storage for current deposition.
    st%gx = gx
    st%gy = gy
    st%gz = gz

    st%cell_x1 = cell_x1
    st%cell_y1 = cell_y1
    st%cell_z1 = cell_z1

  END SUBROUTINE calc_stdata

END MODULE current_deposition
