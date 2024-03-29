MODULE utilities

  USE constants
  USE shared_data

  IMPLICIT NONE

  INTERFACE grow_array
    MODULE PROCEDURE grow_real_array, grow_integer_array, grow_string_array
  END INTERFACE grow_array

  PRIVATE :: grow_real_array, grow_integer_array, grow_string_array

CONTAINS

  ! Need this for general grad B drift, calculating B*, etc.
  ! For the moment, aiming for code that is easy to modify/test rather
  ! than max performance.
  SUBROUTINE calc_Btens(Btens,hdx,hdy,hdz,gdx,gdy,gdz,&
       & cell_x1,cell_x2,cell_y1,cell_y2,cell_z1,cell_z2)

    INTEGER, INTENT(IN) :: cell_x1,cell_x2,cell_y1,cell_y2,cell_z1,cell_z2
    REAL(num), DIMENSION(3,3), INTENT(OUT) :: Btens
    REAL(num), DIMENSION(sf_min-1:sf_max+1,2), INTENT(IN) :: gdx, gdy, gdz
    REAL(num), DIMENSION(sf_min-1:sf_max+1,2), INTENT(IN) :: hdx, hdy, hdz
    INTEGER :: cello_x,cello_y,cello_z
    INTEGER :: xp,yp,zp !Keep track of which derivative is being taken.
    INTEGER :: ii

    Btens = 0.0_num
    ! Calculate grad-B tensor
    DO cello_x = -2,2
      DO cello_y = -2,2
        DO cello_z = -2,2
          DO ii=1,3
            xp = 1+kronecker_delta(ii,1)
            yp = 1+kronecker_delta(ii,2)
            zp = 1+kronecker_delta(ii,3)
            Btens(ii,1) = Btens(ii,1) + &
                &  gdx(cello_x,xp)*hdy(cello_y,yp)*hdz(cello_z,zp)  &
                & *bx(cell_x1+cello_x,cell_y2+cello_y,cell_z2+cello_z)
            Btens(ii,2) = Btens(ii,2) + &
                &  hdx(cello_x,xp)*gdy(cello_y,yp)*hdz(cello_z,zp)  &
                & *by(cell_x2+cello_x,cell_y1+cello_y,cell_z2+cello_z)
            Btens(ii,3) = Btens(ii,3) + &
                &  hdx(cello_x,xp)*hdy(cello_y,yp)*gdz(cello_z,zp)  &
                & *bz(cell_x2+cello_x,cell_y2+cello_y,cell_z1+cello_z)
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE calc_Btens



  PURE INTEGER FUNCTION kronecker_delta(a,b)
    INTEGER, INTENT(IN) :: a,b
    IF (a==b) THEN
       kronecker_delta=1
    ELSE
       kronecker_delta=0
    END IF
  END FUNCTION KRONECKER_DELTA



  ! Find derivative of weight array with respect to position.
  ! This is stored in array(:,2)
  ! array(:,1) is the weight array itself.
  ! cell_frac is negative position in cell relative to centre?
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



  !****************************************************************************
  !> \brief Evaluate Electric and Magnetic fields at a point
  !>
  !> Subroutine for calculating magnetic and electric fields at a given
  !> position, using particle-to-grid weight functions.
  !>
  !> @param[in] pos 3D position (x,y,z) to evaluate fields at
  !> @param[out] bvec Vector magnetic field (Bx, By, Bz) at pos
  !> @param[out] evec Vector magnetic field (Ex, Ey, Ez) at pos
  !> @param[out] st TYPE(fields_eval_tmps) to store weighting values (OPTIONAL)
  !> @param[out] btens REAL(num) (3,3) to store grad B (OPTIONAL)
  !****************************************************************************
  SUBROUTINE get_fields_at_point(pos,bvec,evec,st,btens)

    REAL(num), DIMENSION(3), INTENT(IN) :: pos
    REAL(num), DIMENSION(3),   INTENT(OUT) :: bvec,evec
    TYPE(fields_eval_tmps), INTENT(OUT), OPTIONAL :: st
    REAL(num), DIMENSION(3,3), INTENT(OUT), OPTIONAL :: btens
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
    REAL(num) :: idx, idy, idz
    REAL(num), PARAMETER :: fac = (1.0_num / 24.0_num)**c_ndims

    idx = 1.0_num / dx
    idy = 1.0_num / dy
    idz = 1.0_num / dz

    gx = 0.0_num
    gy = 0.0_num
    gz = 0.0_num

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

    dcellx = 0.0_num
    dcelly = 0.0_num
    dcellz = 0.0_num
    ! NOTE: These weights require an additional multiplication factor!
#include "bspline3/hx_dcell.inc"
    ! These are the electric and magnetic fields interpolated to the
    ! particle position. They have been checked and are correct.
    ! Actually checking this is messy.
#include "bspline3/e_part.inc"
#include "bspline3/b_part.inc"
    Evec(1) = ex_part * fac
    Evec(2) = ey_part * fac
    Evec(3) = ez_part * fac
    Bvec(1) = Bx_part * fac
    Bvec(2) = By_part * fac
    Bvec(3) = Bz_part * fac

    IF (PRESENT(st)) THEN
      ! Temporary storage for further use (e.g. current deposition)
      st%gx = gx
      st%gy = gy
      st%gz = gz

      st%cell_x1 = cell_x1
      st%cell_y1 = cell_y1
      st%cell_z1 = cell_z1
    END IF

    IF(PRESENT(Btens)) THEN
      CALL h_derivs(gdx,gx,0,cell_frac_x,idx)
      CALL h_derivs(gdy,gy,0,cell_frac_y,idy)
      CALL h_derivs(gdz,gz,0,cell_frac_z,idz)
      CALL h_derivs(hdx,hx,dcellx,cell_frac_x,idx)
      CALL h_derivs(hdy,hy,dcelly,cell_frac_y,idy)
      CALL h_derivs(hdz,hz,dcellz,cell_frac_z,idz)
      CALL calc_Btens(Btens,hdx,hdy,hdz,gdx,gdy,gdz, &
          & cell_x1,cell_x2,cell_y1,cell_y2,cell_z1,cell_z2)
    END IF

  END SUBROUTINE get_fields_at_point

  SUBROUTINE grow_real_array(array, idx)

    REAL(num), DIMENSION(:), POINTER :: array
    INTEGER, INTENT(IN) :: idx
    REAL(num), DIMENSION(:), ALLOCATABLE :: tmp_array
    INTEGER :: old_size, new_size, i

    old_size = SIZE(array)
    IF (idx <= old_size) RETURN

    ALLOCATE(tmp_array(old_size))
    DO i = 1, old_size
      tmp_array(i) = array(i)
    ENDDO

    new_size = 2 * old_size
    DEALLOCATE(array)
    ALLOCATE(array(new_size))

    DO i = 1, old_size
      array(i) = tmp_array(i)
    ENDDO

    DEALLOCATE(tmp_array)

  END SUBROUTINE grow_real_array



  SUBROUTINE grow_integer_array(array, idx)

    INTEGER, DIMENSION(:), POINTER :: array
    INTEGER, INTENT(IN) :: idx
    INTEGER, DIMENSION(:), ALLOCATABLE :: tmp_array
    INTEGER :: old_size, new_size, i

    old_size = SIZE(array)
    IF (idx <= old_size) RETURN

    ALLOCATE(tmp_array(old_size))
    DO i = 1, old_size
      tmp_array(i) = array(i)
    ENDDO

    new_size = 2 * old_size
    DEALLOCATE(array)
    ALLOCATE(array(new_size))

    DO i = 1, old_size
      array(i) = tmp_array(i)
    ENDDO

    DEALLOCATE(tmp_array)

  END SUBROUTINE grow_integer_array



  SUBROUTINE grow_string_array(array, idx)

    CHARACTER(LEN=c_max_string_length), DIMENSION(:), POINTER :: array
    INTEGER, INTENT(IN) :: idx
    CHARACTER(LEN=c_max_string_length), DIMENSION(:), ALLOCATABLE :: tmp_array
    INTEGER :: old_size, new_size, i

    old_size = SIZE(array)
    IF (idx <= old_size) RETURN

    ALLOCATE(tmp_array(old_size))
    DO i = 1, old_size
      tmp_array(i) = array(i)
    ENDDO

    new_size = 2 * old_size
    DEALLOCATE(array)
    ALLOCATE(array(new_size))

    DO i = 1, old_size
      array(i) = tmp_array(i)
    ENDDO

    DEALLOCATE(tmp_array)

  END SUBROUTINE grow_string_array

END MODULE utilities
