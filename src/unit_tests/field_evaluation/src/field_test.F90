PROGRAM field_test

  USE shared_data
  USE utilities

  IMPLICIT NONE

  REAL(num), DIMENSION(3) :: pos, evec, bvec
  REAL(num) :: xc, xs, yc, ys, zc, zs, soln
  INTEGER :: i, j, k
  REAL(num), PARAMETER :: tolerance = 1e-14_num
  REAL(num), PARAMETER :: ex_fixed = 1.0_num
  REAL(num), PARAMETER :: ey_fixed = 2.0_num
  REAL(num), PARAMETER :: ez_fixed = 3.0_num
  REAL(num), PARAMETER :: bx_fixed = 4.0_num
  REAL(num), PARAMETER :: by_fixed = 5.0_num
  REAL(num), PARAMETER :: bz_fixed = 6.0_num

  ! Initialise grid sizes
  nx = 10
  ny = 10
  nz = 10

  ! For testing, assume field defined on unit square
  dx = 1.0_num / REAL(nx, num)
  dy = 1.0_num / REAL(ny, num)
  dz = 1.0_num / REAL(nz, num)

  ! Initialise grid - these are cell centred variables in EPOCH
  x_grid_min_local = 0.5_num * dx
  y_grid_min_local = 0.5_num * dy
  z_grid_min_local = 0.5_num * dz


  ALLOCATE(ex(-2:nx+3, -2:ny+3, -2:nz+3))
  ALLOCATE(ey(-2:nx+3, -2:ny+3, -2:nz+3))
  ALLOCATE(ez(-2:nx+3, -2:ny+3, -2:nz+3))
  ALLOCATE(bx(-2:nx+3, -2:ny+3, -2:nz+3))
  ALLOCATE(by(-2:nx+3, -2:ny+3, -2:nz+3))
  ALLOCATE(bz(-2:nx+3, -2:ny+3, -2:nz+3))

  ! Uniform fields
  ex = ex_fixed
  ey = ey_fixed
  ez = ez_fixed
  bx = bx_fixed
  by = by_fixed
  bz = bz_fixed

  pos = (/0.27_num, 0.52_num, 0.43_num/)

  CALL get_fields_at_point(pos, bvec, evec)

  ! Test evaluation for a constant field.
  ! Should be reproduced exactly at the particle
  CALL assert(ABS(evec(1) - ex_fixed) / ex_fixed <= tolerance, &
      'Uniform ex calculation failed!')

  CALL assert(ABS(evec(2) - ey_fixed) / ey_fixed <= tolerance, &
      'Uniform ey calculation failed!')

  CALL assert(ABS(evec(3) - ez_fixed) / ez_fixed <= tolerance, &
       'Uniform ez calculation failed!')

  CALL assert(ABS(bvec(1) - bx_fixed) / bx_fixed <= tolerance, &
      'Uniform bx calculation failed!')

  CALL assert(ABS(bvec(2) - by_fixed) / by_fixed <= tolerance, &
      'Uniform by calculation failed!')

  CALL assert(ABS(bvec(3) - bz_fixed) / bz_fixed <= tolerance, &
       'Uniform bz calculation failed!')

  ! Set-up a linearly varying field.
  ! Again, should be reproducible at the particle position
  DO k = -2, nz + 3
    zs = REAL(k, num) * dz
    zc = zs - 0.5_num * dz
    DO j = -2, ny + 3
      ys = REAL(j, num) * dy
      yc = ys - 0.5_num * dy
      DO i = -2, nx + 3
        xs = REAL(i, num) * dx
        xc = xs - 0.5_num * dx
        ex(i,j,k) = linear_field(xs, yc, zc)
        ey(i,j,k) = linear_field(xc, ys, zc)
        ez(i,j,k) = linear_field(xc, yc, zs)
        bx(i,j,k) = linear_field(xc, ys, zs)
        by(i,j,k) = linear_field(xs, yc, zs)
        bz(i,j,k) = linear_field(xs, ys, zc)
      END DO
    END DO
  END DO

  CALL get_fields_at_point(pos, bvec, evec)

  soln = linear_field(pos(1), pos(2), pos(3))

  CALL assert(ABS(evec(1) - soln) / ABS(soln) <= tolerance, &
      'Linear ex calculation failed!')

  CALL assert(ABS(evec(2) - soln) / ABS(soln) <= tolerance, &
      'Linear ey calculation failed!')

  CALL assert(ABS(evec(3) - soln) / ABS(soln) <= tolerance, &
       'Linear ez calculation failed!')

  CALL assert(ABS(bvec(1) - soln) / ABS(soln) <= tolerance, &
      'Linear bx calculation failed!')

  CALL assert(ABS(bvec(2) - soln) / ABS(soln) <= tolerance, &
      'Linear by calculation failed!')

  CALL assert(ABS(bvec(3) - soln) / ABS(soln) <= tolerance, &
       'Linear bz calculation failed!')

  DEALLOCATE(ex, ey, ez, bx, by, bz)

  PRINT*,'All field evauluation tests passed!'

CONTAINS

  SUBROUTINE assert(check, msg)

    LOGICAL, INTENT(IN) :: check
    CHARACTER(LEN=*), INTENT(IN) :: msg

    IF (.NOT. check) THEN
      PRINT*,TRIM(ADJUSTL(msg))
      STOP 1
    END IF

  END SUBROUTINE assert

  ! Field value = a*x + b*y + c*z
  REAL(num) FUNCTION linear_field(x, y, z)

    REAL(num), INTENT(IN) :: x, y, z
    REAL(num), PARAMETER :: a = 0.8_num
    REAL(num), PARAMETER :: b = -2.9_num
    REAL(num), PARAMETER :: c = 1.3_num

    linear_field = a * x + b * y + c * z

  END FUNCTION linear_field

END PROGRAM field_test
