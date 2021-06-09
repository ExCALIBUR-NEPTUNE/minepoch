PROGRAM field_test

  USE shared_data
  USE utilities

  IMPLICIT NONE

  REAL(num), DIMENSION(3) :: pos, evec, bvec
  REAL(num), PARAMETER :: fac = (1.0_num / 24.0_num)**c_ndims

  ! Initialise grid sizes
  nx = 10
  ny = 10
  nz = 10

  ! Initialise grid
  x_grid_min_local = 0.0_num
  y_grid_min_local = 0.0_num
  z_grid_min_local = 0.0_num

  ! For testing, assume field defined on unit square
  dx = 1.0_num / REAL(nx, num)
  dy = 1.0_num / REAL(ny, num)
  dz = 1.0_num / REAL(nz, num)

  ALLOCATE(ex(-2:nx+3, -2:ny+3, -2:nz+3))
  ALLOCATE(ey(-2:nx+3, -2:ny+3, -2:nz+3))
  ALLOCATE(ez(-2:nx+3, -2:ny+3, -2:nz+3))
  ALLOCATE(bx(-2:nx+3, -2:ny+3, -2:nz+3))
  ALLOCATE(by(-2:nx+3, -2:ny+3, -2:nz+3))
  ALLOCATE(bz(-2:nx+3, -2:ny+3, -2:nz+3))

  ex = 3.0_num
  ey = 0.0_num
  ez = 0.0_num
  bx = 0.0_num
  by = 0.0_num
  bz = 0.0_num

  pos = (/0.2_num, 0.5_num, 0.4_num/)

  CALL get_fields_at_point(pos, bvec, evec)

  PRINT*,'evec = ', evec
  PRINT*,'evec * fac = ', evec * fac

  DEALLOCATE(ex, ey, ez, bx, by, bz)

END PROGRAM field_test
