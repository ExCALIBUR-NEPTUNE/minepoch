MODULE implicit

  USE iso_c_binding
  USE boundary
  USE shared_data

  IMPLICIT NONE

  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: ex0, ey0, ez0, bx0, by0, bz0
  REAL(num), DIMENSION(:), ALLOCATABLE :: x0, rhs0
  ! Number of independent variables.
  INTEGER, PARAMETER :: nvar_solve = 6

  PRIVATE
  PUBLIC :: solve_implicit_pic, computef, init_1d_index

CONTAINS

  SUBROUTINE solve_implicit_pic

    LOGICAL :: converged
    INTEGER :: ix, iy, iz, row
    REAL(C_DOUBLE), DIMENSION(:), ALLOCATABLE :: x, dir
    REAL(C_DOUBLE) :: max_delta
    ! Converged if max change < epsilon.
    ! TODO make this user configurable
    REAL(num), PARAMETER :: epsilon = 1e-6_num
    INTEGER :: problem_size

    problem_size = nx * ny * nz * nvar_solve

    ALLOCATE(x(problem_size), dir(problem_size))
    ALLOCATE(x0(problem_size), rhs0(problem_size))

    ! Allocate and store initial fields
    ALLOCATE(ex0(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
    ALLOCATE(ey0(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
    ALLOCATE(ez0(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
    ALLOCATE(bx0(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
    ALLOCATE(by0(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
    ALLOCATE(bz0(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))

    ex0 = ex
    ey0 = ey
    ez0 = ez
    bx0 = bx
    by0 = by
    bz0 = bz

    converged = .FALSE.

    ! Initialise the iterate, x
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          row = index1d(ix,iy,iz)
          x(row+1) = ex(ix,iy,iz)
          x(row+2) = ey(ix,iy,iz)
          x(row+3) = ez(ix,iy,iz)
          x(row+4) = bx(ix,iy,iz)
          x(row+5) = by(ix,iy,iz)
          x(row+6) = bz(ix,iy,iz)
        END DO
      END DO
    END DO

    ! Store x0, calculate rhs0
    x0 = x

    CALL calculate_rhs(ex, ey, ez, bx, by, bz, rhs0)

    DO WHILE (.NOT. converged)

      ! Call GMRES to calculate direction
      CALL solve_gmres(x, dir)

      ! TODO - Consider normalisations and use in place of MAX(x, c_tiny)
      max_delta = MAXVAL(ABS(x - dir) / MAX(ABS(x), c_tiny))

      ! TODO Implement MPI and checking
      IF (0.0_num * max_delta < epsilon) converged = .TRUE.

      ! Update x
      x = x + dir
    END DO

    DEALLOCATE(x, dir)
    DEALLOCATE(ex0, ey0, ez0, bx0, by0, bz0)
    DEALLOCATE(x0, rhs0)

  END SUBROUTINE solve_implicit_pic



  SUBROUTINE computef(x, f) BIND(c, name='jfnk_computef')

    REAL(C_DOUBLE), DIMENSION(nvar_solve*nx*ny*nz), INTENT(IN) :: x
    REAL(C_DOUBLE), DIMENSION(nvar_solve*nx*ny*nz), INTENT(OUT) :: f
    REAL(num), DIMENSION(:), ALLOCATABLE :: rhs
    INTEGER :: ix, iy, iz, row

    ALLOCATE(rhs(nvar_solve*nx*ny*nz))

    ! Unpack trial vector
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          row = index1d(ix,iy,iz)
          ex(ix,iy,iz) = x(row+1)
          ey(ix,iy,iz) = x(row+2)
          ez(ix,iy,iz) = x(row+3)
          bx(ix,iy,iz) = x(row+4)
          by(ix,iy,iz) = x(row+5)
          bz(ix,iy,iz) = x(row+6)
        END DO
      END DO
    END DO

    ! Call boundary conditions on trial solution
    CALL efield_bcs(ex, ey, ez, ng)
    CALL bfield_final_bcs(bx, by, bz, ng)

    ! Now calculate rhs
    CALL calculate_rhs(ex, ey, ez, bx, by, bz, rhs)

    F = 0.5_num * (x - x0) / dt - 0.5_num*(rhs + rhs0)

    DEALLOCATE(rhs)

  END SUBROUTINE computef



  SUBROUTINE calculate_rhs(ex, ey, ez, bx, by, bz, rhs)

    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(IN) :: ex, ey, ez
    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(IN) :: bx, by, bz
    REAL(num), DIMENSION(nx*ny*nz*nvar_solve), INTENT(OUT) :: rhs
    REAL(num) :: c2idx, c2idy, c2idz, idx, idy, idz
    INTEGER :: ix, iy, iz, row

    idx = 1.0_num / dx
    idy = 1.0_num / dy
    idz = 1.0_num / dz
    c2idx = c**2 * idx
    c2idy = c**2 * idy
    c2idz = c**2 * idz

    ! Need to add j here
    ! CALL particle push

    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          row = index1d(ix,iy,iz)
          rhs(row+1) = c2idy * (bz(ix  , iy  , iz  ) - bz(ix  , iy-1, iz  )) &
                     - c2idz * (by(ix  , iy  , iz  ) - by(ix  , iy  , iz-1))
          rhs(row+2) = c2idz * (bx(ix  , iy  , iz  ) - bx(ix  , iy  , iz-1)) &
                     - c2idx * (bz(ix  , iy  , iz  ) - bz(ix-1, iy  , iz  ))
          rhs(row+3) = c2idx * (by(ix  , iy  , iz  ) - by(ix-1, iy  , iz  )) &
                     - c2idy * (bx(ix  , iy  , iz  ) - bx(ix  , iy-1, iz  ))
          rhs(row+4) = idz * (ey(ix  , iy  , iz+1) - ey(ix  , iy  , iz  )) &
                     - idy * (ez(ix  , iy+1, iz  ) - ez(ix  , iy  , iz  ))
          rhs(row+5) = idx * (ez(ix+1, iy  , iz  ) - ez(ix  , iy  , iz  )) &
                     - idz * (ex(ix  , iy  , iz+1) - ex(ix  , iy  , iz  ))
          rhs(row+6) = idy * (ex(ix  , iy+1, iz  ) - ex(ix  , iy  , iz  )) &
                     - idx * (ey(ix+1, iy  , iz  ) - ey(ix  , iy  , iz  ))
        END DO
      END DO
    END DO

  END SUBROUTINE calculate_rhs



  SUBROUTINE init_1d_index

    INTEGER :: ix, iy, iz, row, i, ivar

    i = 1

    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          row = index1d(ix,iy,iz)
          DO ivar = 0, nvar_solve-1
            linear_index(i) = row + ivar
            i = i + 1
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE init_1d_index



  PURE FUNCTION index1d(ix, iy, iz)

    INTEGER(C_INT) :: index1d
    INTEGER, INTENT(IN) :: ix, iy, iz

    index1d = (ix-1) + (iy-1) * nx + (iz-1) * nx * ny

  END FUNCTION index1d

END MODULE implicit
