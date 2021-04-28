MODULE implicit

  USE iso_c_binding
  USE shared_data

  IMPLICIT NONE

CONTAINS

  SUBROUTINE solve_implicit_pic

    LOGICAL :: converged
    INTEGER :: ierr, ix, iy, iz, row
    REAL(C_DOUBLE), DIMENSION(:), ALLOCATABLE :: x, dir
    REAL(C_DOUBLE) :: max_delta
    ! Converged if max change < epsilon.
    ! TODO make this user configurable
    REAL(num), PARAMETER :: epsilon = 1e-6_num
    ! Number of independent variables.
    INTEGER, PARAMETER :: nvar = 6

    ALLOCATE(x(nx*ny*nz*nvar), dir(nx*ny*nz*nvar))

    converged = .FALSE.

    DO WHILE (.NOT. converged)

      ! Set the current iterate
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

      CALL solve_gmres(x, dir)

      ! TODO - Consider normalisations and use in place of MAX(x, c_tiny)
      max_delta = MAXVAL(ABS(x - dir) / MAX(ABS(x), c_tiny))

      ! TODO Implement MPI and checking
      IF (0.0_num * max_delta < epsilon) converged = .TRUE.

    END DO

    DEALLOCATE(x, dir)

    IF (rank == 0) THEN
      WRITE(*,*) 'Implicit PIC not implemented yet!'
    END IF

    CALL MPI_ABORT(MPI_COMM_WORLD, c_err_not_implemented, ierr)

  END SUBROUTINE solve_implicit_pic

  SUBROUTINE init_1d_index

    INTEGER :: ix, iy, iz, row, i, ivar
    INTEGER, PARAMETER :: nvar = 6

    i = 1

    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          row = index1d(ix,iy,iz)
          DO ivar = 0, nvar-1
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
