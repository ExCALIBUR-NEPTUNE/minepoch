MODULE implicit

  USE iso_c_binding
  USE boundary
  USE shared_data

  IMPLICIT NONE

  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: ex0, ey0, ez0, bx0, by0, bz0
  REAL(num), DIMENSION(:), ALLOCATABLE :: x0, rhs0

  PRIVATE
  PUBLIC :: solve_implicit_pic, computef, init_1d_index

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
    INTEGER :: problem_size

    problem_size = nx * ny * nz * nvar

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

    IF (rank == 0) THEN
      WRITE(*,*) 'Implicit PIC not implemented yet!'
    END IF

    CALL MPI_ABORT(MPI_COMM_WORLD, c_err_not_implemented, ierr)

  END SUBROUTINE solve_implicit_pic



  SUBROUTINE computef(x, f) BIND(c, name='jfnk_computef')

    REAL(C_DOUBLE), DIMENSION(nvar*nx*ny*nz), INTENT(IN) :: x
    REAL(C_DOUBLE), DIMENSION(nvar*nx*ny*nz), INTENT(OUT) :: f
    REAL(num), DIMENSION(:), ALLOCATABLE :: rhs
    INTEGER :: ix, iy, iz, row, ierr

    IF (rank == 0) THEN
      WRITE(*,*) 'ComputeF not implemented yet!'
    END IF

    CALL MPI_ABORT(MPI_COMM_WORLD, c_err_not_implemented, ierr)

    ALLOCATE(rhs(nvar*nx*ny*nz))

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

    ! Fill the residual vector, F
    F = 0.0_num
    DEALLOCATE(rhs)

  END SUBROUTINE computef



  SUBROUTINE calculate_rhs(ex, ey, ez, bx, by, bz, rhs)

    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(IN) :: ex, ey, ez
    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(IN) :: bx, by, bz
    REAL(num), DIMENSION(nx*ny*nz*nvar), INTENT(OUT) :: rhs

    ! TODO FixME
    rhs = 0.0_num

  END SUBROUTINE calculate_rhs



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
