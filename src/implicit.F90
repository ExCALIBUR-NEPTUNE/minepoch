MODULE implicit

  USE iso_c_binding
  USE boundary
  USE calc_df
  USE particles
  USE shared_data

  IMPLICIT NONE

  REAL(num), DIMENSION(:), ALLOCATABLE :: x0
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: rho0, rho
  ! Number of independent variables.
  INTEGER, PARAMETER :: nvar_solve = 6

  PRIVATE
  PUBLIC :: solve_implicit_pic, computef, init_1d_index

CONTAINS

  SUBROUTINE solve_implicit_pic
#ifdef TRILINOS
    LOGICAL :: converged, converged_local
    INTEGER :: ix, iy, iz, row, iters
    INTEGER, PARAMETER :: max_iters = 40
    REAL(C_DOUBLE), DIMENSION(:), ALLOCATABLE :: x, dir, f
    REAL(C_DOUBLE) :: max_delta, ep_conv
    INTEGER :: problem_size, errcode

    problem_size = nx * ny * nz * nvar_solve

    ALLOCATE(x(problem_size), dir(problem_size))
    ALLOCATE(x0(problem_size))
    ALLOCATE(f(problem_size))

    ! If using pseudo current correction to Gauss' law, calculate
    ! initial charge density
    IF (use_pseudo_current) THEN
      ALLOCATE(rho0(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
      ALLOCATE(rho(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
      CALL calc_charge_density(rho0)
      CALL processor_summation_bcs(rho0, ng)
    END IF

    converged = .FALSE.
    iters = 0

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

    ! Store x0
    x0 = x

    CALL computef(x, f)

    ep_conv = nonlinear_tolerance + nonlinear_tolerance * SQRT(DOT_PRODUCT(f, f))

    DO WHILE (.NOT. converged)

      ! Call GMRES to calculate direction
      CALL solve_gmres(x, dir, linear_tolerance)

      ! Update x
      x = x - dir

      CALL computef(x, f)
      max_delta = SQRT(DOT_PRODUCT(f, f))

      converged_local = max_delta < ep_conv
      CALL MPI_ALLREDUCE(converged_local, converged, 1, MPI_LOGICAL, &
          MPI_LAND, comm, errcode)

      iters = iters + 1
      IF (iters > max_iters) THEN
        PRINT*,'Too many iterations'
        PRINT*,max_delta, ep_conv
        STOP
      END IF
    END DO

    IF (verbose_solver) THEN
      PRINT*, '*** Non-linear solver took ', iters, ' iterations.***'
    END IF

    ! Set time centred field for final particle update
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          row = index1d(ix,iy,iz)
          ex(ix,iy,iz) = 0.5_num * (x(row+1) + x0(row+1))
          ey(ix,iy,iz) = 0.5_num * (x(row+2) + x0(row+2))
          ez(ix,iy,iz) = 0.5_num * (x(row+3) + x0(row+3))
          bx(ix,iy,iz) = 0.5_num * (x(row+4) + x0(row+4))
          by(ix,iy,iz) = 0.5_num * (x(row+5) + x0(row+5))
          bz(ix,iy,iz) = 0.5_num * (x(row+6) + x0(row+6))
        END DO
      END DO
    END DO

    ! Call boundary conditions on trial solution
    CALL efield_bcs(ex, ey, ez, ng)
    CALL bfield_final_bcs(bx, by, bz, ng)

    ! Update particle positions
    CALL push_particles(.TRUE.)

    ! Reset field to final value
    CALL unpack_vector(x, ex, ey, ez, bx, by, bz)

    DEALLOCATE(x, dir)
    DEALLOCATE(x0)
    DEALLOCATE(f)
    IF (use_pseudo_current) THEN
      DEALLOCATE(rho0, rho)
    END IF

#endif
  END SUBROUTINE solve_implicit_pic



  SUBROUTINE unpack_vector(x, ex, ey, ez, bx, by, bz)

    REAL(C_DOUBLE), DIMENSION(nvar_solve*nx*ny*nz), INTENT(IN) :: x
    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(OUT) :: ex, ey, ez
    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(OUT) :: bx, by, bz
    INTEGER :: ix, iy, iz, row

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

    CALL efield_bcs(ex, ey, ez, ng)
    CALL bfield_final_bcs(bx, by, bz, ng)

  END SUBROUTINE unpack_vector



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
          ex(ix,iy,iz) = 0.5_num * (x(row+1) + x0(row+1))
          ey(ix,iy,iz) = 0.5_num * (x(row+2) + x0(row+2))
          ez(ix,iy,iz) = 0.5_num * (x(row+3) + x0(row+3))
          bx(ix,iy,iz) = 0.5_num * (x(row+4) + x0(row+4))
          by(ix,iy,iz) = 0.5_num * (x(row+5) + x0(row+5))
          bz(ix,iy,iz) = 0.5_num * (x(row+6) + x0(row+6))
        END DO
      END DO
    END DO

    ! Call boundary conditions on trial solution
    CALL efield_bcs(ex, ey, ez, ng)
    CALL bfield_final_bcs(bx, by, bz, ng)

    ! Now calculate rhs
    CALL calculate_rhs(ex, ey, ez, bx, by, bz, rhs)

    f = (x - x0) / dt - rhs

    DEALLOCATE(rhs)

  END SUBROUTINE computef



  SUBROUTINE calculate_rhs(ex, ey, ez, bx, by, bz, rhs)

    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(IN) :: ex, ey, ez
    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(IN) :: bx, by, bz
    REAL(num), DIMENSION(nx*ny*nz*nvar_solve), INTENT(OUT) :: rhs
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: fcorr
    REAL(num) :: c2idx, c2idy, c2idz, idx, idy, idz
    REAL(num) :: div_e
    INTEGER :: ix, iy, iz, row

    idx = 1.0_num / dx
    idy = 1.0_num / dy
    idz = 1.0_num / dz
    c2idx = c**2 * idx
    c2idy = c**2 * idy
    c2idz = c**2 * idz

    IF (use_pseudo_current) THEN
      CALL push_particles(.FALSE., rho)
      CALL processor_summation_bcs(rho, ng)
    ELSE
      CALL push_particles(.FALSE.)
    END IF

    ! Calculate RHS according to Maxwell's equations
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          row = index1d(ix,iy,iz)
          rhs(row+1) = c2idy * (bz(ix  , iy  , iz  ) - bz(ix  , iy-1, iz  )) &
                     - c2idz * (by(ix  , iy  , iz  ) - by(ix  , iy  , iz-1)) &
                     - jx(ix, iy, iz) / epsilon0
          rhs(row+2) = c2idz * (bx(ix  , iy  , iz  ) - bx(ix  , iy  , iz-1)) &
                     - c2idx * (bz(ix  , iy  , iz  ) - bz(ix-1, iy  , iz  )) &
                     - jy(ix, iy, iz) / epsilon0
          rhs(row+3) = c2idx * (by(ix  , iy  , iz  ) - by(ix-1, iy  , iz  )) &
                     - c2idy * (bx(ix  , iy  , iz  ) - bx(ix  , iy-1, iz  )) &
                     - jz(ix, iy, iz) / epsilon0
          rhs(row+4) = idz * (ey(ix  , iy  , iz+1) - ey(ix  , iy  , iz  )) &
                     - idy * (ez(ix  , iy+1, iz  ) - ez(ix  , iy  , iz  ))
          rhs(row+5) = idx * (ez(ix+1, iy  , iz  ) - ez(ix  , iy  , iz  )) &
                     - idz * (ex(ix  , iy  , iz+1) - ex(ix  , iy  , iz  ))
          rhs(row+6) = idy * (ex(ix  , iy+1, iz  ) - ex(ix  , iy  , iz  )) &
                     - idx * (ey(ix+1, iy  , iz  ) - ey(ix  , iy  , iz  ))
        END DO
      END DO
    END DO

    ! Now loop over cells calculating pseudocurrent
    IF (use_pseudo_current) THEN
      ALLOCATE(fcorr(1:nx+1, 1:ny+1, 1:nz+1))
      fcorr = 0.0_num
      DO iz = 1, nz+1
        DO iy = 1, ny+1
          DO ix = 1, nx+1
            ! E field divergence
            div_e = (ex(ix,iy,iz) - ex(ix-1,iy,iz)) / dx &
                  + (ey(ix,iy,iz) - ey(ix,iy-1,iz)) / dy &
                  + (ez(ix,iy,iz) - ez(ix,iy,iz-1)) / dz
            fcorr(ix,iy,iz) = (div_e &
                - 0.5_num * (rho(ix,iy,iz) + rho0(ix,iy,iz)) / epsilon0) &
                * pseudo_current_fac
          END DO
        END DO
      END DO
      ! Clamp F to zero at boundary
      IF (x_min_boundary) fcorr(1,:,:) = 0.0_num
      IF (x_max_boundary) fcorr(nx+1,:,:) = 0.0_num
      IF (y_min_boundary) fcorr(:,1,:) = 0.0_num
      IF (y_max_boundary) fcorr(:,ny+1,:) = 0.0_num
      IF (z_min_boundary) fcorr(:,:,1) = 0.0_num
      IF (z_max_boundary) fcorr(:,:,nz+1) = 0.0_num
      ! Now add contribution to RHS
      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx
            row = index1d(ix,iy,iz)
            rhs(row+1) = rhs(row+1) + (fcorr(ix+1,iy,iz) - fcorr(ix,iy,iz)) / epsilon0 / dx
            rhs(row+2) = rhs(row+2) + (fcorr(ix,iy+1,iz) - fcorr(ix,iy,iz)) / epsilon0 / dy
            rhs(row+3) = rhs(row+3) + (fcorr(ix,iy,iz+1) - fcorr(ix,iy,iz)) / epsilon0 / dz
          END DO
        END DO
      END DO
      DEALLOCATE(fcorr)
    END IF

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

    index1d = ((ix-1) + (iy-1) * nx + (iz-1) * nx * ny) * nvar_solve

  END FUNCTION index1d

END MODULE implicit
