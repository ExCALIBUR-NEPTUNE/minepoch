MODULE boundary

  USE partlist
  USE particle_temperature
  USE laser
  USE mpi_subtype_control

  IMPLICIT NONE

CONTAINS

  SUBROUTINE setup_particle_boundaries

    INTEGER :: i, ierr
    LOGICAL :: error
    CHARACTER(LEN=5), DIMENSION(2*c_ndims) :: &
        boundary = (/ 'x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max' /)

    ! For some types of boundary, fields and particles are treated in
    ! different ways, deal with that here

    DO i = 1, 2*c_ndims
      IF (bc_particle(i) == c_bc_other) bc_particle(i) = c_bc_reflect
      IF (bc_field(i) == c_bc_other) bc_field(i) = c_bc_clamp
      IF (bc_field(i) == c_bc_simple_laser) add_laser(i) = .TRUE.
    ENDDO

    ! Note, for laser bcs to work, the main bcs must be set IN THE CODE to
    ! simple_laser (or outflow) and the field bcs to c_bc_clamp. Particles
    ! can then be set separately. IN THE DECK, laser bcs are chosen either
    ! by seting the main bcs OR by setting the field bcs to simple_laser
    ! (or outflow).

    ! Laser boundaries assume open particles unless otherwise specified.
    DO i = 1, 2*c_ndims
      IF (bc_particle(i) == c_bc_simple_laser &
          .OR. bc_particle(i) == c_bc_simple_outflow) &
              bc_particle(i) = c_bc_open
    ENDDO

    ! Note: reflecting EM boundaries not yet implemented.
    DO i = 1, 2*c_ndims
      IF (bc_field(i) == c_bc_reflect) bc_field(i) = c_bc_clamp
      IF (bc_field(i) == c_bc_open) bc_field(i) = c_bc_simple_outflow
    ENDDO

    ! Sanity check on particle boundaries
    error = .FALSE.
    DO i = 1, 2*c_ndims
      IF (bc_particle(i) == c_bc_periodic &
          .OR. bc_particle(i) == c_bc_reflect &
          .OR. bc_particle(i) == c_bc_thermal &
          .OR. bc_particle(i) == c_bc_open) CYCLE
      IF (rank == 0) THEN
        WRITE(*,*)
        WRITE(*,*) '*** ERROR ***'
        WRITE(*,*) 'Unrecognised particle boundary condition on "', &
            boundary(i), '" boundary.'
      ENDIF
      error = .TRUE.
    ENDDO

    IF (error) CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)

  END SUBROUTINE setup_particle_boundaries



  ! Exchanges field values at processor boundaries and applies field
  ! boundary conditions
  SUBROUTINE field_bc(field, ng)

    INTEGER, INTENT(IN) :: ng
    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(INOUT) :: field

    CALL do_field_mpi_with_lengths(field, ng, nx, ny, nz)

  END SUBROUTINE field_bc



  SUBROUTINE do_field_mpi_with_lengths_slice(field, direction, ng, n1_local, &
      n2_local)

    INTEGER, INTENT(IN) :: direction, ng
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(INOUT) :: field
    INTEGER, INTENT(IN) :: n1_local, n2_local
    INTEGER :: proc1_min, proc1_max
    INTEGER :: proc2_min, proc2_max
    INTEGER, DIMENSION(c_ndims-1) :: sizes, subsizes, starts
    INTEGER :: subarray, basetype, sz, szmax, i, j, n
    REAL(num), ALLOCATABLE :: temp(:)

    basetype = mpireal

    IF (direction == c_dir_x) THEN
      proc1_min = proc_y_min
      proc1_max = proc_y_max
      proc2_min = proc_z_min
      proc2_max = proc_z_max
    ELSE IF (direction == c_dir_y) THEN
      proc1_min = proc_x_min
      proc1_max = proc_x_max
      proc2_min = proc_z_min
      proc2_max = proc_z_max
    ELSE
      proc1_min = proc_x_min
      proc1_max = proc_x_max
      proc2_min = proc_y_min
      proc2_max = proc_y_max
    ENDIF

    sizes(1) = n1_local + 2 * ng
    sizes(2) = n2_local + 2 * ng
    starts = 1

    szmax = sizes(1) * ng
    sz = sizes(2) * ng
    IF (sz > szmax) szmax = sz

    ALLOCATE(temp(szmax))

    subsizes(1) = ng
    subsizes(2) = sizes(2)

    sz = subsizes(1) * subsizes(2)

    subarray = create_2d_array_subtype(basetype, subsizes, sizes, starts)

    CALL MPI_SENDRECV(field(1,1-ng), 1, subarray, proc1_min, &
        tag, temp, sz, basetype, proc1_max, tag, comm, status, errcode)

    IF (proc1_max /= MPI_PROC_NULL) THEN
      n = 1
      DO j = 1-ng, subsizes(2)-ng
      DO i = n1_local+1, subsizes(1)+n1_local
        field(i,j) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_SENDRECV(field(n1_local+1-ng,1-ng), 1, subarray, proc1_max, &
        tag, temp, sz, basetype, proc1_min, tag, comm, status, errcode)

    IF (proc1_min /= MPI_PROC_NULL) THEN
      n = 1
      DO j = 1-ng, subsizes(2)-ng
      DO i = 1-ng, subsizes(1)-ng
        field(i,j) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_TYPE_FREE(subarray, errcode)

    subsizes(1) = sizes(1)
    subsizes(2) = ng

    sz = subsizes(1) * subsizes(2)

    subarray = create_2d_array_subtype(basetype, subsizes, sizes, starts)

    CALL MPI_SENDRECV(field(1-ng,1), 1, subarray, proc2_min, &
        tag, temp, sz, basetype, proc2_max, tag, comm, status, errcode)

    IF (proc2_max /= MPI_PROC_NULL) THEN
      n = 1
      DO j = n2_local+1, subsizes(2)+n2_local
      DO i = 1-ng, subsizes(1)-ng
        field(i,j) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_SENDRECV(field(1-ng,n2_local+1-ng), 1, subarray, proc2_max, &
        tag, temp, sz, basetype, proc2_min, tag, comm, status, errcode)

    IF (proc2_min /= MPI_PROC_NULL) THEN
      n = 1
      DO j = 1-ng, subsizes(2)-ng
      DO i = 1-ng, subsizes(1)-ng
        field(i,j) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_TYPE_FREE(subarray, errcode)

    DEALLOCATE(temp)

  END SUBROUTINE do_field_mpi_with_lengths_slice



  SUBROUTINE do_field_mpi_with_lengths(field, ng, nx_local, ny_local, &
      nz_local)

    INTEGER, INTENT(IN) :: ng
    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(INOUT) :: field
    INTEGER, INTENT(IN) :: nx_local, ny_local, nz_local
    INTEGER, DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER :: subarray, basetype, sz, szmax, i, j, k, n
    REAL(num), ALLOCATABLE :: temp(:)

    basetype = mpireal

    sizes(1) = nx_local + 2 * ng
    sizes(2) = ny_local + 2 * ng
    sizes(3) = nz_local + 2 * ng
    starts = 1

    szmax = sizes(1) * sizes(2) * ng
    sz = sizes(1) * sizes(3) * ng
    IF (sz > szmax) szmax = sz
    sz = sizes(2) * sizes(3) * ng
    IF (sz > szmax) szmax = sz

    ALLOCATE(temp(szmax))

    subsizes(1) = ng
    subsizes(2) = sizes(2)
    subsizes(3) = sizes(3)

    sz = subsizes(1) * subsizes(2) * subsizes(3)

    subarray = create_3d_array_subtype(basetype, subsizes, sizes, starts)

    CALL MPI_SENDRECV(field(1,1-ng,1-ng), 1, subarray, proc_x_min, &
        tag, temp, sz, basetype, proc_x_max, tag, comm, status, errcode)

    IF (proc_x_max /= MPI_PROC_NULL) THEN
      n = 1
      DO k = 1-ng, subsizes(3)-ng
      DO j = 1-ng, subsizes(2)-ng
      DO i = nx_local+1, subsizes(1)+nx_local
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_SENDRECV(field(nx_local+1-ng,1-ng,1-ng), 1, subarray, proc_x_max, &
        tag, temp, sz, basetype, proc_x_min, tag, comm, status, errcode)

    IF (proc_x_min /= MPI_PROC_NULL) THEN
      n = 1
      DO k = 1-ng, subsizes(3)-ng
      DO j = 1-ng, subsizes(2)-ng
      DO i = 1-ng, subsizes(1)-ng
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_TYPE_FREE(subarray, errcode)

    subsizes(1) = sizes(1)
    subsizes(2) = ng
    subsizes(3) = sizes(3)

    sz = subsizes(1) * subsizes(2) * subsizes(3)

    subarray = create_3d_array_subtype(basetype, subsizes, sizes, starts)

    CALL MPI_SENDRECV(field(1-ng,1,1-ng), 1, subarray, proc_y_min, &
        tag, temp, sz, basetype, proc_y_max, tag, comm, status, errcode)

    IF (proc_y_max /= MPI_PROC_NULL) THEN
      n = 1
      DO k = 1-ng, subsizes(3)-ng
      DO j = ny_local+1, subsizes(2)+ny_local
      DO i = 1-ng, subsizes(1)-ng
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_SENDRECV(field(1-ng,ny_local+1-ng,1-ng), 1, subarray, proc_y_max, &
        tag, temp, sz, basetype, proc_y_min, tag, comm, status, errcode)

    IF (proc_y_min /= MPI_PROC_NULL) THEN
      n = 1
      DO k = 1-ng, subsizes(3)-ng
      DO j = 1-ng, subsizes(2)-ng
      DO i = 1-ng, subsizes(1)-ng
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_TYPE_FREE(subarray, errcode)

    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = ng

    sz = subsizes(1) * subsizes(2) * subsizes(3)

    subarray = create_3d_array_subtype(basetype, subsizes, sizes, starts)

    CALL MPI_SENDRECV(field(1-ng,1-ng,1), 1, subarray, proc_z_min, &
        tag, temp, sz, basetype, proc_z_max, tag, comm, status, errcode)

    IF (proc_z_max /= MPI_PROC_NULL) THEN
      n = 1
      DO k = nz_local+1, subsizes(3)+nz_local
      DO j = 1-ng, subsizes(2)-ng
      DO i = 1-ng, subsizes(1)-ng
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_SENDRECV(field(1-ng,1-ng,nz_local+1-ng), 1, subarray, proc_z_max, &
        tag, temp, sz, basetype, proc_z_min, tag, comm, status, errcode)

    IF (proc_z_min /= MPI_PROC_NULL) THEN
      n = 1
      DO k = 1-ng, subsizes(3)-ng
      DO j = 1-ng, subsizes(2)-ng
      DO i = 1-ng, subsizes(1)-ng
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_TYPE_FREE(subarray, errcode)

    DEALLOCATE(temp)

  END SUBROUTINE do_field_mpi_with_lengths



  SUBROUTINE do_field_mpi_with_lengths_r4(field, ng, nx_local, ny_local, &
      nz_local)

    INTEGER, INTENT(IN) :: ng
    REAL(r4), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(INOUT) :: field
    INTEGER, INTENT(IN) :: nx_local, ny_local, nz_local
    INTEGER, DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER :: subarray, basetype, sz, szmax, i, j, k, n
    REAL(r4), ALLOCATABLE :: temp(:)

    basetype = MPI_REAL4

    sizes(1) = nx_local + 2 * ng
    sizes(2) = ny_local + 2 * ng
    sizes(3) = nz_local + 2 * ng
    starts = 1

    szmax = sizes(1) * sizes(2) * ng
    sz = sizes(1) * sizes(3) * ng
    IF (sz > szmax) szmax = sz
    sz = sizes(2) * sizes(3) * ng
    IF (sz > szmax) szmax = sz

    ALLOCATE(temp(szmax))

    subsizes(1) = ng
    subsizes(2) = sizes(2)
    subsizes(3) = sizes(3)

    sz = subsizes(1) * subsizes(2) * subsizes(3)

    subarray = create_3d_array_subtype(basetype, subsizes, sizes, starts)

    CALL MPI_SENDRECV(field(1,1-ng,1-ng), 1, subarray, proc_x_min, &
        tag, temp, sz, basetype, proc_x_max, tag, comm, status, errcode)

    IF (proc_x_max /= MPI_PROC_NULL) THEN
      n = 1
      DO k = 1-ng, subsizes(3)-ng
      DO j = 1-ng, subsizes(2)-ng
      DO i = nx_local+1, subsizes(1)+nx_local
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_SENDRECV(field(nx_local+1-ng,1-ng,1-ng), 1, subarray, proc_x_max, &
        tag, temp, sz, basetype, proc_x_min, tag, comm, status, errcode)

    IF (proc_x_min /= MPI_PROC_NULL) THEN
      n = 1
      DO k = 1-ng, subsizes(3)-ng
      DO j = 1-ng, subsizes(2)-ng
      DO i = 1-ng, subsizes(1)-ng
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_TYPE_FREE(subarray, errcode)

    subsizes(1) = sizes(1)
    subsizes(2) = ng
    subsizes(3) = sizes(3)

    sz = subsizes(1) * subsizes(2) * subsizes(3)

    subarray = create_3d_array_subtype(basetype, subsizes, sizes, starts)

    CALL MPI_SENDRECV(field(1-ng,1,1-ng), 1, subarray, proc_y_min, &
        tag, temp, sz, basetype, proc_y_max, tag, comm, status, errcode)

    IF (proc_y_max /= MPI_PROC_NULL) THEN
      n = 1
      DO k = 1-ng, subsizes(3)-ng
      DO j = ny_local+1, subsizes(2)+ny_local
      DO i = 1-ng, subsizes(1)-ng
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_SENDRECV(field(1-ng,ny_local+1-ng,1-ng), 1, subarray, proc_y_max, &
        tag, temp, sz, basetype, proc_y_min, tag, comm, status, errcode)

    IF (proc_y_min /= MPI_PROC_NULL) THEN
      n = 1
      DO k = 1-ng, subsizes(3)-ng
      DO j = 1-ng, subsizes(2)-ng
      DO i = 1-ng, subsizes(1)-ng
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_TYPE_FREE(subarray, errcode)

    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = ng

    sz = subsizes(1) * subsizes(2) * subsizes(3)

    subarray = create_3d_array_subtype(basetype, subsizes, sizes, starts)

    CALL MPI_SENDRECV(field(1-ng,1-ng,1), 1, subarray, proc_z_min, &
        tag, temp, sz, basetype, proc_z_max, tag, comm, status, errcode)

    IF (proc_z_max /= MPI_PROC_NULL) THEN
      n = 1
      DO k = nz_local+1, subsizes(3)+nz_local
      DO j = 1-ng, subsizes(2)-ng
      DO i = 1-ng, subsizes(1)-ng
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_SENDRECV(field(1-ng,1-ng,nz_local+1-ng), 1, subarray, proc_z_max, &
        tag, temp, sz, basetype, proc_z_min, tag, comm, status, errcode)

    IF (proc_z_min /= MPI_PROC_NULL) THEN
      n = 1
      DO k = 1-ng, subsizes(3)-ng
      DO j = 1-ng, subsizes(2)-ng
      DO i = 1-ng, subsizes(1)-ng
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_TYPE_FREE(subarray, errcode)

    DEALLOCATE(temp)

  END SUBROUTINE do_field_mpi_with_lengths_r4



  SUBROUTINE field_zero_gradient(field, stagger_type, boundary)

    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(INOUT) :: field
    INTEGER, INTENT(IN) :: stagger_type, boundary
    INTEGER :: i, nn

    IF (bc_field(boundary) == c_bc_periodic) RETURN

    IF (boundary == c_bd_x_min .AND. x_min_boundary) THEN
      IF (stagger(c_dir_x,stagger_type)) THEN
        DO i = 1, ng
          field(i-ng,:,:) = field(ng-i,:,:)
        ENDDO
      ELSE
        DO i = 1, ng
          field(i-ng,:,:) = field(ng+1-i,:,:)
        ENDDO
      ENDIF
    ELSE IF (boundary == c_bd_x_max .AND. x_max_boundary) THEN
      nn = nx
      IF (stagger(c_dir_x,stagger_type)) THEN
        DO i = 1, ng
          field(nn+i,:,:) = field(nn-i,:,:)
        ENDDO
      ELSE
        DO i = 1, ng
          field(nn+i,:,:) = field(nn+1-i,:,:)
        ENDDO
      ENDIF

    ELSE IF (boundary == c_bd_y_min .AND. y_min_boundary) THEN
      IF (stagger(c_dir_y,stagger_type)) THEN
        DO i = 1, ng
          field(:,i-ng,:) = field(:,ng-i,:)
        ENDDO
      ELSE
        DO i = 1, ng
          field(:,i-ng,:) = field(:,ng+1-i,:)
        ENDDO
      ENDIF
    ELSE IF (boundary == c_bd_y_max .AND. y_max_boundary) THEN
      nn = ny
      IF (stagger(c_dir_y,stagger_type)) THEN
        DO i = 1, ng
          field(:,nn+i,:) = field(:,nn-i,:)
        ENDDO
      ELSE
        DO i = 1, ng
          field(:,nn+i,:) = field(:,nn+1-i,:)
        ENDDO
      ENDIF

    ELSE IF (boundary == c_bd_z_min .AND. z_min_boundary) THEN
      IF (stagger(c_dir_z,stagger_type)) THEN
        DO i = 1, ng
          field(:,:,i-ng) = field(:,:,ng-i)
        ENDDO
      ELSE
        DO i = 1, ng
          field(:,:,i-ng) = field(:,:,ng+1-i)
        ENDDO
      ENDIF
    ELSE IF (boundary == c_bd_z_max .AND. z_max_boundary) THEN
      nn = nz
      IF (stagger(c_dir_z,stagger_type)) THEN
        DO i = 1, ng
          field(:,:,nn+i) = field(:,:,nn-i)
        ENDDO
      ELSE
        DO i = 1, ng
          field(:,:,nn+i) = field(:,:,nn+1-i)
        ENDDO
      ENDIF
    ENDIF

  END SUBROUTINE field_zero_gradient



  SUBROUTINE field_clamp_zero(field, ng, stagger_type, boundary)

    INTEGER, INTENT(IN) :: ng, stagger_type, boundary
    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(INOUT) :: field
    INTEGER :: i, nn

    IF (bc_field(boundary) == c_bc_periodic) RETURN

    IF (boundary == c_bd_x_min .AND. x_min_boundary) THEN
      IF (stagger(c_dir_x,stagger_type)) THEN
        DO i = 1, ng-1
          field(i-ng,:,:) = -field(ng-i,:,:)
        ENDDO
        field(0,:,:) = 0.0_num
      ELSE
        DO i = 1, ng
          field(i-ng,:,:) = -field(ng+1-i,:,:)
        ENDDO
      ENDIF
    ELSE IF (boundary == c_bd_x_max .AND. x_max_boundary) THEN
      nn = nx
      IF (stagger(c_dir_x,stagger_type)) THEN
        field(nn,:,:) = 0.0_num
        DO i = 1, ng-1
          field(nn+i,:,:) = -field(nn-i,:,:)
        ENDDO
      ELSE
        DO i = 1, ng
          field(nn+i,:,:) = -field(nn+1-i,:,:)
        ENDDO
      ENDIF

    ELSE IF (boundary == c_bd_y_min .AND. y_min_boundary) THEN
      IF (stagger(c_dir_y,stagger_type)) THEN
        DO i = 1, ng-1
          field(:,i-ng,:) = -field(:,ng-i,:)
        ENDDO
        field(:,0,:) = 0.0_num
      ELSE
        DO i = 1, ng
          field(:,i-ng,:) = -field(:,ng+1-i,:)
        ENDDO
      ENDIF
    ELSE IF (boundary == c_bd_y_max .AND. y_max_boundary) THEN
      nn = ny
      IF (stagger(c_dir_y,stagger_type)) THEN
        field(:,nn,:) = 0.0_num
        DO i = 1, ng-1
          field(:,nn+i,:) = -field(:,nn-i,:)
        ENDDO
      ELSE
        DO i = 1, ng
          field(:,nn+i,:) = -field(:,nn+1-i,:)
        ENDDO
      ENDIF

    ELSE IF (boundary == c_bd_z_min .AND. z_min_boundary) THEN
      IF (stagger(c_dir_z,stagger_type)) THEN
        DO i = 1, ng-1
          field(:,:,i-ng) = -field(:,:,ng-i)
        ENDDO
        field(:,:,0) = 0.0_num
      ELSE
        DO i = 1, ng
          field(:,:,i-ng) = -field(:,:,ng+1-i)
        ENDDO
      ENDIF
    ELSE IF (boundary == c_bd_z_max .AND. z_max_boundary) THEN
      nn = nz
      IF (stagger(c_dir_z,stagger_type)) THEN
        field(:,:,nn) = 0.0_num
        DO i = 1, ng-1
          field(:,:,nn+i) = -field(:,:,nn-i)
        ENDDO
      ELSE
        DO i = 1, ng
          field(:,:,nn+i) = -field(:,:,nn+1-i)
        ENDDO
      ENDIF
    ENDIF

  END SUBROUTINE field_clamp_zero



  SUBROUTINE processor_summation_bcs(array, ng, flip_direction)

    INTEGER, INTENT(IN) :: ng
    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(INOUT) :: array
    INTEGER, INTENT(IN), OPTIONAL :: flip_direction
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: temp
    INTEGER, DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER :: subarray, nn, sz, i, flip_dir = 0

    IF (PRESENT(flip_direction)) flip_dir = flip_direction

    sizes(1) = nx + 2 * ng
    sizes(2) = ny + 2 * ng
    sizes(3) = nz + 2 * ng
    starts = 1

    subsizes(1) = ng
    subsizes(2) = sizes(2)
    subsizes(3) = sizes(3)
    nn = nx

    subarray = create_3d_array_subtype(mpireal, subsizes, sizes, starts)

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp(subsizes(1), subsizes(2), subsizes(3)))

    temp = 0.0_num
    CALL MPI_SENDRECV(array(nn+1,1-ng,1-ng), 1, subarray, &
        neighbour( 1,0,0), tag, temp, sz, mpireal, &
        neighbour(-1,0,0), tag, comm, status, errcode)

    ! Deal with reflecting boundaries differently
    IF ((bc_particle(c_bd_x_min) == c_bc_reflect .AND. x_min_boundary)) THEN
      IF (flip_dir == c_dir_x) THEN
        ! Currents get reversed in the direction of the boundary
        DO i = 1, ng-1
          array(i,:,:) = array(i,:,:) - array(-i,:,:)
        ENDDO
      ELSE
        DO i = 1, ng-1
          array(i,:,:) = array(i,:,:) + array(1-i,:,:)
        ENDDO
      ENDIF
    ELSE
      array(1:ng,:,:) = array(1:ng,:,:) + temp
    ENDIF

    temp = 0.0_num
    CALL MPI_SENDRECV(array(1-ng,1-ng,1-ng), 1, subarray, &
        neighbour(-1,0,0), tag, temp, sz, mpireal, &
        neighbour( 1,0,0), tag, comm, status, errcode)

    ! Deal with reflecting boundaries differently
    IF ((bc_particle(c_bd_x_max) == c_bc_reflect .AND. x_max_boundary)) THEN
      IF (flip_dir == c_dir_x) THEN
        ! Currents get reversed in the direction of the boundary
        DO i = 1, ng
          array(nn-i,:,:) = array(nn-i,:,:) - array(nn+i,:,:)
        ENDDO
      ELSE
        DO i = 1, ng
          array(nn+1-i,:,:) = array(nn+1-i,:,:) + array(nn+i,:,:)
        ENDDO
      ENDIF
    ELSE
      array(nn+1-ng:nn,:,:) = array(nn+1-ng:nn,:,:) + temp
    ENDIF

    DEALLOCATE(temp)
    CALL MPI_TYPE_FREE(subarray, errcode)

    subsizes(1) = sizes(1)
    subsizes(2) = ng
    subsizes(3) = sizes(3)
    nn = ny

    subarray = create_3d_array_subtype(mpireal, subsizes, sizes, starts)

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp(subsizes(1), subsizes(2), subsizes(3)))

    temp = 0.0_num
    CALL MPI_SENDRECV(array(1-ng,nn+1,1-ng), 1, subarray, &
        neighbour(0, 1,0), tag, temp, sz, mpireal, &
        neighbour(0,-1,0), tag, comm, status, errcode)

    ! Deal with reflecting boundaries differently
    IF ((bc_particle(c_bd_y_min) == c_bc_reflect .AND. y_min_boundary)) THEN
      IF (flip_dir == c_dir_y) THEN
        ! Currents get reversed in the direction of the boundary
        DO i = 1, ng-1
          array(:,i,:) = array(:,i,:) - array(:,-i,:)
        ENDDO
      ELSE
        DO i = 1, ng-1
          array(:,i,:) = array(:,i,:) + array(:,1-i,:)
        ENDDO
      ENDIF
    ELSE
      array(:,1:ng,:) = array(:,1:ng,:) + temp
    ENDIF

    temp = 0.0_num
    CALL MPI_SENDRECV(array(1-ng,1-ng,1-ng), 1, subarray, &
        neighbour(0,-1,0), tag, temp, sz, mpireal, &
        neighbour(0, 1,0), tag, comm, status, errcode)

    ! Deal with reflecting boundaries differently
    IF ((bc_particle(c_bd_y_max) == c_bc_reflect .AND. y_max_boundary)) THEN
      IF (flip_dir == c_dir_y) THEN
        ! Currents get reversed in the direction of the boundary
        DO i = 1, ng
          array(:,nn-i,:) = array(:,nn-i,:) - array(:,nn+i,:)
        ENDDO
      ELSE
        DO i = 1, ng
          array(:,nn+1-i,:) = array(:,nn+1-i,:) + array(:,nn+i,:)
        ENDDO
      ENDIF
    ELSE
      array(:,nn+1-ng:nn,:) = array(:,nn+1-ng:nn,:) + temp
    ENDIF

    DEALLOCATE(temp)
    CALL MPI_TYPE_FREE(subarray, errcode)

    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = ng
    nn = nz

    subarray = create_3d_array_subtype(mpireal, subsizes, sizes, starts)

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp(subsizes(1), subsizes(2), subsizes(3)))

    temp = 0.0_num
    CALL MPI_SENDRECV(array(1-ng,1-ng,nn+1), 1, subarray, &
        neighbour(0,0, 1), tag, temp, sz, mpireal, &
        neighbour(0,0,-1), tag, comm, status, errcode)

    ! Deal with reflecting boundaries differently
    IF ((bc_particle(c_bd_z_min) == c_bc_reflect .AND. z_min_boundary)) THEN
      IF (flip_dir == c_dir_z) THEN
        ! Currents get reversed in the direction of the boundary
        DO i = 1, ng-1
          array(:,:,i) = array(:,:,i) - array(:,:,-i)
        ENDDO
      ELSE
        DO i = 1, ng-1
          array(:,:,i) = array(:,:,i) + array(:,:,1-i)
        ENDDO
      ENDIF
    ELSE
      array(:,:,1:ng) = array(:,:,1:ng) + temp
    ENDIF

    temp = 0.0_num
    CALL MPI_SENDRECV(array(1-ng,1-ng,1-ng), 1, subarray, &
        neighbour(0,0,-1), tag, temp, sz, mpireal, &
        neighbour(0,0, 1), tag, comm, status, errcode)

    ! Deal with reflecting boundaries differently
    IF ((bc_particle(c_bd_z_max) == c_bc_reflect .AND. z_max_boundary)) THEN
      IF (flip_dir == c_dir_z) THEN
        ! Currents get reversed in the direction of the boundary
        DO i = 1, ng
          array(:,:,nn-i) = array(:,:,nn-i) - array(:,:,nn+i)
        ENDDO
      ELSE
        DO i = 1, ng
          array(:,:,nn+1-i) = array(:,:,nn+1-i) + array(:,:,nn+i)
        ENDDO
      ENDIF
    ELSE
      array(:,:,nn+1-ng:nn) = array(:,:,nn+1-ng:nn) + temp
    ENDIF

    DEALLOCATE(temp)
    CALL MPI_TYPE_FREE(subarray, errcode)

    CALL field_bc(array, ng)

  END SUBROUTINE processor_summation_bcs



  SUBROUTINE efield_bcs(ex, ey, ez, ng)

    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(INOUT) :: ex, ey, ez
    INTEGER, INTENT(IN) :: ng
    INTEGER :: i

    ! These are the MPI boundaries
    CALL field_bc(ex, ng)
    CALL field_bc(ey, ng)
    CALL field_bc(ez, ng)

    ! Perfectly conducting boundaries
    DO i = c_bd_x_min, c_bd_x_max, c_bd_x_max - c_bd_x_min
      IF (bc_field(i) == c_bc_conduct) THEN
        CALL field_clamp_zero(ey, ng, c_stagger_ey, i)
        CALL field_clamp_zero(ez, ng, c_stagger_ez, i)
      ENDIF
    ENDDO

    DO i = c_bd_y_min, c_bd_y_max, c_bd_y_max - c_bd_y_min
      IF (bc_field(i) == c_bc_conduct) THEN
        CALL field_clamp_zero(ex, ng, c_stagger_ex, i)
        CALL field_clamp_zero(ez, ng, c_stagger_ez, i)
      ENDIF
    ENDDO

    DO i = c_bd_z_min, c_bd_z_max, c_bd_z_max - c_bd_z_min
      IF (bc_field(i) == c_bc_conduct) THEN
        CALL field_clamp_zero(ex, ng, c_stagger_ex, i)
        CALL field_clamp_zero(ez, ng, c_stagger_ey, i)
      ENDIF
    ENDDO

    DO i = 1, 2*c_ndims
      ! These apply zero field boundary conditions on the edges
      IF (bc_field(i) == c_bc_clamp &
          .OR. bc_field(i) == c_bc_simple_laser &
          .OR. bc_field(i) == c_bc_simple_outflow) THEN
        CALL field_clamp_zero(ex, ng, c_stagger_ex, i)
        CALL field_clamp_zero(ey, ng, c_stagger_ey, i)
        CALL field_clamp_zero(ez, ng, c_stagger_ez, i)
      ENDIF

      ! These apply zero gradient boundary conditions on the edges
      IF (bc_field(i) == c_bc_zero_gradient) THEN
        CALL field_zero_gradient(ex, c_stagger_ex, i)
        CALL field_zero_gradient(ey, c_stagger_ey, i)
        CALL field_zero_gradient(ez, c_stagger_ez, i)
      ENDIF
    ENDDO

  END SUBROUTINE efield_bcs



  SUBROUTINE bfield_bcs(bx, by, bz, ng, mpi_only)

    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(INOUT) :: bx, by, bz
    INTEGER, INTENT(IN) :: ng
    LOGICAL, INTENT(IN) :: mpi_only
    INTEGER :: i

    ! These are the MPI boundaries
    CALL field_bc(bx, ng)
    CALL field_bc(by, ng)
    CALL field_bc(bz, ng)

    IF (mpi_only) RETURN

    ! Perfectly conducting boundaries
    DO i = c_bd_x_min, c_bd_x_max, c_bd_x_max - c_bd_x_min
      IF (bc_field(i) == c_bc_conduct) THEN
        CALL field_clamp_zero(bx, ng, c_stagger_bx, i)
        CALL field_zero_gradient(by, c_stagger_by, i)
        CALL field_zero_gradient(bz, c_stagger_bz, i)
      ENDIF
    ENDDO

    DO i = c_bd_y_min, c_bd_y_max, c_bd_y_max - c_bd_y_min
      IF (bc_field(i) == c_bc_conduct) THEN
        CALL field_clamp_zero(by, ng, c_stagger_by, i)
        CALL field_zero_gradient(bx, c_stagger_bx, i)
        CALL field_zero_gradient(bz, c_stagger_bz, i)
      ENDIF
    ENDDO

    DO i = c_bd_z_min, c_bd_z_max, c_bd_z_max - c_bd_z_min
      IF (bc_field(i) == c_bc_conduct) THEN
        CALL field_clamp_zero(bz, ng, c_stagger_bz, i)
        CALL field_zero_gradient(bx, c_stagger_bx, i)
        CALL field_zero_gradient(by, c_stagger_by, i)
      ENDIF
    ENDDO

    DO i = 1, 2*c_ndims
      ! These apply zero field boundary conditions on the edges
      IF (bc_field(i) == c_bc_clamp &
          .OR. bc_field(i) == c_bc_simple_laser &
          .OR. bc_field(i) == c_bc_simple_outflow) THEN
        CALL field_clamp_zero(bx, ng, c_stagger_bx, i)
        CALL field_clamp_zero(by, ng, c_stagger_by, i)
        CALL field_clamp_zero(bz, ng, c_stagger_bz, i)
      ENDIF

      ! These apply zero gradient boundary conditions on the edges
      IF (bc_field(i) == c_bc_zero_gradient) THEN
        CALL field_zero_gradient(bx, c_stagger_bx, i)
        CALL field_zero_gradient(by, c_stagger_by, i)
        CALL field_zero_gradient(bz, c_stagger_bz, i)
      ENDIF
    ENDDO

  END SUBROUTINE bfield_bcs



  SUBROUTINE bfield_final_bcs(bx, by, bz, ng)

    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(INOUT) :: bx, by, bz
    INTEGER, INTENT(IN) :: ng
    INTEGER :: i

    CALL bfield_bcs(bx, by, bz, ng, .FALSE.)

    IF (x_min_boundary) THEN
      i = c_bd_x_min
      IF (add_laser(i) .OR. bc_field(i) == c_bc_simple_outflow) &
          CALL outflow_bcs_x_min
    ENDIF

    IF (x_max_boundary) THEN
      i = c_bd_x_max
      IF (add_laser(i) .OR. bc_field(i) == c_bc_simple_outflow) &
          CALL outflow_bcs_x_max
    ENDIF

    IF (y_min_boundary) THEN
      i = c_bd_y_min
      IF (add_laser(i) .OR. bc_field(i) == c_bc_simple_outflow) &
          CALL outflow_bcs_y_min
    ENDIF

    IF (y_max_boundary) THEN
      i = c_bd_y_max
      IF (add_laser(i) .OR. bc_field(i) == c_bc_simple_outflow) &
          CALL outflow_bcs_y_max
    ENDIF

    IF (z_min_boundary) THEN
      i = c_bd_z_min
      IF (add_laser(i) .OR. bc_field(i) == c_bc_simple_outflow) &
          CALL outflow_bcs_z_min
    ENDIF

    IF (z_max_boundary) THEN
      i = c_bd_z_max
      IF (add_laser(i) .OR. bc_field(i) == c_bc_simple_outflow) &
          CALL outflow_bcs_z_max
    ENDIF

    CALL bfield_bcs(bx, by, bz, ng, .TRUE.)

  END SUBROUTINE bfield_final_bcs



  SUBROUTINE particle_bcs

    TYPE(particle), POINTER :: cur, next
    TYPE(particle_list), DIMENSION(-1:1,-1:1,-1:1) :: send, recv
    INTEGER :: xbd, ybd, zbd
    INTEGER(i8) :: ixp, iyp, izp
    LOGICAL :: out_of_bounds
    INTEGER :: ispecies, i, ix, iy, iz
    INTEGER :: cell_x, cell_y, cell_z
    REAL(num), DIMENSION(-1:1) :: gx, gy, gz
    REAL(num) :: cell_x_r, cell_frac_x
    REAL(num) :: cell_y_r, cell_frac_y
    REAL(num) :: cell_z_r, cell_frac_z
    REAL(num) :: cf2, temp(3)
    REAL(num) :: part_pos
    TYPE(particle_species), POINTER :: species, next_species

    next_species => species_list
    DO ispecies = 1, n_species
      species => next_species
      next_species => species%next

      cur => species%attached_list%head

      DO iz = -1, 1
        DO iy = -1, 1
          DO ix = -1, 1
            IF (ABS(ix) + ABS(iy) + ABS(iz) == 0) CYCLE
            CALL create_empty_partlist(send(ix, iy, iz))
            CALL create_empty_partlist(recv(ix, iy, iz))
          ENDDO
        ENDDO
      ENDDO

      DO WHILE (ASSOCIATED(cur))
        next => cur%next

        xbd = 0
        ybd = 0
        zbd = 0
        out_of_bounds = .FALSE.

        part_pos = cur%part_pos(1)
        ! Particle has left this processor
        IF (part_pos < x_min_local) THEN
          xbd = -1
          ! Particle has left the system
          IF (x_min_boundary) THEN
            xbd = 0
            IF (bc_particle(c_bd_x_min) == c_bc_reflect) THEN
              cur%part_pos(1) = 2.0_num * x_min - part_pos
              cur%part_p(1) = -cur%part_p(1)
            ELSE IF (bc_particle(c_bd_x_min) == c_bc_thermal) THEN
              ! Always use the triangle particle weighting for simplicity
              cell_y_r = (cur%part_pos(2) - y_grid_min_local) / dy
              cell_y = FLOOR(cell_y_r + 0.5_num)
              cell_frac_y = REAL(cell_y, num) - cell_y_r
              cell_y = cell_y + 1

              cell_z_r = (cur%part_pos(3) - z_grid_min_local) / dz
              cell_z = FLOOR(cell_z_r + 0.5_num)
              cell_frac_z = REAL(cell_z, num) - cell_z_r
              cell_z = cell_z + 1

              cf2 = cell_frac_y**2
              gy(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_y)
              gy( 0) = 0.75_num - cf2
              gy( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_y)

              cf2 = cell_frac_z**2
              gz(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_z)
              gz( 0) = 0.75_num - cf2
              gz( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_z)

              DO i = 1, 3
                temp(i) = 0.0_num
                DO iz = -1, 1
                  DO iy = -1, 1
                    temp(i) = temp(i) + gy(iy) * gz(iz) &
                        * species%ext_temp_x_min(cell_y+iy, cell_z+iz, i)
                  ENDDO
                ENDDO
              ENDDO

              ! x-direction
              i = 1
              cur%part_p(i) = ABS(momentum_from_temperature(&
                  species%mass, temp(i), 0.0_num))

              ! y-direction
              i = 2
              cur%part_p(i) = momentum_from_temperature(&
                  species%mass, temp(i), 0.0_num)

              ! z-direction
              i = 3
              cur%part_p(i) = momentum_from_temperature(&
                  species%mass, temp(i), 0.0_num)

              cur%part_pos(1) = 2.0_num * x_min - part_pos

            ELSE IF (bc_particle(c_bd_x_min) == c_bc_periodic) THEN
              xbd = -1
              cur%part_pos(1) = part_pos + length_x
            ELSE
              ! Default to open boundary conditions - remove particle
              out_of_bounds = .TRUE.
            ENDIF
          ENDIF
        ENDIF

        ! Particle has left this processor
        IF (part_pos >= x_max_local) THEN
          xbd = 1
          ! Particle has left the system
          IF (x_max_boundary) THEN
            xbd = 0
            IF (bc_particle(c_bd_x_max) == c_bc_reflect) THEN
              cur%part_pos(1) = 2.0_num * x_max - part_pos
              cur%part_p(1) = -cur%part_p(1)
            ELSE IF (bc_particle(c_bd_x_max) == c_bc_thermal) THEN
              ! Always use the triangle particle weighting for simplicity
              cell_y_r = (cur%part_pos(2) - y_grid_min_local) / dy
              cell_y = FLOOR(cell_y_r + 0.5_num)
              cell_frac_y = REAL(cell_y, num) - cell_y_r
              cell_y = cell_y + 1

              cell_z_r = (cur%part_pos(3) - z_grid_min_local) / dz
              cell_z = FLOOR(cell_z_r + 0.5_num)
              cell_frac_z = REAL(cell_z, num) - cell_z_r
              cell_z = cell_z + 1

              cf2 = cell_frac_y**2
              gy(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_y)
              gy( 0) = 0.75_num - cf2
              gy( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_y)

              cf2 = cell_frac_z**2
              gz(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_z)
              gz( 0) = 0.75_num - cf2
              gz( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_z)

              DO i = 1, 3
                temp(i) = 0.0_num
                DO iz = -1, 1
                  DO iy = -1, 1
                    temp(i) = temp(i) + gy(iy) * gz(iz) &
                        * species%ext_temp_x_max(cell_y+iy, cell_z+iz, i)
                  ENDDO
                ENDDO
              ENDDO

              ! x-direction
              i = 1
              cur%part_p(i) = -ABS(momentum_from_temperature(&
                  species%mass, temp(i), 0.0_num))

              ! y-direction
              i = 2
              cur%part_p(i) = momentum_from_temperature(&
                  species%mass, temp(i), 0.0_num)

              ! z-direction
              i = 3
              cur%part_p(i) = momentum_from_temperature(&
                  species%mass, temp(i), 0.0_num)

              cur%part_pos(1) = 2.0_num * x_max - part_pos

            ELSE IF (bc_particle(c_bd_x_max) == c_bc_periodic) THEN
              xbd = 1
              cur%part_pos(1) = part_pos - length_x
            ELSE
              ! Default to open boundary conditions - remove particle
              out_of_bounds = .TRUE.
            ENDIF
          ENDIF
        ENDIF

        part_pos = cur%part_pos(2)
        ! Particle has left this processor
        IF (part_pos < y_min_local) THEN
          ybd = -1
          ! Particle has left the system
          IF (y_min_boundary) THEN
            ybd = 0
            IF (bc_particle(c_bd_y_min) == c_bc_reflect) THEN
              cur%part_pos(2) = 2.0_num * y_min - part_pos
              cur%part_p(2) = -cur%part_p(2)
            ELSE IF (bc_particle(c_bd_y_min) == c_bc_thermal) THEN
              ! Always use the triangle particle weighting for simplicity
              cell_x_r = (cur%part_pos(1) - x_grid_min_local) / dx
              cell_x = FLOOR(cell_x_r + 0.5_num)
              cell_frac_x = REAL(cell_x, num) - cell_x_r
              cell_x = cell_x + 1

              cell_z_r = (cur%part_pos(3) - z_grid_min_local) / dz
              cell_z = FLOOR(cell_z_r + 0.5_num)
              cell_frac_z = REAL(cell_z, num) - cell_z_r
              cell_z = cell_z + 1

              cf2 = cell_frac_x**2
              gx(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_x)
              gx( 0) = 0.75_num - cf2
              gx( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_x)

              cf2 = cell_frac_z**2
              gz(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_z)
              gz( 0) = 0.75_num - cf2
              gz( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_z)

              DO i = 1, 3
                temp(i) = 0.0_num
                DO iz = -1, 1
                  DO ix = -1, 1
                    temp(i) = temp(i) + gx(ix) * gz(iz) &
                        * species%ext_temp_y_min(cell_x+ix, cell_z+iz, i)
                  ENDDO
                ENDDO
              ENDDO

              ! x-direction
              i = 1
              cur%part_p(i) = momentum_from_temperature(&
                  species%mass, temp(i), 0.0_num)

              ! y-direction
              i = 2
              cur%part_p(i) = ABS(momentum_from_temperature(&
                  species%mass, temp(i), 0.0_num))

              ! z-direction
              i = 3
              cur%part_p(i) = momentum_from_temperature(&
                  species%mass, temp(i), 0.0_num)

              cur%part_pos(2) = 2.0_num * y_min - part_pos

            ELSE IF (bc_particle(c_bd_y_min) == c_bc_periodic) THEN
              ybd = -1
              cur%part_pos(2) = part_pos + length_y
            ELSE
              ! Default to open boundary conditions - remove particle
              out_of_bounds = .TRUE.
            ENDIF
          ENDIF
        ENDIF

        ! Particle has left this processor
        IF (part_pos >= y_max_local) THEN
          ybd = 1
          ! Particle has left the system
          IF (y_max_boundary) THEN
            ybd = 0
            IF (bc_particle(c_bd_y_max) == c_bc_reflect) THEN
              cur%part_pos(2) = 2.0_num * y_max - part_pos
              cur%part_p(2) = -cur%part_p(2)
            ELSE IF (bc_particle(c_bd_y_max) == c_bc_thermal) THEN
              ! Always use the triangle particle weighting for simplicity
              cell_x_r = (cur%part_pos(1) - x_grid_min_local) / dx
              cell_x = FLOOR(cell_x_r + 0.5_num)
              cell_frac_x = REAL(cell_x, num) - cell_x_r
              cell_x = cell_x + 1

              cell_z_r = (cur%part_pos(3) - z_grid_min_local) / dz
              cell_z = FLOOR(cell_z_r + 0.5_num)
              cell_frac_z = REAL(cell_z, num) - cell_z_r
              cell_z = cell_z + 1

              cf2 = cell_frac_x**2
              gx(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_x)
              gx( 0) = 0.75_num - cf2
              gx( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_x)

              cf2 = cell_frac_z**2
              gz(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_z)
              gz( 0) = 0.75_num - cf2
              gz( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_z)

              DO i = 1, 3
                temp(i) = 0.0_num
                DO iz = -1, 1
                  DO ix = -1, 1
                    temp(i) = temp(i) + gx(ix) * gz(iz) &
                        * species%ext_temp_y_max(cell_x+ix, cell_z+iz, i)
                  ENDDO
                ENDDO
              ENDDO

              ! x-direction
              i = 1
              cur%part_p(i) = momentum_from_temperature(&
                  species%mass, temp(i), 0.0_num)

              ! y-direction
              i = 2
              cur%part_p(i) = -ABS(momentum_from_temperature(&
                  species%mass, temp(i), 0.0_num))

              ! z-direction
              i = 3
              cur%part_p(i) = momentum_from_temperature(&
                  species%mass, temp(i), 0.0_num)

              cur%part_pos(2) = 2.0_num * y_max - part_pos

            ELSE IF (bc_particle(c_bd_y_max) == c_bc_periodic) THEN
              ybd = 1
              cur%part_pos(2) = part_pos - length_y
            ELSE
              ! Default to open boundary conditions - remove particle
              out_of_bounds = .TRUE.
            ENDIF
          ENDIF
        ENDIF

        part_pos = cur%part_pos(3)
        ! Particle has left this processor
        IF (part_pos < z_min_local) THEN
          zbd = -1
          ! Particle has left the system
          IF (z_min_boundary) THEN
            zbd = 0
            IF (bc_particle(c_bd_z_min) == c_bc_reflect) THEN
              cur%part_pos(3) = 2.0_num * z_min - part_pos
              cur%part_p(3) = -cur%part_p(3)
            ELSE IF (bc_particle(c_bd_z_min) == c_bc_thermal) THEN
              ! Always use the triangle particle weighting for simplicity
              cell_x_r = (cur%part_pos(1) - x_grid_min_local) / dx
              cell_x = FLOOR(cell_x_r + 0.5_num)
              cell_frac_x = REAL(cell_x, num) - cell_x_r
              cell_x = cell_x + 1

              cell_y_r = (cur%part_pos(2) - y_grid_min_local) / dy
              cell_y = FLOOR(cell_y_r + 0.5_num)
              cell_frac_y = REAL(cell_y, num) - cell_y_r
              cell_y = cell_y + 1

              cf2 = cell_frac_x**2
              gx(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_x)
              gx( 0) = 0.75_num - cf2
              gx( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_x)

              cf2 = cell_frac_y**2
              gy(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_y)
              gy( 0) = 0.75_num - cf2
              gy( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_y)

              DO i = 1, 3
                temp(i) = 0.0_num
                DO iy = -1, 1
                  DO ix = -1, 1
                    temp(i) = temp(i) + gx(ix) * gy(iy) &
                        * species%ext_temp_z_min(cell_x+ix, cell_y+iy, i)
                  ENDDO
                ENDDO
              ENDDO

              ! x-direction
              i = 1
              cur%part_p(i) = momentum_from_temperature(&
                  species%mass, temp(i), 0.0_num)

              ! y-direction
              i = 2
              cur%part_p(i) = momentum_from_temperature(&
                  species%mass, temp(i), 0.0_num)

              ! z-direction
              i = 3
              cur%part_p(i) = ABS(momentum_from_temperature(&
                  species%mass, temp(i), 0.0_num))

              cur%part_pos(3) = 2.0_num * z_min - part_pos

            ELSE IF (bc_particle(c_bd_z_min) == c_bc_periodic) THEN
              zbd = -1
              cur%part_pos(3) = part_pos + length_z
            ELSE
              ! Default to open boundary conditions - remove particle
              out_of_bounds = .TRUE.
            ENDIF
          ENDIF
        ENDIF

        ! Particle has left this processor
        IF (part_pos >= z_max_local) THEN
          zbd = 1
          ! Particle has left the system
          IF (z_max_boundary) THEN
            zbd = 0
            IF (bc_particle(c_bd_z_max) == c_bc_reflect) THEN
              cur%part_pos(3) = 2.0_num * z_max - part_pos
              cur%part_p(3) = -cur%part_p(3)
            ELSE IF (bc_particle(c_bd_z_max) == c_bc_thermal) THEN
              ! Always use the triangle particle weighting for simplicity
              cell_x_r = (cur%part_pos(1) - x_grid_min_local) / dx
              cell_x = FLOOR(cell_x_r + 0.5_num)
              cell_frac_x = REAL(cell_x, num) - cell_x_r
              cell_x = cell_x + 1

              cell_y_r = (cur%part_pos(2) - y_grid_min_local) / dy
              cell_y = FLOOR(cell_y_r + 0.5_num)
              cell_frac_y = REAL(cell_y, num) - cell_y_r
              cell_y = cell_y + 1

              cf2 = cell_frac_x**2
              gx(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_x)
              gx( 0) = 0.75_num - cf2
              gx( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_x)

              cf2 = cell_frac_y**2
              gy(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_y)
              gy( 0) = 0.75_num - cf2
              gy( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_y)

              DO i = 1, 3
                temp(i) = 0.0_num
                DO iy = -1, 1
                  DO ix = -1, 1
                    temp(i) = temp(i) + gx(ix) * gy(iy) &
                        * species%ext_temp_z_max(cell_x+ix, cell_y+iy, i)
                  ENDDO
                ENDDO
              ENDDO

              ! x-direction
              i = 1
              cur%part_p(i) = momentum_from_temperature(&
                  species%mass, temp(i), 0.0_num)

              ! y-direction
              i = 2
              cur%part_p(i) = momentum_from_temperature(&
                  species%mass, temp(i), 0.0_num)

              ! z-direction
              i = 3
              cur%part_p(i) = -ABS(momentum_from_temperature(&
                  species%mass, temp(i), 0.0_num))

              cur%part_pos(3) = 2.0_num * z_max - part_pos

            ELSE IF (bc_particle(c_bd_z_max) == c_bc_periodic) THEN
              zbd = 1
              cur%part_pos(3) = part_pos - length_z
            ELSE
              ! Default to open boundary conditions - remove particle
              out_of_bounds = .TRUE.
            ENDIF
          ENDIF
        ENDIF

        IF (out_of_bounds) THEN
          ! Particle has gone forever
          CALL remove_particle_from_partlist(species%attached_list, cur)
          DEALLOCATE(cur)
        ELSE IF (ABS(xbd) + ABS(ybd) + ABS(zbd) > 0) THEN
          ! Particle has left processor, send it to its neighbour
          CALL remove_particle_from_partlist(species%attached_list, cur)
          CALL add_particle_to_partlist(send(xbd, ybd, zbd), cur)
        ENDIF

        ! Move to next particle
        cur => next
      ENDDO

      ! swap Particles
      DO iz = -1, 1
        DO iy = -1, 1
          DO ix = -1, 1
            IF (ABS(ix) + ABS(iy) + ABS(iz) == 0) CYCLE
            ixp = -ix
            iyp = -iy
            izp = -iz
            CALL partlist_sendrecv(send(ix, iy, iz), recv(ixp, iyp, izp), &
                neighbour(ix, iy, iz), neighbour(ixp, iyp, izp))
            CALL append_partlist(species%attached_list, recv(ixp, iyp, izp))
          ENDDO
        ENDDO
      ENDDO

      DO iz = -1, 1
        DO iy = -1, 1
          DO ix = -1, 1
            IF (ABS(ix) + ABS(iy) + ABS(iz) == 0) CYCLE
            CALL destroy_partlist(send(ix, iy, iz))
            CALL destroy_partlist(recv(ix, iy, iz))
          ENDDO
        ENDDO
      ENDDO

    ENDDO

  END SUBROUTINE particle_bcs



  SUBROUTINE current_bcs

    INTEGER :: i

    ! domain is decomposed. Just add currents at edges
    CALL processor_summation_bcs(jx, jng, c_dir_x)
    CALL processor_summation_bcs(jy, jng, c_dir_y)
    CALL processor_summation_bcs(jz, jng, c_dir_z)

    DO i = 1, 2*c_ndims
      IF (bc_particle(i) == c_bc_reflect) THEN
        CALL field_clamp_zero(jx, jng, c_stagger_jx, i)
        CALL field_clamp_zero(jy, jng, c_stagger_jy, i)
        CALL field_clamp_zero(jz, jng, c_stagger_jz, i)
      ENDIF
    ENDDO

  END SUBROUTINE current_bcs

END MODULE boundary
