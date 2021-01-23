PROGRAM pic

  ! EPOCH3D is a Birdsall and Langdon type PIC code derived from the PSC
  ! written by Hartmut Ruhl.

  ! The particle pusher (particles.F90) and the field solver (fields.f90) are
  ! almost exact copies of the equivalent routines from PSC, modified slightly
  ! to allow interaction with the changed portions of the code and for
  ! readability. The MPI routines are exactly equivalent to those in PSC, but
  ! are completely rewritten in a form which is easier to extend with arbitrary
  ! fields and particle properties. The support code is entirely new and is not
  ! equivalent to PSC.

  ! EPOCH3D written by C.S.Brady, Centre for Fusion, Space and Astrophysics,
  ! University of Warwick, UK
  ! PSC written by Hartmut Ruhl

  USE balance
  USE deck
  USE diagnostics
  USE fields
  USE helper
  USE ic_module
  USE mpi_routines
  USE particles
  USE setup
  USE shared_data
  USE problem_setup_module
  USE finish
  USE welcome
#ifdef PAT_DEBUG
  USE pat_mpi_lib
#endif
  IMPLICIT NONE

  INTEGER :: ispecies, ierr
  LOGICAL :: push = .TRUE.
  CHARACTER(LEN=*), PARAMETER :: data_dir_file = 'USE_DATA_DIRECTORY'
  CHARACTER(LEN=64) :: timestring
#ifdef PAT_DEBUG
  CHARACTER(LEN=17) :: patc_out_fn = "patc_epoch3d.out"//CHAR(0)
#endif

  REAL(num) :: runtime
  TYPE(particle_species), POINTER :: species, next_species

  step = 0
  time = 0.0_num

  CALL mpi_minimal_init ! mpi_routines.f90
  real_walltime_start = MPI_WTIME()
  CALL minimal_init     ! setup.f90
  CALL timer_init
  CALL setup_partlists  ! partlist.f90
  CALL welcome_message  ! welcome.f90

  IF (rank == 0) THEN
    OPEN(unit=lu, status='OLD', file=TRIM(data_dir_file), iostat=ierr)
    IF (ierr == 0) THEN
      READ(lu,'(A)') data_dir
      CLOSE(lu)
      PRINT*, 'Using data directory "' // TRIM(data_dir) // '"'
    ELSE
      data_dir = 'Data'
    ENDIF
  ENDIF

  CALL MPI_BCAST(data_dir, 64, MPI_CHARACTER, 0, comm, errcode)

  ! Read deck here to tell which set-up
  CALL read_deck
  CALL problem_setup(c_ds_first, problem)
  ! Now re-read to override default values
  CALL read_deck

  CALL setup_particle_boundaries ! boundary.f90
  CALL mpi_initialise  ! mpi_routines.f90
  CALL after_control   ! setup.f90

  CALL problem_setup(c_ds_last, problem)
  CALL after_deck_last

  ! auto_load particles
  CALL auto_load
  time = 0.0_num

  CALL manual_load
  CALL set_dt

  ! Use initial conditions to setup fluid equations.
  CALL setup_fluid

  npart_global = 0

  next_species => species_list
  DO ispecies = 1, n_species
    species => next_species
    next_species => species%next

    npart_global = npart_global + species%count
  ENDDO

  ! .TRUE. to over_ride balance fraction check
  IF (npart_global > 0) CALL balance_workload(.TRUE.)

  IF (explicit_pic) THEN
    CALL particle_bcs
    CALL efield_bcs(ex, ey, ez, ng)
    CALL bfield_final_bcs(bx, by, bz, ng)
  ELSE
    CALL init_trilinos(local_elements, linear_index, comm)
  END IF

  IF (rank == 0) PRINT *, 'Equilibrium set up OK, running code'

  walltime_start = MPI_WTIME()
  CALL output_routines(step) ! diagnostics.f90

  IF (timer_collect) CALL timer_start(c_timer_step)

#ifdef PAT_DEBUG
  CALL pat_mpi_open(patc_out_fn)
#endif

  ! Set-up particle push module
  CALL setup_particle_push

  DO
    IF (timer_collect) THEN
      CALL timer_stop(c_timer_step)
      CALL timer_reset
      timer_first(c_timer_step) = timer_walltime
    ENDIF
    push = (time >= particle_push_start_time)

    IF (explicit_pic) THEN
      CALL update_eb_fields_half

      IF (push) THEN
        ! .FALSE. this time to use load balancing threshold
        IF (use_balance) CALL balance_workload(.FALSE.)
        CALL push_particles
      ENDIF

#ifdef PAT_DEBUG
      CALL pat_mpi_monitor(step,1)
#endif
      step = step + 1
      time = time + dt
 
      CALL update_eb_fields_final
      ! At this point, do the second substep of the push if there are any drift-kinetic particles
      IF (drift_kinetic_species_exist) THEN
        step=step-1
        time=time-dt
        IF (push) THEN
          CALL push_particles_2ndstep
        END IF
        ! Rewind the fields half a step
        CALL rewind_fields_halfstep
        ! Then update using corrected current.
        step = step + 1
        time = time + dt
        CALL update_eb_fields_final      
      END IF
    ELSE
      IF (rank == 0) THEN
        WRITE(*,*) 'Implicit PIC not implemented yet!'
      END IF
      CALL MPI_ABORT(MPI_COMM_WORLD, c_err_not_implemented, ierr)
    END IF
   
    ! Output any diagnostics
    CALL output_routines(step)
    
    IF ((step >= nsteps .AND. nsteps >= 0) .OR. (time >= t_end)) EXIT

    ! This section ensures that the particle count for the species_list
    ! objects is accurate. This makes some things easier, but increases
    ! communication
#ifdef PARTICLE_COUNT_UPDATE
    next_species => species_list
    DO ispecies = 1, n_species
      species => next_species
      next_species => species%next
      CALL MPI_ALLREDUCE(species%attached_list%count, species%count, 1, &
          MPI_INTEGER8, MPI_SUM, comm, errcode)
      species%count_update_step = step
    ENDDO
#endif

  ENDDO

#ifdef PAT_DEBUG
  CALL pat_mpi_close()
#endif

  IF (rank == 0) runtime = MPI_WTIME() - walltime_start

  CALL output_routines(step)

  IF (rank == 0) THEN
    CALL create_full_timestring(runtime, timestring)
    WRITE(*,*) 'Final runtime of core = ' // TRIM(timestring)
  ENDIF

  IF (.NOT. explicit_pic) THEN
    CALL end_trilinos
  END IF

  CALL finalise

END PROGRAM pic
