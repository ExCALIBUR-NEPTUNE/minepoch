MODULE diagnostics

  USE calc_df
  USE shared_data
  USE strings
  USE timer

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: output_routines, create_full_timestring

CONTAINS

  SUBROUTINE output_routines(step, force_write)   ! step = step index

    INTEGER, INTENT(INOUT) :: step
    LOGICAL, INTENT(IN), OPTIONAL :: force_write
    REAL(num) :: elapsed_time
    CHARACTER(LEN=16) :: timestring
    INTEGER, SAVE :: last_step = -1
    LOGICAL :: fw
#ifdef NO_IO
    RETURN
#endif
    IF(PRESENT(force_write)) THEN
       fw = force_write
    ELSE
       fw=.false.
    END IF

    timer_walltime = -1.0_num
    IF (step /= last_step) THEN
      last_step = step
      IF (rank == 0 .AND. ((stdout_frequency > 0 &
          .AND. MOD(step, stdout_frequency) == 0) .OR. fw)) THEN
        timer_walltime = MPI_WTIME()
        elapsed_time = timer_walltime - walltime_start
        CALL create_timestring(elapsed_time, timestring)
        WRITE(*, '(''Time'', g20.12, '' and iteration'', i12, '' after'', a)') &
            time, step, timestring
      ENDIF
    ENDIF

    IF (write_time_history) CALL time_history(step)

    IF (timer_collect) CALL timer_stop(c_timer_io)

    IF (n_field_probes > 0) CALL write_field_probes(step)

  END SUBROUTINE output_routines



  SUBROUTINE write_field_probes(step)

    INTEGER, INTENT(IN) :: step
    INTEGER :: iprobe
    CHARACTER(LEN=c_max_string_length) :: filename
    CHARACTER(LEN=11+data_dir_max_length) :: full_filename
    LOGICAL :: first_call = .TRUE.
    INTEGER :: errcode
    INTEGER :: ixl, iyl, izl
    INTEGER, PARAMETER :: fu = 68
    LOGICAL, PARAMETER :: do_flush = .TRUE.
    INTEGER, PARAMETER :: dump_frequency = 10

    IF (first_call) THEN
      DO iprobe = 1, n_field_probes
        ! Check if on proc
        IF (.NOT. on_proc(field_probes(1, iprobe), &
            field_probes(2, iprobe), field_probes(3, iprobe))) CYCLE

        ! Write file header
        WRITE(filename, '(''/Field_Probe_'', i2.2, ''.dat'')') iprobe
        full_filename = TRIM(data_dir) // filename
        OPEN(unit=fu, status='REPLACE', file=TRIM(full_filename), &
            iostat=errcode)
        IF (errcode /= 0) THEN
          PRINT*, 'Failed to open file: ', TRIM(full_filename)
          CALL MPI_ABORT(MPI_COMM_WORLD, c_err_io, errcode)
        END IF

        ! Local coordinates of probe
        ixl = field_probes(1, iprobe) - n_global_min(1)
        iyl = field_probes(2, iprobe) - n_global_min(2)
        izl = field_probes(3, iprobe) - n_global_min(3)
        WRITE(fu,'(3a5, 6a23)') 'i', 'j', 'k', &
            'x(ix)', 'y(iy)', 'z(iz)', &
            'dx', 'dy', 'dz'
        WRITE(fu, '(3i5,6e23.14)') field_probes(1:3, iprobe), &
            x(ixl), y(iyl), z(izl), dx, dy, dz
        WRITE(fu,'(a5, 99a23)') 'time', 'ex', 'ey', 'ez', 'bx', 'by', 'bz'
        IF (do_flush) CLOSE(unit=fu)
      END DO
      first_call = .FALSE.
    END IF

    IF (MOD(step, dump_frequency) /= 0) RETURN

    DO iprobe = 1, n_field_probes
      ! Check if on proc
      IF (.NOT. on_proc(field_probes(1, iprobe), &
          field_probes(2, iprobe), field_probes(3, iprobe))) CYCLE
      WRITE(filename, '(''/Field_Probe_'', i2.2, ''.dat'')') iprobe
      full_filename = TRIM(data_dir) // filename
      OPEN(unit=fu, status='OLD', position='APPEND', file=TRIM(full_filename), &
          iostat=errcode)
      IF (errcode /= 0) THEN
        PRINT*, 'Failed to open file: ', TRIM(full_filename)
        CALL MPI_ABORT(MPI_COMM_WORLD, c_err_io, errcode)
      END IF
      ixl = field_probes(1, iprobe) - n_global_min(1)
      iyl = field_probes(2, iprobe) - n_global_min(2)
      izl = field_probes(3, iprobe) - n_global_min(3)
      WRITE(fu,'(99e23.14)') time, ex(ixl,iyl,izl), ey(ixl,iyl,izl), ez(ixl,iyl,izl), &
          bx(ixl,iyl,izl), by(ixl,iyl,izl), bz(ixl,iyl,izl)
      IF (do_flush) CLOSE(unit=fu)
    END DO

  END SUBROUTINE write_field_probes


  PURE LOGICAL FUNCTION on_proc(ix, iy, iz)

    INTEGER, INTENT(IN) :: ix, iy, iz

    on_proc = ix >= n_global_min(1) .AND. ix <= n_global_max(1) &
        .AND. iy >= n_global_min(2) .AND. iy <= n_global_max(2) &
        .AND. iz >= n_global_min(3) .AND. iz <= n_global_max(3)

  END FUNCTION on_proc



  SUBROUTINE create_timestring(time, timestring)

    REAL(num), INTENT(IN) :: time
    CHARACTER(LEN=*), INTENT(INOUT) :: timestring ! length at least 15
    INTEGER :: days, hours, minutes, seconds, frac_seconds

    days = INT(time) / 60 / 60  / 24
    hours = INT(time) / 60 / 60 - days * 24
    minutes = INT(time) / 60 - (days * 24 + hours) * 60
    seconds = INT(time) - ((days * 24 + hours) * 60 + minutes) * 60
    frac_seconds = FLOOR((time - INT(time)) * 100)

    WRITE(timestring, '(i3,'':'',i2.2,'':'',i2.2,'':'',i2.2,''.'',i2.2)') &
        days, hours, minutes, seconds, frac_seconds

  END SUBROUTINE create_timestring



  SUBROUTINE create_full_timestring(time, timestring)

    REAL(num), INTENT(IN) :: time
    CHARACTER(LEN=*), INTENT(INOUT) :: timestring ! length at least 48
    INTEGER :: days, hours, minutes, seconds, frac_seconds, var
    CHARACTER(LEN=8) :: varstring
    CHARACTER(LEN=4) :: intstring, fracstring
    LOGICAL :: string_started

    days = INT(time) / 60 / 60  / 24
    hours = INT(time) / 60 / 60 - days * 24
    minutes = INT(time) / 60 - (days * 24 + hours) * 60
    seconds = INT(time) - ((days * 24 + hours) * 60 + minutes) * 60
    frac_seconds = FLOOR((time - INT(time)) * 100)

    timestring = ''
    string_started = .FALSE.

    var = days
    varstring = ' day'
    IF (var > 0) THEN
      CALL integer_as_string(var, intstring)
      IF (string_started) THEN
        timestring = TRIM(timestring) // ', ' // TRIM(intstring) &
            // TRIM(varstring)
      ELSE
        timestring = TRIM(timestring) // TRIM(intstring) // TRIM(varstring)
      ENDIF
      IF (var > 1) timestring = TRIM(timestring) // 's'
      string_started = .TRUE.
    ENDIF

    var = hours
    varstring = ' hour'
    IF (var > 0) THEN
      CALL integer_as_string(var, intstring)
      IF (string_started) THEN
        timestring = TRIM(timestring) // ', ' // TRIM(intstring) &
            // TRIM(varstring)
      ELSE
        timestring = TRIM(timestring) // TRIM(intstring) // TRIM(varstring)
      ENDIF
      IF (var > 1) timestring = TRIM(timestring) // 's'
      string_started = .TRUE.
    ENDIF

    var = minutes
    varstring = ' minute'
    IF (var > 0) THEN
      CALL integer_as_string(var, intstring)
      IF (string_started) THEN
        timestring = TRIM(timestring) // ', ' // TRIM(intstring) &
            // TRIM(varstring)
      ELSE
        timestring = TRIM(timestring) // TRIM(intstring) // TRIM(varstring)
      ENDIF
      IF (var > 1) timestring = TRIM(timestring) // 's'
      string_started = .TRUE.
    ENDIF

    var = seconds
    varstring = ' seconds'
    IF (var > 0 .OR. frac_seconds > 0 .OR. .NOT.string_started) THEN
      CALL integer_as_string(var, intstring)
      WRITE(fracstring, '(i2.2)') frac_seconds
      IF (string_started) THEN
        timestring = TRIM(timestring) // ', ' // TRIM(intstring) // '.' &
            // TRIM(fracstring) // TRIM(varstring)
      ELSE
        timestring = TRIM(timestring) // TRIM(intstring) // '.' &
            // TRIM(fracstring) // TRIM(varstring)
      ENDIF
    ENDIF

  END SUBROUTINE create_full_timestring



  SUBROUTINE time_history(step)   ! step = step index

    INTEGER, INTENT(INOUT) :: step
    CHARACTER(LEN=11+data_dir_max_length) :: full_filename
    LOGICAL, SAVE :: first_call = .TRUE.
    INTEGER :: errcode
    INTEGER, PARAMETER :: fu = 67
    LOGICAL, PARAMETER :: do_flush = .TRUE.
    INTEGER, PARAMETER :: dump_frequency = 10

    full_filename = TRIM(data_dir) // '/output.dat'

    IF (first_call) THEN
      ! open the file
      IF (rank == 0) THEN
        OPEN(unit=fu, status='REPLACE', file=TRIM(full_filename), &
            iostat=errcode)
        IF (errcode /= 0) THEN
          PRINT*, 'Failed to open file: ', TRIM(full_filename)
          CALL MPI_ABORT(MPI_COMM_WORLD, c_err_io, errcode)
        END IF
        WRITE(fu,'(''# '',a5,99a23)') 'step', 'time', 'dt', &
            'laser_injected', 'laser_absorbed', 'total_particle_energy', &
            'total_field_energy'
        IF (do_flush) CLOSE(unit=fu)
      ENDIF
      first_call = .FALSE.
    ENDIF

    IF (MOD(step, dump_frequency) /= 0) RETURN

    IF (timer_collect) THEN
      IF (timer_walltime < 0.0_num) THEN
        CALL timer_start(c_timer_io)
      ELSE
        CALL timer_start(c_timer_io, .TRUE.)
      ENDIF
    ENDIF

    io_list => species_list

    CALL MPI_ALLREDUCE(laser_absorb_local, laser_absorbed, 1, mpireal, &
        MPI_SUM, comm, errcode)
    CALL MPI_ALLREDUCE(laser_inject_local, laser_injected, 1, mpireal, &
        MPI_SUM, comm, errcode)
    IF (laser_injected > 0.0_num) THEN
      laser_absorbed = laser_absorbed / laser_injected
    ELSE
      laser_absorbed = 0.0_num
    ENDIF

    CALL calc_total_energy_sum

    IF (rank == 0) THEN
      IF (do_flush) THEN
        OPEN(unit=fu, status='OLD', position='APPEND', &
            file=TRIM(full_filename), iostat=errcode)
      ENDIF

      WRITE(fu,'(i7,99e23.14)') step, time, dt, laser_injected, laser_absorbed, &
          total_particle_energy, total_field_energy

      IF (do_flush) CLOSE(unit=fu)
    ENDIF

  END SUBROUTINE time_history

END MODULE diagnostics
