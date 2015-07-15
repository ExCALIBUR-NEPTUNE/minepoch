MODULE diagnostics

  USE calc_df
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

#ifdef NO_IO
    RETURN
#endif

    timer_walltime = -1.0_num
    IF (step /= last_step) THEN
      last_step = step
      IF (rank == 0 .AND. stdout_frequency > 0 &
          .AND. MOD(step, stdout_frequency) == 0) THEN
        timer_walltime = MPI_WTIME()
        elapsed_time = timer_walltime - walltime_start
        CALL create_timestring(elapsed_time, timestring)
        WRITE(*, '(''Time'', g20.12, '' and iteration'', i12, '' after'', a)') &
            time, step, timestring
      ENDIF
    ENDIF

    CALL time_history(step)

    IF (timer_collect) CALL timer_stop(c_timer_io)

  END SUBROUTINE output_routines



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
