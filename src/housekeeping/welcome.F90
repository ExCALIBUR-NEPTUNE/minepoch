!******************************************************************************
! Welcome message routines
!******************************************************************************

MODULE welcome

  USE version_data
  USE shared_data

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: welcome_message, create_ascii_header

CONTAINS

  !****************************************************************************
  ! This routine prints the welcome message, MPI status
  !****************************************************************************

  SUBROUTINE welcome_message

    IF (rank /= 0) RETURN

    WRITE(*,*)
    WRITE(*,*)
    WRITE(*,*)
    WRITE(*,*)

    WRITE(*,'(A)') '        d########P  d########b        .######b          ' &
        // 'd#######  d##P      d##P'
    WRITE(*,'(A)') '       d########P  d###########    d###########     .###' &
        // '#######  d##P      d##P '
    WRITE(*,'(A)') '      ----        ----     ----  -----     ----   ----- ' &
        // '        ----      -- P  '
    WRITE(*,'(A)') '     d########P  d####,,,####P ####.      .#### d###P   ' &
        // '       d############P   '
    WRITE(*,'(A)') '    d########P  d#########P   ####       .###P ####.    ' &
        // '      d############P    '
    WRITE(*,'(A)') '   d##P        d##P           ####     d####   ####.    ' &
        // '     d##P      d##P     '
    WRITE(*,'(A)') '  d########P  d##P            ###########P     #########' &
        // '#P  d##P      d##P      '
    WRITE(*,'(A)') ' d########P  d##P              d######P          #######' &
        // 'P  d##P      d##P       '
    WRITE(*,*)

    CALL create_ascii_header
    CALL compiler_directives

    WRITE(*,*)
    WRITE(*,*) 'Welcome to ', TRIM(c_code_name), ' version ', &
        TRIM(version_string) // '   (commit ' // TRIM(c_commit_id) // ')'
    WRITE(*,*)

    CALL mpi_status_message

  END SUBROUTINE welcome_message



  SUBROUTINE compiler_directives

    WRITE(*,*) 'The code was compiled with the following compile time options'
    WRITE(*,*) '*************************************************************'
#ifdef PARSER_DEBUG
    WRITE(*,*) 'Particle Debug information -DPARSER_DEBUG'
#endif
#ifdef PARTICLE_COUNT_UPDATE
    WRITE(*,*) 'Global particle counting -DPARTICLE_COUNT_UPDATE'
#endif
#ifdef MPI_DEBUG
    WRITE(*,*) 'MPI error handling -DMPI_DEBUG'
#endif
#ifdef NO_IO
    ! There is no need to add a c_def for this since no I/O occurs.
    WRITE(*,*) 'Perform no I/O -DNO_IO'
#endif
    WRITE(*,*) '*************************************************************'

  END SUBROUTINE compiler_directives


  !****************************************************************************
  ! This routine prints the mpi status information
  !****************************************************************************

  SUBROUTINE mpi_status_message

    CHARACTER(LEN=8) :: string

    CALL integer_as_string(nproc, string)

    WRITE(*,*) 'Code is running on ', TRIM(string), ' processing elements'
    WRITE(*,*)
    WRITE(*,*) 'particles_uniformly_distributed = ', particles_uniformly_distributed
    WRITE(*,*)

  END SUBROUTINE mpi_status_message



  SUBROUTINE create_ascii_header

    CHARACTER(LEN=4) :: str
    INTEGER :: i, strmin, strmax, strlen

    ! Parse commit string to get version number
    ! Commit ID begins with the string v[0-9].[0-9].[0-9]-
    strlen = LEN_TRIM(c_commit_id)
    strmin = 2
    strmax = strmin + 4

    ! Version
    DO i = strmin, MIN(strmax,strlen)
      IF (c_commit_id(i:i) == '.') THEN
        str = c_commit_id(strmin:i-1)
        READ(str, '(i9)') c_version
        strmin = i + 1
        strmax = strmin + 4
        EXIT
      ENDIF
    ENDDO

    ! Revision
    DO i = strmin, MIN(strmax,strlen)
      IF (c_commit_id(i:i) == '.') THEN
        str = c_commit_id(strmin:i-1)
        READ(str, '(i9)') c_revision
        strmin = i + 1
        strmax = strmin + 4
        EXIT
      ENDIF
    ENDDO

    ! Minor revision
    DO i = strmin, MIN(strmax,strlen)
      IF (c_commit_id(i:i) == '-') THEN
        str = c_commit_id(strmin:i-1)
        READ(str, '(i9)') c_minor_rev
        strmax = i - 1
        EXIT
      ENDIF
    ENDDO

    version_string = c_commit_id(2:strmax)

    ascii_header = c_code_name // ' v' // TRIM(version_string) // '   ' &
        // c_commit_id

  END SUBROUTINE create_ascii_header



  SUBROUTINE integer_as_string(int_in, string)

    INTEGER, INTENT(IN) :: int_in
    CHARACTER(LEN=*), INTENT(OUT) :: string

    INTEGER :: n_nums
    CHARACTER(LEN=9) :: numfmt

    IF (int_in == 0) THEN
      n_nums = 1
    ELSE
      n_nums = 1 + INT(LOG10(REAL(ABS(int_in), num)))
    ENDIF
    WRITE(numfmt, '(''(I'', I6.6, '')'')') n_nums
    WRITE(string, numfmt) int_in

  END SUBROUTINE integer_as_string

END MODULE welcome
