PROGRAM particle_load_test

  USE particle_loading
  USE random_generator
  USE shared_data

  IMPLICIT NONE

  TYPE(particle), POINTER :: current
  INTEGER, PARAMETER :: nsamples = 100, fu = 70
  REAL(num), PARAMETER :: x0 = 0.0_num, y0 = 2.0_num, z0 = 4.0_num
  INTEGER :: i

  ! Set dx, dy, dz
  dx = 1.0_num
  dy = 2.0_num
  dz = 3.0_num

  ALLOCATE(current)

  ! Dump data to file to check in python
  OPEN(unit=fu, status='REPLACE', file='particle_pos.dat', iostat=errcode)

  DO i = 1, nsamples
    CALL init_particle_position(current, x0, y0, z0, dx, dy, dz)
    WRITE(fu,*) current%part_pos
  END DO

  CLOSE(fu)

  DEALLOCATE(current)

END PROGRAM particle_load_test
