PROGRAM particle_momentum_test

  USE particle_init
  USE random_generator
  USE shared_data

  IMPLICIT NONE

  TYPE(particle), POINTER :: current
  INTEGER, PARAMETER :: nsamples = 100, fu = 70
  REAL(num), PARAMETER :: ev_k = q0 / kb
  REAL(num), PARAMETER :: tempx = 10.0_num * ev_k
  REAL(num), PARAMETER :: tempy = ev_k
  REAL(num), PARAMETER :: tempz = 200.0_num * ev_k
  REAL(num), PARAMETER :: driftx = 0.01_num * c * m0
  REAL(num), PARAMETER :: drifty = 0.1_num * c * m0
  REAL(num), PARAMETER :: driftz = 0.2_num * c * m0
  INTEGER :: i

  ALLOCATE(current)

  ! Dump data to file to check in python
  OPEN(unit=fu, status='REPLACE', file='particle_mom.dat', iostat=errcode)

  DO i = 1, nsamples
    current%part_p(1) = momentum_from_temperature(m0, tempx, driftx)
    current%part_p(2) = momentum_from_temperature(m0, tempy, drifty)
    current%part_p(3) = momentum_from_temperature(m0, tempz, driftz)
    WRITE(fu,*) current%part_p
  END DO

  CLOSE(fu)

  DEALLOCATE(current)

END PROGRAM particle_momentum_test
