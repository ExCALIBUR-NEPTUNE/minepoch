MODULE particle_init

  USE shared_data
  USE random_generator

  IMPLICIT NONE

CONTAINS

  SUBROUTINE init_particle_position(part, x, y, z, dx, dy, dz)

    TYPE(particle), POINTER, INTENT(INOUT) :: part
    REAL(num), INTENT(IN) :: x, y, z, dx, dy, dz

    part%part_pos(1) = x + (random() - 0.5_num) * dx
    part%part_pos(2) = y + (random() - 0.5_num) * dy
    part%part_pos(3) = z + (random() - 0.5_num) * dz

  END SUBROUTINE init_particle_position


  FUNCTION momentum_from_temperature(mass, temperature, drift)

    REAL(num), INTENT(IN) :: mass, temperature, drift
    REAL(num) :: momentum_from_temperature

    REAL(num) :: stdev
    REAL(num) :: rand1, rand2, w
    REAL(num), SAVE :: val
    LOGICAL, SAVE :: cached = .FALSE.

        ! This is a basic polar Box-Muller transform
    ! It generates gaussian distributed random numbers
    ! The standard deviation (stdev) is related to temperature

    stdev = SQRT(temperature * kb * mass)

    IF (cached) THEN
      cached = .FALSE.
      momentum_from_temperature = val * stdev + drift
    ELSE
      cached = .TRUE.

      DO
        rand1 = random()
        rand2 = random()

        rand1 = 2.0_num * rand1 - 1.0_num
        rand2 = 2.0_num * rand2 - 1.0_num

        w = rand1**2 + rand2**2

        IF (w > c_tiny .AND. w < 1.0_num) EXIT
      ENDDO

      w = SQRT((-2.0_num * LOG(w)) / w)

      momentum_from_temperature = rand1 * w * stdev + drift
      val = rand2 * w
    ENDIF

  END FUNCTION momentum_from_temperature

END MODULE particle_init
