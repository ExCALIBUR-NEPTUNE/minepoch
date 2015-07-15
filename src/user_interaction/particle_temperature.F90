MODULE particle_temperature

  USE shared_data
  USE random_generator

  IMPLICIT NONE

CONTAINS

  ! Subroutine to initialise a thermal particle distribution
  ! Assumes linear interpolation of temperature between cells
  SUBROUTINE setup_particle_temperature(temperature, direction, part_species, &
      drift)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(IN) :: temperature
    INTEGER, INTENT(IN) :: direction
    TYPE(particle_species), POINTER :: part_species
    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(IN) :: drift
    TYPE(particle_list), POINTER :: partlist
    REAL(num) :: mass, temp_local, drift_local
    TYPE(particle), POINTER :: current
    INTEGER(i8) :: ipart
    INTEGER :: ix, iy, iz
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, gy, gz
    REAL(num) :: cell_x_r, cell_frac_x
    REAL(num) :: cell_y_r, cell_frac_y
    REAL(num) :: cell_z_r, cell_frac_z
    INTEGER :: cell_x, cell_y, cell_z
    REAL(num) :: cx2, cy2, cz2
    REAL(num), PARAMETER :: third = 1.0_num / 3.0_num
    REAL(num), PARAMETER :: fac1 = 0.125_num * third
    REAL(num), PARAMETER :: fac2 = 0.5_num * third
    REAL(num), PARAMETER :: fac3 = 7.1875_num * third

    partlist => part_species%attached_list
    current => partlist%head
    ipart = 0
    DO WHILE(ipart < partlist%count)
      mass = current%mass

      ! Assume that temperature is cell centred
      cell_x_r = (current%part_pos(1) - x_grid_min_local) / dx
      cell_y_r = (current%part_pos(2) - y_grid_min_local) / dy
      cell_z_r = (current%part_pos(3) - z_grid_min_local) / dz
      cell_x = FLOOR(cell_x_r + 0.5_num)
      cell_y = FLOOR(cell_y_r + 0.5_num)
      cell_z = FLOOR(cell_z_r + 0.5_num)
      cell_frac_x = REAL(cell_x, num) - cell_x_r
      cell_frac_y = REAL(cell_y, num) - cell_y_r
      cell_frac_z = REAL(cell_z, num) - cell_z_r
      cell_x = cell_x + 1
      cell_y = cell_y + 1
      cell_z = cell_z + 1

      cx2 = cell_frac_x**2
      gx(-2) = fac1 * (0.5_num + cell_frac_x)**4
      gx(-1) = fac2 * (1.1875_num + 2.75_num * cell_frac_x &
          + cx2 * (1.5_num - cell_frac_x - cx2))
      gx( 0) = 0.25_num * (fac3 + cx2 * (cx2 - 2.5_num))
      gx( 1) = fac2 * (1.1875_num - 2.75_num * cell_frac_x &
          + cx2 * (1.5_num + cell_frac_x - cx2))
      gx( 2) = fac1 * (0.5_num - cell_frac_x)**4

      cy2 = cell_frac_y**2
      gy(-2) = fac1 * (0.5_num + cell_frac_y)**4
      gy(-1) = fac2 * (1.1875_num + 2.75_num * cell_frac_y &
          + cy2 * (1.5_num - cell_frac_y - cy2))
      gy( 0) = 0.25_num * (fac3 + cy2 * (cy2 - 2.5_num))
      gy( 1) = fac2 * (1.1875_num - 2.75_num * cell_frac_y &
          + cy2 * (1.5_num + cell_frac_y - cy2))
      gy( 2) = fac1 * (0.5_num - cell_frac_y)**4

      cz2 = cell_frac_z**2
      gz(-2) = fac1 * (0.5_num + cell_frac_z)**4
      gz(-1) = fac2 * (1.1875_num + 2.75_num * cell_frac_z &
          + cz2 * (1.5_num - cell_frac_z - cz2))
      gz( 0) = 0.25_num * (fac3 + cz2 * (cz2 - 2.5_num))
      gz( 1) = fac2 * (1.1875_num - 2.75_num * cell_frac_z &
          + cz2 * (1.5_num + cell_frac_z - cz2))
      gz( 2) = fac1 * (0.5_num - cell_frac_z)**4

      temp_local = 0.0_num
      drift_local = 0.0_num
      DO iz = sf_min, sf_max
        DO iy = sf_min, sf_max
          DO ix = sf_min, sf_max
            temp_local = temp_local + gx(ix) * gy(iy) * gz(iz) &
                * temperature(cell_x+ix, cell_y+iy, cell_z+iz)
            drift_local = drift_local + gx(ix) * gy(iy) * gz(iz) &
                * drift(cell_x+ix, cell_y+iy, cell_z+iz)
          ENDDO
        ENDDO
      ENDDO

      IF (direction == c_dir_x) current%part_p(1) = &
          momentum_from_temperature(mass, temp_local, drift_local)

      IF (direction == c_dir_y) current%part_p(2) = &
          momentum_from_temperature(mass, temp_local, drift_local)

      IF (direction == c_dir_z) current%part_p(3) = &
          momentum_from_temperature(mass, temp_local, drift_local)

      current => current%next
      ipart = ipart + 1
    ENDDO

  END SUBROUTINE setup_particle_temperature



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

END MODULE particle_temperature
