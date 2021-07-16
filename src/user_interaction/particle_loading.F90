MODULE particle_loading

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

END MODULE particle_loading
