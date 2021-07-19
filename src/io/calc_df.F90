MODULE calc_df

  USE shared_data

  IMPLICIT NONE

CONTAINS

  SUBROUTINE calc_total_energy_sum

    REAL(num) :: particle_energy, field_energy
    REAL(num) :: part_ux, part_uy, part_uz
    REAL(num) :: part_mc, part_w, fac, gamma
    REAL(num) :: sum_out(2), sum_in(2)
    REAL(num), PARAMETER :: c2 = c**2
    INTEGER :: ispecies, i, j, k
    TYPE(particle), POINTER :: current
    TYPE(particle_species), POINTER :: species, next_species

    particle_energy = 0.0_num

    ! Sum over all particles to calculate total kinetic energy
    next_species => species_list
    DO ispecies = 1, n_species
      species => next_species
      next_species => species%next

      current => species%attached_list%head
      part_mc = c * species%mass
      part_w = species%weight
      fac = part_mc * part_w * c

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
        part_mc = c * current%mass
        IF (particles_uniformly_distributed) THEN
          part_w = current%weight
        ENDIF
        fac = part_mc * part_w * c

        part_ux = current%part_p(1) / part_mc
        part_uy = current%part_p(2) / part_mc
        part_uz = current%part_p(3) / part_mc
        gamma = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)
        particle_energy = particle_energy + (gamma - 1.0_num) * fac

        current => current%next
      ENDDO
    ENDDO

    ! EM field energy
    field_energy = 0.0_num
    DO k = 1, nz
    DO j = 1, ny
    DO i = 1, nx
      field_energy = field_energy + ex(i,j,k)**2 + ey(i,j,k)**2 &
          + ez(i,j,k)**2 + c2 * (bx(i,j,k)**2 + by(i,j,k)**2 + bz(i,j,k)**2)
    ENDDO
    ENDDO
    ENDDO
    field_energy = 0.5_num * epsilon0 * field_energy * dx * dy * dz

    sum_out(1) = particle_energy
    sum_out(2) = field_energy
    CALL MPI_REDUCE(sum_out, sum_in, 2, mpireal, MPI_SUM, 0, comm, errcode)
    total_particle_energy = sum_in(1)
    total_field_energy = sum_in(2)

  END SUBROUTINE calc_total_energy_sum

END MODULE calc_df
