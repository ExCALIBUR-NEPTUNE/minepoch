MODULE finish

  USE partlist
  USE laser

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: finalise

CONTAINS

  SUBROUTINE finalise

    CALL deallocate_memory
    CALL MPI_FINALIZE(errcode)

  END SUBROUTINE finalise



  SUBROUTINE deallocate_memory

    INTEGER :: i, stat
    TYPE(particle_species), POINTER :: species, next_species
    TYPE(particle_species), POINTER :: iolist, next_iolist

    DEALLOCATE(x, x_global, xb_global)
    DEALLOCATE(y, y_global, yb_global)
    DEALLOCATE(z, z_global, zb_global)
    DEALLOCATE(ex, ey, ez, bx, by, bz, jx, jy, jz)

    DEALLOCATE(npart_each_rank)
    DEALLOCATE(x_grid_mins, x_grid_maxs, cell_x_min, cell_x_max)
    DEALLOCATE(y_grid_mins, y_grid_maxs, cell_y_min, cell_y_max)
    DEALLOCATE(z_grid_mins, z_grid_maxs, cell_z_min, cell_z_max)

    DEALLOCATE(ex_x_min, ex_x_max, ey_x_min, ey_x_max, ez_x_min, ez_x_max)
    DEALLOCATE(bx_x_min, bx_x_max, by_x_min, by_x_max, bz_x_min, bz_x_max)
    DEALLOCATE(ex_y_min, ex_y_max, ey_y_min, ey_y_max, ez_y_min, ez_y_max)
    DEALLOCATE(bx_y_min, bx_y_max, by_y_min, by_y_max, bz_y_min, bz_y_max)
    DEALLOCATE(ex_z_min, ex_z_max, ey_z_min, ey_z_max, ez_z_min, ez_z_max)
    DEALLOCATE(bx_z_min, bx_z_max, by_z_min, by_z_max, bz_z_min, bz_z_max)

    next_species => species_list
    next_iolist => io_list_data
    DO i = 1, n_species
      species => next_species
      next_species => species%next

      CALL destroy_partlist(species%attached_list)
      DEALLOCATE(species%ext_temp_x_min, STAT=stat)
      DEALLOCATE(species%ext_temp_x_max, STAT=stat)
      DEALLOCATE(species%ext_temp_y_min, STAT=stat)
      DEALLOCATE(species%ext_temp_y_max, STAT=stat)
      DEALLOCATE(species%ext_temp_z_min, STAT=stat)
      DEALLOCATE(species%ext_temp_z_max, STAT=stat)
      DEALLOCATE(species%density, STAT=stat)
      DEALLOCATE(species%temp, STAT=stat)
      DEALLOCATE(species%drift, STAT=stat)
      DEALLOCATE(species, STAT=stat)

      iolist => next_iolist
      next_iolist => iolist%next
      DEALLOCATE(iolist, STAT=stat)
    ENDDO

    CALL deallocate_lasers

    CALL MPI_COMM_FREE(comm, errcode)

  END SUBROUTINE deallocate_memory

END MODULE finish
