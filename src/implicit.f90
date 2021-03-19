MODULE implicit

  USE iso_c_binding
  USE shared_data

  IMPLICIT NONE

CONTAINS

  SUBROUTINE init_1d_index

    INTEGER :: ix, iy, iz, row, i, ivar
    INTEGER, PARAMETER :: nvar = 6

    i = 1

    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          row = index1d(ix,iy,iz)
          DO ivar = 0, nvar-1
            linear_index(i) = row + ivar
            i = i + 1
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE init_1d_index



  PURE FUNCTION index1d(ix, iy, iz)

    INTEGER(C_INT) :: index1d
    INTEGER, INTENT(IN) :: ix, iy, iz

    index1d = (ix-1) + (iy-1) * nx + (iz-1) * nx * ny

  END FUNCTION index1d

END MODULE implicit
