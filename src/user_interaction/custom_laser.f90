MODULE custom_laser

  USE shared_data

  IMPLICIT NONE

CONTAINS

  FUNCTION custom_laser_time_profile(laser)

    TYPE(laser_block), INTENT(IN) :: laser
    REAL(num) :: custom_laser_time_profile

    custom_laser_time_profile = 1.0_num

  END FUNCTION custom_laser_time_profile



  SUBROUTINE custom_laser_update_phase(laser)

    TYPE(laser_block), POINTER :: laser
    !INTEGER :: i, j, k
    !REAL(num) :: xx, yy, zz

    ! use_phase_function only needs setting to .TRUE. if the phase is
    ! time-varying
    !laser%use_phase_function = .TRUE.
    !SELECT CASE(laser%boundary)
    !  CASE(c_bd_x_min, c_bd_x_max)
    !    DO k = 1,nz
    !    zz = z(k)
    !    DO j = 1,ny
    !      yy = y(j)
    !      laser%phase(j,k) = phase_function_x(yy,zz,time)
    !    ENDDO
    !    ENDDO
    !  CASE(c_bd_y_min, c_bd_y_max)
    !    DO k = 1,nz
    !    zz = z(k)
    !    DO i = 1,nx
    !      xx = x(i)
    !      laser%phase(i,k) = phase_function_y(xx,zz,time)
    !    ENDDO
    !    ENDDO
    !  CASE(c_bd_z_min, c_bd_z_max)
    !    DO j = 1,ny
    !    yy = y(j)
    !    DO i = 1,nx
    !      xx = x(i)
    !      laser%phase(i,j) = phase_function_z(xx,yy,time)
    !    ENDDO
    !    ENDDO
    !END SELECT

  END SUBROUTINE custom_laser_update_phase



  SUBROUTINE custom_laser_update_profile(laser)

    TYPE(laser_block), POINTER :: laser
    !INTEGER :: i, j, k
    !REAL(num) :: xx, yy, zz

    ! use_profile_function only needs setting to .TRUE. if the profile is
    ! time-varying
    !laser%use_profile_function = .TRUE.
    !SELECT CASE(laser%boundary)
    !  CASE(c_bd_x_min, c_bd_x_max)
    !    DO k = 1,nz
    !    zz = z(k)
    !    DO j = 1,ny
    !      yy = y(j)
    !      laser%profile(j,k) = profile_function_x(yy,zz,time)
    !    ENDDO
    !    ENDDO
    !  CASE(c_bd_y_min, c_bd_y_max)
    !    DO k = 1,nz
    !    zz = z(k)
    !    DO i = 1,nx
    !      xx = x(i)
    !      laser%profile(i,k) = profile_function_y(xx,zz,time)
    !    ENDDO
    !    ENDDO
    !  CASE(c_bd_z_min, c_bd_z_max)
    !    DO j = 1,ny
    !    yy = y(j)
    !    DO i = 1,nx
    !      xx = x(i)
    !      laser%profile(i,j) = profile_function_z(xx,yy,time)
    !    ENDDO
    !    ENDDO
    !END SELECT

  END SUBROUTINE custom_laser_update_profile

END MODULE custom_laser
