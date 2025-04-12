module condiciones_frontera
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use variables_globales
  use utilidades
  implicit none

  public :: calcular_velocidades_fronteras,actualizar_presion_fronteras, &
       temperaturas_fronteras
  private

  contains

  subroutine calcular_velocidades_fronteras(phi,cf_vel)
    real(dp), dimension(nx,ny), intent(inout) :: phi
    real(dp) :: cf_vel

    do i=1,nx
      phi(i,ny)=cf_vel ! Frontera norte
      phi(i,1)=0.0_dp ! Frontera sur
    end do

    do j=2,ny-1
      phi(nx,j)=0.0_dp ! Frontera este
      phi(1,j)=0.0_dp ! Frontera oeste
    end do

  end subroutine calcular_velocidades_fronteras

  subroutine actualizar_presion_fronteras(P)
    real(dp), dimension(nx,ny), intent(inout) :: P

    do j=2,ny-1
      ! P(nx,j)=P(nx-1,j) ! Frontera este

      P(nx,j)=P(nx-1,j)+(P(nx,j)-0.5_dp*(P(nx-1,j)+P(nx-2,j)))/deltax(j)&
      *(x(nx)-x(nx-1))

      ! P(1,j)=P(2,j) ! Frontera oeste

      P(1,j)=P(2,j)+(0.5_dp*(P(3,j)+P(2,j))-P(1,j))/deltax(j)*(x(2)-x(1))
    end do

    ! do i=1,nx
    do i=2,nx-1
      ! P(i,ny)=P(i,ny-1) ! Frontera norte

      P(i,ny)=P(i,ny-1)+(P(i,ny)-0.5_dp*(P(i,ny-1)+P(i,ny-2)))/deltay(i)&
      *(y(ny)-y(ny-1))

      ! P(i,1)=P(i,2) ! Frontera sur

      P(i,1)=P(i,2)+(0.5_dp*(P(i,3)+P(i,2))-P(i,1))/deltay(i)*(y(2)-y(1))
    end do

  end subroutine actualizar_presion_fronteras

  subroutine temperaturas_fronteras()

    do i=2,nx-1
      T(i,ny)=T(i,ny-1) ! Frontera norte
      T(i,1)=T(i,2) ! Frontera sur
    end do

    do j=2,ny-1
      T(nx,j)=T_F ! Frontera este
      T(1,j)=T_C ! Frontera oeste
    end do

  end subroutine temperaturas_fronteras

end module condiciones_frontera
