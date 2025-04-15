module conductancia_difusiva
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use variables_globales

  implicit none

  public :: conductancia_difusiva_
  private

  contains

  subroutine conductancia_difusiva_(gamma,FluxFe_dif,FluxFw_dif,FluxFn_dif &
    ,FluxFs_dif)

    real(dp), intent(in) :: gamma ! Coeificiente difusivo
    real(dp), dimension(nx,ny), intent(out) :: FluxFe_dif,FluxFw_dif,FluxFn_dif &
    ,FluxFs_dif

    do j=2,ny-1
      do i=2,nx-1
        FluxFe_dif(i,j)=gamma*deltay(j)/(x_r(i+1+nx*(j-1))-x_r(i+nx*(j-1)))
        FluxFw_dif(i,j)=gamma*deltay(j)/(x_r(i+nx*(j-1))-x_r(i-1+nx*(j-1)))
        FluxFn_dif(i,j)=gamma*deltax(i)/(y_r(i+nx*(j))-y_r(i+nx*(j-1)))
        FluxFs_dif(i,j)=gamma*deltax(i)/(y_r(i+nx*(j-1))-y_r(i+nx*(j-2)))
      end do
    end do

  end subroutine conductancia_difusiva_

end module conductancia_difusiva
