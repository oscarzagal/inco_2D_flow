module ecuacion_energia
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use variables_globales
  use gauss_seidel
  use ADI
  use condiciones_frontera

  implicit none

  public :: ecuacion_energia_
  private

  contains

  subroutine ecuacion_energia_()

    do j=2,ny-1
       do i=2,nx-1

          ae_T(i,j)=-FluxFe_dif_T(i,j)-max(-me_star(i,j),0.0_dp)
          aw_T(i,j)=-FluxFw_dif_T(i,j)-max(-mw_star(i,j),0.0_dp)
          an_T(i,j)=-FluxFn_dif_T(i,j)-max(-mn_star(i,j),0.0_dp)
          as_T(i,j)=-FluxFs_dif_T(i,j)-max(-ms_star(i,j),0.0_dp)

          ap_T(i,j)=(max(me_star(i,j),0.0_dp)+max(mw_star(i,j),0.0_dp) &
               +max(mn_star(i,j),0.0_dp)+max(ms_star(i,j),0.0_dp) &
               +FluxFe_dif_T(i,j)+FluxFw_dif_T(i,j) &
               +FluxFn_dif_T(i,j)+FluxFs_dif_T(i,j) &
               +vol(i,j)/deltat)/lambdaT

          b_T(i,j)=(1.0_dp-lambdaT)*ap_T(i,j)*T(i,j)+T(i,j)*vol(i,j)/deltat

       end do
    end do

    ! call gauss_seidel_(ap_T,ae_T,aw_T,an_T,as_T,b_T,T)
    call ADI_(ap_T,ae_T,aw_T,an_T,as_T,b_T,T)

    ! Actualizar condiciones de frontera
    ! call temperaturas_fronteras()

  end subroutine ecuacion_energia_

end module ecuacion_energia
