module ecuacion_energia
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use variables_globales
  use gauss_seidel
  use condiciones_frontera

  implicit none

  public :: ecuacion_energia_
  private

  contains

  ! NOTE: subrutina provisional, se debe de refactorizar el archivo
  ! "esquea_upwind.f90"

  ! subroutine esquema_upwind_energia()

  !   do j=2,ny-1
  !      do i=2,nx-1

  !         ! Coeficiente central
  !         FluxCe_conv(i,j)=max(me_star(i,j),0.0_dp)
  !         FluxCw_conv(i,j)=max(mw_star(i,j),0.0_dp)
  !         FluxCn_conv(i,j)=max(mn_star(i,j),0.0_dp)
  !         FluxCs_conv(i,j)=max(ms_star(i,j),0.0_dp)

  !         ! Demas coeficientes
  !         FluxFe_conv(i,j)=-max(-me_star(i,j),0.0_dp)
  !         FluxFw_conv(i,j)=-max(-mw_star(i,j),0.0_dp)
  !         FluxFn_conv(i,j)=-max(-mn_star(i,j),0.0_dp)
  !         FluxFs_conv(i,j)=-max(-ms_star(i,j),0.0_dp)

  !      end do
  !   end do

  ! end subroutine esquema_upwind_energia

  subroutine ecuacion_energia_()

      ! call esquema_upwind_energia()

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
               ! )/lambdaT
               +vol(i,j)/deltat)/lambdaT

          ! b_T(i,j)=0.0_dp
          b_T(i,j)=(1.0_dp-lambdaT)*ap_T(i,j)*T(i,j)+T(i,j)*vol(i,j)/deltat

       end do
    end do

    call gauss_seidel_(ap_T,ae_T,aw_T,an_T,as_T,b_T,T)

    ! Actualizar condiciones de frontera
    call temperaturas_fronteras()

  end subroutine ecuacion_energia_

end module ecuacion_energia
