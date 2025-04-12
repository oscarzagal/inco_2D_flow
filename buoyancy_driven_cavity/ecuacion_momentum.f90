module ecuacion_momentum
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use variables_globales
  use esquema_upwind, only:esquema_upwind_
  use gauss_seidel, only:gauss_seidel_
  use gradiente_explicito
  implicit none

  public :: ecuacion_momentum_
  private

  contains

  ! Subrutina para el calculo de las fuerzas de cuerpo (fuerza de flotacion)
  subroutine fuerzas_cuerpo_()

    do j=2,ny-1
       do i=2,nx-1
        fx_flotacion(i,j)=g*sin(theta)*(beta*(T(i,j)-T_inf))*vol(i,j)
        fy_flotacion(i,j)=g*cos(theta)*(beta*(T(i,j)-T_inf))*vol(i,j)
       end do
    end do

  end subroutine fuerzas_cuerpo_

  subroutine ecuacion_momentum_()

    ! Calculo del gradiente de presion
    call gradiente_explicito_()

    ! Calculo de la fuerza de flotacion
    call fuerzas_cuerpo_()

    call esquema_upwind_(ap_u,ae_u,aw_u,an_u,as_u,b_u,gPstar_u_vol &
         ,fx_flotacion,u_star,lambdau)
    call esquema_upwind_(ap_v,ae_v,aw_v,an_v,as_v,b_v,gPstar_v_vol &
         ,fy_flotacion,v_star,lambdav)


    ! Resolucion de la ecuacion de momentum en "u"
    call gauss_seidel_(ap_u,ae_u,aw_u,an_u,as_u,b_u,u_star)

    ! Resolucion de la ecuacion de momentum en "v"
    call gauss_seidel_(ap_v,ae_v,aw_v,an_v,as_v,b_v,v_star)


  end subroutine ecuacion_momentum_

end module ecuacion_momentum
