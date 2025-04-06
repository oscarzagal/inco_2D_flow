module ecuacion_momentum
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use variables_globales
  use esquema_lineal, only:esquema_lineal_
  use esquema_upwind, only:esquema_upwind_
  use esquema_ley_de_potencia, only:esquema_ley_de_potencia_
  use gauss_seidel, only:gauss_seidel_
  use gradiente_explicito
  implicit none

  public :: ecuacion_momentum_
  private

  contains

  ! TODO: esquema Euler implicito (pag 534 PDF Moukalled)

  subroutine ecuacion_momentum_()

    ! Calculo del gradiente de presion
    call gradiente_explicito_()

    select case (esquema)
      case (1)
        call esquema_upwind_(ap_u,ae_u,aw_u,an_u,as_u,b_u,gPstar_u_vol,u_star &
        ,lambdau)
        call esquema_upwind_(ap_v,ae_v,aw_v,an_v,as_v,b_v,gPstar_v_vol,v_star &
        ,lambdav)
      case (2)
        call esquema_lineal_(ap_u,ae_u,aw_u,an_u,as_u,b_u,gPstar_u_vol,u_star &
        ,lambdau)
        call esquema_lineal_(ap_v,ae_v,aw_v,an_v,as_v,b_v,gPstar_v_vol,v_star &
        ,lambdav)
      case (3)
        call esquema_ley_de_potencia_(ap_u,ae_u,aw_u,an_u,as_u,b_u,gPstar_u_vol &
        ,u_star,lambdau)
        call esquema_ley_de_potencia_(ap_v,ae_v,aw_v,an_v,as_v,b_v,gPstar_v_vol &
        ,v_star,lambdav)
      case default
    end select


    ! Resolucion de la ecuacion de momentum en "u"
    call gauss_seidel_(ap_u,ae_u,aw_u,an_u,as_u,b_u,u_star)

    ! Resolucion de la ecuacion de momentum en "v"
    call gauss_seidel_(ap_v,ae_v,aw_v,an_v,as_v,b_v,v_star)


  end subroutine ecuacion_momentum_

end module ecuacion_momentum
