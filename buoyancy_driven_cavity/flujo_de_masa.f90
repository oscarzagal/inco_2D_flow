module flujo_de_masa
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use variables_globales
  use utilidades
  use gradiente_explicito
  use malla
  implicit none

  ! Variables de apoyo
  real(dp) :: gPstar_e_i,gPstar_w_i,gPstar_n_i,gPstar_s_i
  real(dp) :: de_i,dw_i,dn_i,ds_i

  public :: actualizar_flujo_de_masa
  private

  contains

  ! Esta subrutina calcula el coeficiente "d"
  subroutine coef_d()

    do j=2,ny-1
       do i=2,nx-1
          dP_u(i,j)=vol(i,j)/ap_u(i,j)
          dE_u(i,j)=vol(i+1,j)/ap_u(i+1,j)
          dW_u(i,j)=vol(i-1,j)/ap_u(i-1,j)
          dP_v(i,j)=vol(i,j)/ap_v(i,j)
          dN_v(i,j)=vol(i,j+1)/ap_v(i,j+1)
          dS_v(i,j)=vol(i,j-1)/ap_v(i,j-1)
       end do
    end do

  end subroutine coef_d

  ! TODO: añadir el termino para las fuerzas de cuerpo

  subroutine actualizar_flujo_de_masa()

    call coef_d()

    do j=2,ny-1
      do i=2,nx-1

        ! Frontera este de Rhie-Chow
        if (i==nx-1) then
         ue(i,j)=0.0_dp
        else

         ! Velocidad interpolada
         ue_i(i,j)=ge(i,j)*u_star(i+1,j)+(1.0_dp-ge(i,j))*u_star(i,j)

         ! Coeficiente "d" interpolado
         de_i=ge(i,j)*dE_u(i,j)+(1.0_dp-ge(i,j))*dP_u(i,j)

         ! Gradiente de presion en la cara interpolado
         gPstar_e_i=ge(i,j)*gPstar_u_vol(i+1,j)/vol(i+1,j)+(1.0_dp-ge(i,j))&
         *gPstar_u_vol(i,j)/vol(i,j)

         ue(i,j)=ue_i(i,j)-de_i*((Pstar(i+1,j)-Pstar(i,j))/(x(i+1)-x(i))&
              -gPstar_e_i)+(1.0_dp-lambdau)*(ue_n(i,j)-ue_i_n(i,j))

        end if

        ! Frontera oeste de Rhie-Chow
        if (i==2) then
          uw(i,j)=0.0_dp
        else

           ! Velocidad interpolada
           uw_i(i,j)=gw(i,j)*u_star(i-1,j)+(1.0_dp-gw(i,j))*u_star(i,j)

           ! Coeficiente "d" interpolado
           dw_i=gw(i,j)*dW_u(i,j)+(1.0_dp-gw(i,j))*dP_u(i,j)

           ! Gradiente de presion en la cara interpolado
           gPstar_w_i=gw(i,j)*gPstar_u_vol(i-1,j)/vol(i-1,j)+(1.0_dp-gw(i,j))&
           *gPstar_u_vol(i,j)/vol(i,j)

           uw(i,j)=uw_i(i,j)-dw_i*((Pstar(i,j)-Pstar(i-1,j))/(x(i)-x(i-1))&
                -gPstar_w_i)+(1.0_dp-lambdau)*(uw_n(i,j)-uw_i_n(i,j))

        end if

        ! Frontera norte de Rhie-Chow
        if (j==ny-1) then
          vn(i,j)=0.0_dp
        else

           ! Velocidad interpolada
           vn_i(i,j)=gn(i,j)*v_star(i,j+1)+(1.0_dp-gn(i,j))*v_star(i,j)

           ! Coeficiente "d" interpolado
           dn_i=gn(i,j)*dN_v(i,j)+(1.0_dp-gn(i,j))*dP_v(i,j)

           ! Gradiente de presion en la cara interpolado
           gPstar_n_i=gn(i,j)*gPstar_v_vol(i,j+1)/vol(i,j+1)+(1.0_dp-gn(i,j))&
           *gPstar_v_vol(i,j)/vol(i,j)

           vn(i,j)=vn_i(i,j)-dn_i*((Pstar(i,j+1)-Pstar(i,j))/(y(j+1)-y(j))&
                -gPstar_n_i)+(1.0_dp-lambdav)*(vn_n(i,j)-vn_i_n(i,j))

        end if

        ! Frontera sur de Rhie-Chow
        if (j==2) then
          vs(i,j)=0.0_dp
        else

           ! Velocidad interpolada
           vs_i(i,j)=gs(i,j)*v_star(i,j-1)+(1.0_dp-gs(i,j))*v_star(i,j)

           ! Coeficiente "d" interpolado
           ds_i=gs(i,j)*dS_v(i,j)+(1.0_dp-gs(i,j))*dP_v(i,j)

           ! Gradiente de presion en la cara interpolado
           gPstar_s_i=gs(i,j)*gPstar_v_vol(i,j-1)/vol(i,j-1)+(1.0_dp-gs(i,j))&
           *gPstar_v_vol(i,j)/vol(i,j)

           vs(i,j)=vs_i(i,j)-ds_i*((Pstar(i,j)-Pstar(i,j-1))/(y(j)-y(j-1))&
                -gPstar_s_i)+(1.0_dp-lambdav)*(vs_n(i,j)-vs_i_n(i,j))

        end if

        me_star(i,j)=ue(i,j)*deltay(j)
        mw_star(i,j)=-uw(i,j)*deltay(j)
        mn_star(i,j)=vn(i,j)*deltax(i)
        ms_star(i,j)=-vs(i,j)*deltax(i)

      end do
    end do

  end subroutine actualizar_flujo_de_masa

end module flujo_de_masa
