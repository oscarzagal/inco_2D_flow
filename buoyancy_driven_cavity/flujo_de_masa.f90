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

  ! Fuerzas de cuerpo en las caras promediadas 1 vez
  real(dp), dimension(nx,ny) :: B_e_i,B_w_i,B_n_i,B_s_i

  ! Fuerzas de cuerpo en el centroide de la celda promediadas 2 veces
  real(dp), dimension(nx,ny) :: B_u_ii,B_v_ii

  ! Fuerzas de cuerpo en las caras promediadas 3 veces
  real(dp), dimension(nx,ny) :: B_e_iii,B_w_iii,B_n_iii,B_s_iii

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

  subroutine fuerzas_flotacion()

      do j=2,ny-1
        do i=2,nx-1
         ! Promediado 1 vez
         B_e_i(i,j)=ge(i,j)*fx_flotacion(i+1,j)+(1.0_dp-ge(i,j))*fx_flotacion(i,j)
         B_w_i(i,j)=gw(i,j)*fx_flotacion(i-1,j)+(1.0_dp-gw(i,j))*fx_flotacion(i,j)
         B_n_i(i,j)=gn(i,j)*fy_flotacion(i,j+1)+(1.0_dp-gn(i,j))*fy_flotacion(i,j)
         B_s_i(i,j)=gs(i,j)*fy_flotacion(i,j-1)+(1.0_dp-gs(i,j))*fy_flotacion(i,j)

         ! Promediado 2 veces
         B_u_ii(i,j)=(-1.0_dp*((1.0_dp-ge(i,j))*(B_e_i(i,j)*(x_r(i+1+nx*(j-1))-x_r(i+nx*(j-1)))*deltay(j))) &
         -((1.0_dp-gw(i,j))*(B_w_i(i,j)*(x_r(i+nx*(j-1))-x_r(i-1+nx*(j-1)))*deltay(j))))/vol(i,j)
         B_v_ii(i,j)=(-1.0_dp*((1.0_dp-gn(i,j))*(B_n_i(i,j)*(y_r(i+nx*(j))-y_r(i+nx*(j-1)))*deltax(i))) &
         -((1.0_dp-gs(i,j))*(B_s_i(i,j)*(y_r(i+nx*(j-1))-y_r(i+nx*(j-2)))*deltax(i))))/vol(i,j)

        end do
      end do

      ! Promediado 3 veces
      do j=2,ny-1
        do i=2,nx-1
         B_e_iii(i,j)=ge(i,j)*B_u_ii(i+1,j)+(1.0_dp-ge(i,j))*B_u_ii(i,j)
         B_w_iii(i,j)=gw(i,j)*B_u_ii(i-1,j)+(1.0_dp-gw(i,j))*B_u_ii(i,j)
         B_n_iii(i,j)=gn(i,j)*B_u_ii(i,j+1)+(1.0_dp-gn(i,j))*B_u_ii(i,j)
         B_s_iii(i,j)=gs(i,j)*B_u_ii(i,j-1)+(1.0_dp-gs(i,j))*B_u_ii(i,j)
        end do
      end do

  end subroutine fuerzas_flotacion

  ! TODO: añadir el termino para las fuerzas de cuerpo

  subroutine actualizar_flujo_de_masa()

    call coef_d()
    call fuerzas_flotacion()

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

         ue(i,j)=ue_i(i,j)-de_i*((Pstar(i+1,j)-Pstar(i,j))/(x_r(i+1+nx*(j-1))-x_r(i+nx*(j-1)))&
              -gPstar_e_i)+(1.0_dp-lambdau)*(ue_n(i,j)-ue_i_n(i,j)) &
              +de_i*(B_e_i(i,j)-B_e_iii(i,j))

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

           uw(i,j)=uw_i(i,j)-dw_i*((Pstar(i,j)-Pstar(i-1,j))/(x_r(i+nx*(j-1))-x_r(i-1+nx*(j-1)))&
                -gPstar_w_i)+(1.0_dp-lambdau)*(uw_n(i,j)-uw_i_n(i,j)) &
                +dw_i*(B_w_i(i,j)-B_w_iii(i,j))

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

           vn(i,j)=vn_i(i,j)-dn_i*((Pstar(i,j+1)-Pstar(i,j))/(y_r(i+nx*(j))-y_r(i+nx*(j-1)))&
                -gPstar_n_i)+(1.0_dp-lambdav)*(vn_n(i,j)-vn_i_n(i,j)) &
                +dn_i*(B_n_i(i,j)-B_n_iii(i,j))

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

           vs(i,j)=vs_i(i,j)-ds_i*((Pstar(i,j)-Pstar(i,j-1))/(y_r(i+nx*(j-1))-y_r(i+nx*(j-2)))&
                -gPstar_s_i)+(1.0_dp-lambdav)*(vs_n(i,j)-vs_i_n(i,j)) &
                +ds_i*(B_s_i(i,j)-B_s_iii(i,j))

        end if

        me_star(i,j)=ue(i,j)*deltay(j)
        mw_star(i,j)=-uw(i,j)*deltay(j)
        mn_star(i,j)=vn(i,j)*deltax(i)
        ms_star(i,j)=-vs(i,j)*deltax(i)

      end do
    end do

  end subroutine actualizar_flujo_de_masa

end module flujo_de_masa
