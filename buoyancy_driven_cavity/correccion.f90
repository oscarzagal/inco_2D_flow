module correccion
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use variables_globales
  use condiciones_frontera
  use utilidades

  implicit none

  ! Variables de apoyo
  real(dp), dimension(nx,ny) :: gPprime_u,gPprime_v

  public :: corregir_presion_velocidad_flujo
  private

  contains

  subroutine gradiente_Pprime_centrado_celda()

    do j=2,ny-1
      do i=2,nx-1

         gPprime_u(i,j)=(ge(i,j)*Pprime(i+1,j)+(1.0_dp-ge(i,j))*Pprime(i,j)-gw(i,j)&
              *Pprime(i-1,j)-(1.0_dp-gw(i,j))*Pprime(i,j))/deltax(i)
         gPprime_v(i,j)=(gn(i,j)*Pprime(i,j+1)+(1.0_dp-gn(i,j))*Pprime(i,j)-gs(i,j)&
              *Pprime(i,j-1)-(1.0_dp-gs(i,j))*Pprime(i,j))/deltay(j)

      end do
    end do

  end subroutine gradiente_Pprime_centrado_celda

  ! TODO: refactorizar la correccion del flujo de masa

  subroutine corregir_presion_velocidad_flujo()

    call gradiente_Pprime_centrado_celda()

    do j=2,ny-1
      do i=2,nx-1
          ! Correccion de velocidades
          u_star(i,j)=u_star(i,j)-dP_u(i,j)*gPprime_u(i,j)
          v_star(i,j)=v_star(i,j)-dP_v(i,j)*gPprime_v(i,j)

          ! Correccion de flujo de masa
          if (i==nx-1) then
             me_star(i,j)=0.0_dp
          else

             me_star(i,j)=me_star(i,j)+0.5_dp*(1.0_dp/ap_u(i,j)&
                  +1.0_dp/ap_u(i+1,j))*(Pprime(i,j)-Pprime(i+1,j))/(x(i+1)-x(i))&
                  *(deltay(j)**2.0_dp)*deltax(i)
          end if

          if (i==2) then
             mw_star(i,j)=0.0_dp
          else

             mw_star(i,j)=mw_star(i,j)+0.5_dp*(1.0_dp/ap_u(i,j)&
                  +1.0_dp/ap_u(i-1,j))*(Pprime(i,j)-Pprime(i-1,j))/(x(i)-x(i-1))&
                  *(deltay(j)**2.0_dp)*deltax(i)
          end if

          if (j==ny-1) then
             mn_star(i,j)=0.0_dp
          else

             mn_star(i,j)=mn_star(i,j)+0.5_dp*(1.0_dp/ap_v(i,j)&
                  +1.0_dp/ap_v(i,j+1))*(Pprime(i,j)-Pprime(i,j+1))/(y(j+1)-y(j))&
                  *(deltax(i)**2.0_dp)*deltay(j)
          end if

          if (j==2) then
             ms_star(i,j)=0.0_dp
          else

             ms_star(i,j)=ms_star(i,j)+0.5_dp*(1.0_dp/ap_v(i,j)&
                  +1.0_dp/ap_v(i,j-1))*(Pprime(i,j)-Pprime(i,j-1))/(y(j)-y(j-1))&
                  *(deltax(i)**2.0_dp)*deltay(j)

          end if

          ! Correcion de presion
          ! if (i==2 .and. j==1) then
          !    Pstar(i,j)=0.0_dp
          ! else
              Pstar(i,j)=Pstar(i,j)+lambdaP*Pprime(i,j)
          ! end if

      end do
    end do

    call actualizar_presion_fronteras(Pstar)

  end subroutine corregir_presion_velocidad_flujo

end module correccion
