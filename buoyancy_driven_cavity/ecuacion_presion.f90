module ecuacion_presion
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use variables_globales
  use gauss_seidel
  use condiciones_frontera
  use utilidades
  implicit none

  ! Variables de apoyo
  real(dp) :: de_i,dw_i,dn_i,ds_i

  public :: ecuacion_presion_
  private

  contains

  subroutine ecuacion_presion_()


    do j=2,ny-1
      do i=2,nx-1

         if (i==nx-1) then
            de_i=0.0_dp
         else
            de_i=0.5_dp*deltax(i)*deltay(j)*(1.0_dp/ap_u(i,j)+1.0_dp/ap_u(i+1,j))
         end if

         if (i==2) then
            dw_i=0.0_dp
         else
            dw_i=0.5_dp*deltax(i)*deltay(j)*(1.0_dp/ap_u(i,j)+1.0_dp/ap_u(i-1,j))
         end if

         if (j==ny-1) then
            dn_i=0.0_dp
         else
            dn_i=0.5_dp*deltax(i)*deltay(j)*(1.0_dp/ap_v(i,j)+1.0_dp/ap_v(i,j+1))
         end if

         if (j==2) then
            ds_i=0.0_dp
         else
            ds_i=0.5_dp*deltax(i)*deltay(j)*(1.0_dp/ap_v(i,j)+1.0_dp/ap_v(i,j-1))
         end if

         ae_p(i,j)=-de_i*deltay(j)/(x(i+1)-x(i))
         aw_p(i,j)=-dw_i*deltay(j)/(x(i)-x(i-1))
         an_p(i,j)=-dn_i*deltax(i)/(y(j+1)-y(j))
         as_p(i,j)=-ds_i*deltax(i)/(y(j)-y(j-1))
         ap_p(i,j)=-(ae_p(i,j)+aw_p(i,j)+an_p(i,j)+as_p(i,j))
         b_p(i,j)=-(me_star(i,j)+mw_star(i,j)+mn_star(i,j)+ms_star(i,j))

      end do
    end do

    call gauss_seidel_(ap_p,ae_p,aw_p,an_p,as_p,b_p,Pprime)

    call actualizar_presion_fronteras(Pprime)

  end subroutine ecuacion_presion_

end module ecuacion_presion
