module esquema_upwind
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use variables_globales
  implicit none

  public :: esquema_upwind_
  private

  contains

  subroutine esquema_upwind_(ap,ae,aw,an,as,b,gPstar,f_flotacion,vel,lambda)
    real(dp), dimension(nx,ny), intent(out) :: ap,ae,aw,an,as,b
    real(dp), dimension(nx,ny), intent(in) :: gPstar,f_flotacion,vel
    real(dp), intent(in) :: lambda

    do j=2,ny-1
      do i=2,nx-1

        ae(i,j)=-FluxFe_dif_V(i,j)-max(0.0_dp,-me_star(i,j))
        aw(i,j)=-FluxFw_dif_V(i,j)-max(-mw_star(i,j),0.0_dp)
        an(i,j)=-FluxFn_dif_V(i,j)-max(0.0_dp,-mn_star(i,j))
        as(i,j)=-FluxFs_dif_V(i,j)-max(-ms_star(i,j),0.0_dp)

        ap(i,j)=(FluxFe_dif_V(i,j)+max(0.0_dp,me_star(i,j)) &
             +FluxFw_dif_V(i,j)+max(mw_star(i,j),0.0_dp) &
             +FluxFn_dif_V(i,j)+max(0.0_dp,mn_star(i,j)) &
             +FluxFs_dif_V(i,j)+max(ms_star(i,j),0.0_dp) &
             +vol(i,j)/deltat)/lambda

        b(i,j)=-gPstar(i,j)+f_flotacion(i,j)+(1.0_dp-lambda)*ap(i,j)*vel(i,j) &
        +vel(i,j)*vol(i,j)/deltat

      end do
    end do

  end subroutine esquema_upwind_

end module esquema_upwind
