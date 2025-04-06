module esquema_upwind
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use variables_globales
  implicit none

  public :: esquema_upwind_
  private

  contains

  subroutine esquema_upwind_(ap,ae,aw,an,as,b,gPstar,vel,lambda)
    real(dp), dimension(nx,ny), intent(out) :: ap,ae,aw,an,as,b
    real(dp), dimension(nx,ny), intent(in) :: gPstar,vel
    real(dp), intent(in) :: lambda

    do j=2,ny-1
      do i=2,nx-1

        ae(i,j)=-De(i,j)-max(0.0_dp,-me_star(i,j))
        aw(i,j)=-Dw(i,j)-max(-mw_star(i,j),0.0_dp)
        an(i,j)=-Dn(i,j)-max(0.0_dp,-mn_star(i,j))
        as(i,j)=-Ds(i,j)-max(-ms_star(i,j),0.0_dp)

        ap(i,j)=(De(i,j)+max(0.0_dp,me_star(i,j)) &
             +Dw(i,j)+max(mw_star(i,j),0.0_dp) &
             +Dn(i,j)+max(0.0_dp,mn_star(i,j)) &
             +Ds(i,j)+max(ms_star(i,j),0.0_dp) &
             )/lambda

        b(i,j)=-gPstar(i,j)+(1.0_dp-lambda)*ap(i,j)*vel(i,j)

      end do
    end do

  end subroutine esquema_upwind_

end module esquema_upwind
