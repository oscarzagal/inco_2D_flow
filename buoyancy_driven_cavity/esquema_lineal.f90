module esquema_lineal
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use variables_globales

  implicit none

  public :: esquema_lineal_
  private

  contains

  subroutine esquema_lineal_(ap,ae,aw,an,as,b,gPstar,vel,lambda)
    real(dp), dimension(nx,ny), intent(out) :: ap,ae,aw,an,as,b
    real(dp), dimension(nx,ny), intent(in) :: gPstar,vel
    real(dp), intent(in) :: lambda

    do j=2,ny-1
      do i=2,nx-1
        ae(i,j)=De(i,j)-0.5_dp*me_star(i,j)
        aw(i,j)=Dw(i,j)+0.5_dp*mw_star(i,j)
        an(i,j)=Dn(i,j)-0.5_dp*mn_star(i,j)
        as(i,j)=Ds(i,j)+0.5_dp*ms_star(i,j)
        ap(i,j)=(ae(i,j)+aw(i,j)+an(i,j)+as(i,j)+me_star(i,j)-mw_star(i,j)&
        +mn_star(i,j)-ms_star(i,j))/lambda
        b(i,j)=gPstar(i,j)+(1.0_dp-lambda)*ap(i,j)*vel(i,j)
      end do
    end do

  end subroutine esquema_lineal_

end module esquema_lineal
