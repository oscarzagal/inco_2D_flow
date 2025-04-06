module esquema_ley_de_potencia
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use variables_globales
  implicit none

  public :: esquema_ley_de_potencia_
  private

  contains

  subroutine esquema_ley_de_potencia_(ap,ae,aw,an,as,b,gPstar,vel,lambda)
    real(dp), dimension(nx,ny), intent(out) :: ap,ae,aw,an,as,b
    real(dp), dimension(nx,ny), intent(in) :: gPstar,vel
    real(dp), intent(in) :: lambda

    do j=2,ny-1
      do i=2,nx-1
        ae(i,j)=De(i,j)*max(0.0_dp,(1.0_dp-0.1_dp*abs(me_star(i,j)/De(i,j))) &
        **5.0_dp)+max(0.0_dp,-me_star(i,j))
        aw(i,j)=Dw(i,j)*max(0.0_dp,(1.0_dp-0.1_dp*abs(mw_star(i,j)/Dw(i,j))) &
        **5.0_dp)+max(mw_star(i,j),0.0_dp)
        an(i,j)=Dn(i,j)*max(0.0_dp,(1.0_dp-0.1_dp*abs(mn_star(i,j)/Dn(i,j))) &
        **5.0_dp)+max(0.0_dp,-mn_star(i,j))
        as(i,j)=Ds(i,j)*max(0.0_dp,(1.0_dp-0.1_dp*abs(ms_star(i,j)/Ds(i,j))) &
        **5.0_dp)+max(ms_star(i,j),0.0_dp)
        ap(i,j)=(ae(i,j)+aw(i,j)+an(i,j)+as(i,j)+me_star(i,j)-mw_star(i,j) &
        +mn_star(i,j)-ms_star(i,j))/lambda
        b(i,j)=gPstar(i,j)+(1.0_dp-lambda)*ap(i,j)*vel(i,j)
      end do
    end do

  end subroutine esquema_ley_de_potencia_

end module esquema_ley_de_potencia
