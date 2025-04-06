!En esta subrutina se calcula el residual de las variables de campo
module criterio_convergencia
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use variables_globales
  implicit none
  real(dp) :: res !residual

  public :: convergencia
  private

  contains

  subroutine convergencia(phinew,phiold,var)
    integer, intent(in) :: var
    real(dp), dimension(nx,ny), intent(in) :: phinew,phiold

    resnew(var)=0.0_dp
    do j=1,ny
      do i=1,nx
        res=abs(phinew(i,j)-phiold(i,j))
        resnew(var)=max(resnew(var),res)
      end do
    end do
    
  end subroutine convergencia

end module criterio_convergencia
