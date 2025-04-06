module gradiente_explicito
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use variables_globales
  use utilidades
  implicit none

  public :: gradiente_explicito_
  private

  contains

  subroutine gradiente_explicito_()

    do j=2,ny-1
      do i=2,nx-1

         gPstar_u_vol(i,j)=(ge(i,j)*Pstar(i+1,j)+(1.0_dp-ge(i,j))*Pstar(i,j)-gw(i,j)&
              *Pstar(i-1,j)-(1.0_dp-gw(i,j))*Pstar(i,j))*deltay(j)
         gPstar_v_vol(i,j)=(gn(i,j)*Pstar(i,j+1)+(1.0_dp-gn(i,j))*Pstar(i,j)-gs(i,j)&
              *Pstar(i,j-1)-(1.0_dp-gs(i,j))*Pstar(i,j))*deltax(i)

      end do
    end do

  end subroutine gradiente_explicito_

end module gradiente_explicito
