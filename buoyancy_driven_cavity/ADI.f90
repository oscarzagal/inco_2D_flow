module ADI
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use variables_globales
  implicit none

  ! coeficiente "p", "q"
  real(dp), dimension(nx,ny) :: P,Q

  ! NOTE: tengo la hipotesis de que se requiere la informacion de los nodos
  ! frontera para que el metodo sea efectivo. En caso contrario produce
  ! cualquier cochinada.

  public :: ADI_
  private

  contains

  subroutine ADI_(ap,ae,aw,an,as,b,phinew)
    real(dp), dimension(nx,ny), intent(in) :: ap,ae,aw,an,as,b
    real(dp), dimension(nx,ny), intent(inout) :: phinew


    ! Barrido por LBL-X
    do j=1,ny
      P(1,j)=ae(1,j)/ap(1,j)
      Q(1,j)=(an(1,j)*phinew(1,j+1)+as(1,j)*phinew(1,j-1)+b(1,j))/ &
      ap(1,j)
      do i=2,nx
        P(i,j)=ae(i,j)/(ap(i,j)-aw(i,j)*P(i-1,j))
        Q(i,j)=((an(i,j)*phinew(i,j+1)+as(i,j)*phinew(i,j-1)+ &
        b(i,j))+aw(i,j)*Q(i-1,j))/(ap(i,j)-aw(i,j)*P(i-1,j))
      end do

      phinew(nx,j)=Q(nx,j)

      do i=nx-1,1,-1
        phinew(i,j)=phinew(i+1,j)*P(i,j)+Q(i,j)
      end do
   end do

   ! Barrido por LBL-Y
    do i=1,nx
      P(i,1)=an(i,1)/ap(i,1)
      Q(i,1)=(ae(i,1)*phinew(i+1,1)+aw(i,1)*phinew(i-1,1)+ &
      b(i,1))/ap(i,1)
      do j=2,ny
        P(i,j)=an(i,j)/(ap(i,j)-as(i,j)*P(i,j-1))
        Q(i,j)=(ae(i,j)*phinew(i+1,j)+aw(i,j)*phinew(i-1,j)+ &
        b(i,j)+as(i,j)*Q(i,j-1))/(ap(i,j)-as(i,j)*P(i,j-1))
      end do

      phinew(i,ny)=Q(i,ny)

      do j=ny-1,1,-1
        phinew(i,j)=phinew(i,j+1)*P(i,j)+Q(i,j)
      end do
    end do

  end subroutine ADI_

end module ADI
