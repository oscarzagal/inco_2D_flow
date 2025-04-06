module gauss_seidel
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use variables_globales
  implicit none

  public :: gauss_seidel_
  private

  contains

  subroutine gauss_seidel_(ap,ae,aw,an,as,b,phinew)
    real(dp), dimension(nx,ny), intent(in) :: ap,ae,aw,an,as,b
    real(dp), dimension(nx,ny), intent(inout) :: phinew

    do j=2,ny-1
      do i=2,nx-1
            phinew(i,j)=(phinew(i+1,j)*(-ae(i,j))+phinew(i-1,j)*(-aw(i,j))+phinew(i &
                 ,j+1)*(-an(i,j))+phinew(i,j-1)*(-as(i,j))+b(i,j))/ap(i,j)
      end do
   end do


  end subroutine gauss_seidel_

end module gauss_seidel
