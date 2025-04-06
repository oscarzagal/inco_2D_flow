module conductancia_difusiva
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use variables_globales

  implicit none

  public :: conductancia_difusiva_
  private

  contains

  subroutine conductancia_difusiva_()

    do j=2,ny-1
      do i=2,nx-1
        De(i,j)=mu*deltay(j)/(x(i+1)-x(i))
        Dw(i,j)=mu*deltay(j)/(x(i)-x(i-1))
        Dn(i,j)=mu*deltax(i)/(y(j+1)-y(j))
        Ds(i,j)=mu*deltax(i)/(y(j)-y(j-1))
      end do
    end do

  end subroutine conductancia_difusiva_

end module conductancia_difusiva
