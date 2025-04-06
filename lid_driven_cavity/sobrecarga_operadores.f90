module sobrecarga_operadores
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use variables_globales
  implicit none

  public :: assignment(=)
  private

  ! Al menos uno de los dos elementos debe de ser de un tipo derivado
  type :: matriz
     real(dp), dimension(nx,ny) :: M
  end type matriz

  interface assignment(=)
     module procedure asignar_matrices
  end interface

  contains

  subroutine asignar_matrices(B,A)
    real(dp), dimension(nx,ny), intent(in) :: A
    type(matriz), dimension(nx,ny), intent(out) :: B

    do j=1,ny
      do i=1,nx
         B(i,j)%M=A(i,j)
      end do
    end do

  end subroutine asignar_matrices

end module sobrecarga_operadores
