!En este modulo se imprimen los resultados para una cierta variable
!"psi" junto con las coordenadas de la respectiva malla en un archivo
!con extensión .dat
module escritura
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use variables_globales
  implicit none

  public :: escritura_
  private

  contains

  subroutine escritura_(id,nombre_archivo,nombre_variable,psi)
    integer, intent(in) :: id !identificador
    !len=* es para un string de tamaño arbitrario
    character(len=*) :: nombre_archivo
    character(len=*) :: nombre_variable
    real(dp), dimension(nx,ny), intent(in) :: psi

    write(*,*) "Escribiendo a "//trim(nombre_archivo)
    open(id,file=nombre_archivo,status='unknown',access='sequential')
    write(id,*)"title=",nombre_variable
    write(id,*)"variables=x,y,",nombre_variable
    write(id,*)"zone i=",nx,",j=",ny
    write(id,*)""
    do j=1,ny
      do i=1,nx
        write(id,*)x(i),y(j),psi(i,j)
      end do
    end do
    close(id)
    
  end subroutine escritura_

end module escritura
