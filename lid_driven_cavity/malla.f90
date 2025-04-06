!En este modulo se construye el dominio discretizado
module malla
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use variables_globales
  implicit none
  public :: mallado,volumenes
  private

  contains

  ! Esta subrutina obtiene los puntos y deltas que conforman la malla
  ! computacional
  subroutine mallado(nn,h,c,delta)
    integer,intent(in) :: nn !dato de entrada
    real(dp),intent(in) :: h !valor que no va a ser modificado
    real(dp),dimension(nn),intent(out) :: c  !dato de salida
    real(dp),dimension(nn),intent(out) :: delta !valor que si puede modificar
    real(dp) :: del

    del=h/(nn-2.0_dp)

    delta(1)=0.0_dp
    delta(nn)=0.0_dp

    !asignacion de las coordenadas al vector c
    c(1)=0.0_dp
    c(nn)=h

    do i=2,nn-1
      delta(i)=del
      c(i)=c(i-1)+(delta(i-1)+delta(i))/2.0_dp
    end do
    
  end subroutine mallado

  ! Esta subrutina calcula el volumen de las celdas computacionales
  subroutine volumenes()

    do j=1,ny
       do i=1,nx
          vol(i,j)=deltax(i)*deltay(j)
       end do
    end do

  end subroutine volumenes

end module malla
