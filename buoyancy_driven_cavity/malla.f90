!En este modulo se construye el dominio discretizado
module malla
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use variables_globales
  implicit none
  public :: mallado,mallado_rotado,volumenes
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

  ! Esta subrutina calcula los puntos de la malla rotada en base al algulo
  ! theta
  subroutine mallado_rotado()

    ! Calculo de la distancia de los puntos al origen
    do j=1,ny
      do i=1,nx
        distancia(i+nx*(j-1))=(x(i)**2.0_dp+y(j)**2.0_dp)**0.5_dp
        ! write(*,*)"distancia(",i,j,")=",distancia(i+nx*(j-1))
      end do
    end do

    ! Calculo de los puntos rotados
    do j=1,ny
      do i=1,nx
        if (i==1 .and. j==1) then
          x_r(i+nx*(j-1))=0.0_dp
          y_r(i+nx*(j-1))=0.0_dp
        else
          xi(i+nx*(j-1))=asin(y(j)/distancia(i+nx*(j-1)))
          x_r(i+nx*(j-1))=cos(theta+xi(i+nx*(j-1)))*distancia(i+nx*(j-1))
          y_r(i+nx*(j-1))=sin(theta+xi(i+nx*(j-1)))*distancia(i+nx*(j-1))
        end if
      end do
    end do

  end subroutine mallado_rotado

  ! Esta subrutina calcula el volumen de las celdas computacionales
  subroutine volumenes()

    do j=1,ny
       do i=1,nx
          vol(i,j)=deltax(i)*deltay(j)
       end do
    end do

  end subroutine volumenes

end module malla
