module utilidades
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use variables_globales
  implicit none

  public :: interpolar,inter_vol
  private

  contains

  ! Esta subrutina calcula los factores de poderacion para una interpolacion
  ! basada en volumen (ver Moukalled pag 165 PDF)
  subroutine inter_vol()

    do j=2,ny-1
       do i=2,nx-1
          ge(i,j)=vol(i,j)/(vol(i,j)+vol(i+1,j))
          gw(i,j)=vol(i,j)/(vol(i,j)+vol(i-1,j))
          gn(i,j)=vol(i,j)/(vol(i,j)+vol(i,j+1))
          gs(i,j)=vol(i,j)/(vol(i,j)+vol(i,j-1))
       end do
    end do

  end subroutine inter_vol

  ! TODO: revisar esta funcion

  ! Interpolacion de variables en las caras
  function interpolar(psi,cara,coefA)
    real(dp) :: interpolar
    real(dp), dimension(nx,ny), intent(in) :: psi
    character(len=*), intent(in) :: cara
    logical :: coefA
    real(dp) :: A,B,lenght,delta_,Ic

    ! Inicializacion
    A=0.0_dp
    B=0.0_dp
    lenght=0.0_dp
    delta_=0.0_dp
    Ic=0.0_dp

    if (cara=='e' .or. cara=='n') then
       if (cara=='e') then
          lenght=x(i+1)-x(i)
          B=psi(i+1,j)
          delta_=deltax(i)
       else if (cara=='n') then
          lenght=y(j+1)-y(j)
          B=psi(i,j+1)
          delta_=deltay(j)
       end if
       A=psi(i,j)
    else if (cara=='w' .or. cara=='s') then
       if (cara=='w') then
          lenght=x(i)-x(i-1)
          A=psi(i-1,j)
          delta_=deltax(i)
       else if (cara=='s') then
          lenght=y(j)-y(j-1)
          A=psi(i,j-1)
          delta_=deltay(j)
       end if
       B=psi(i,j)
    else
       write(*,*)"CARA NO VALIDA"
    end if

    if (coefA .eqv. .TRUE.) then
       A=A**(-1.0_dp)
       B=B**(-1.0_dp)
    end if

    Ic=delta_/(2.0_dp*lenght)

    interpolar=Ic*A+(1.0_dp-Ic)*B

  end function

end module utilidades
