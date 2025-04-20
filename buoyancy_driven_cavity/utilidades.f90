module utilidades
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use variables_globales
  implicit none

  public :: interpolar,inicializacion,inter_vol
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

  ! Inicializacion de campos y variables
  subroutine inicializacion()

    ! Inicializacion de los coeficientes "ap" para evitar errores de punto
    ! flotante
    ap_u(:,:)=1.0_dp
    ap_v(:,:)=1.0_dp
    ap_p(:,:)=1.0_dp
    ap_T(:,:)=1.0_dp

    ! Inicializacion del campo supuesto de presion
    Pstar(:,:)=1.0_dp

    ! Flujos de masa inicializados en cero
    me_star(:,:)=0.0_dp
    mw_star(:,:)=0.0_dp
    mn_star(:,:)=0.0_dp
    ms_star(:,:)=0.0_dp

    ! Inicializacion del campo supuesto de velocidad
    u_star(:,:)=0.0_dp
    v_star(:,:)=0.0_dp

    u_old=u_star
    v_old=v_star

    ! Inicializacion del campo de temperaturas
    T(:,:)=0.0_dp

    ! Coeficientes agrupados
    ! fronteras norte y sur
    do i=1,nx
      ! Ecuacion de momentum
      b_u(i,ny)=0.0_dp
      b_u(i,1)=0.0_dp
      b_v(i,ny)=0.0_dp
      b_v(i,1)=0.0_dp

      ! Ecuacion de presion
      as_p(i,ny)=1.0_dp
      an_p(i,1)=1.0_dp

      ! Ecuacion de energia
      as_T(i,ny)=1.0_dp
      an_T(i,1)=1.0_dp
    end do

    ! fronteras este y oeste
    do j=2,ny-1
      ! Ecuacion de momentum
      b_u(nx,j)=0.0_dp
      b_u(1,j)=0.0_dp
      b_v(nx,j)=0.0_dp
      b_v(1,j)=0.0_dp

      ! Ecuacion de presion
      aw_p(nx,j)=1.0_dp
      ae_p(1,j)=1.0_dp

      ! Ecuacion de energia
      b_T(nx,j)=T_F
      b_T(1,j)=T_C
    end do

  end subroutine inicializacion

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
