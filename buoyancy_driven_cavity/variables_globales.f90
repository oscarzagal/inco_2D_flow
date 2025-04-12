module variables_globales
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  implicit none

  ! Variables para bucles
  integer :: i,j,k

  ! Subrelajacion
  real(dp), parameter :: lambdaP=0.3_dp
  real(dp), parameter :: lambdau=0.6_dp
  real(dp), parameter :: lambdav=0.6_dp
  real(dp), parameter :: lambdaT=0.1_dp

  ! Variables geometricas
  integer, parameter :: nx=55,ny=55
  real(dp), parameter :: hx=0.05_dp,hy=0.05_dp
  real(dp), dimension(nx) :: x,deltax
  real(dp), dimension(ny) :: y,deltay
  real(dp), dimension(nx,ny) :: vol

  ! Paso de tiempo para el falso transitorio
  real(dp), parameter :: deltat=0.35_dp

  ! Interpolacion
  real(dp), dimension(nx,ny) :: ge,gw,gn,gs

  ! Limite de iteraciones
  integer, parameter :: limite=50000

  ! Tolerancia
  real(dp), parameter :: epsilon=1e-5_dp


  ! Propiedades termofisicas

  ! Viscosidad dinamica
  real(dp), parameter :: nu=1.562e-5_dp

  ! Temperatura de referencia
  real(dp), parameter :: T_inf=25.0_dp

  ! Coeficiente de expansividad termica
  real(dp), parameter :: beta=1.0_dp/(25.0_dp+273.15_dp)

  ! Densidad
  real(dp), parameter :: rho=1.0_dp

  ! Calor espefico
  real(dp), parameter :: Cp=7006.071_dp

  ! Conductividad termica
  real(dp), parameter :: k_ter=0.15_dp

  ! Difusividad termica
  real(dp), parameter :: alpha=k_ter/(rho*Cp)



  ! Condiciones de frontera para la ecuacion de la energia
  real(dp), parameter :: T_C=25.0814_dp
  real(dp), parameter :: T_F=25.0_dp

  ! Gravedad
  real(dp), parameter :: g=-9.81_dp

  ! Angulo de rotacion de la cavidad
  real(dp), parameter :: theta=0.0_dp

  ! Fuerza de flotacion
  real(dp), dimension(nx,ny) :: fx_flotacion,fy_flotacion

  ! Numero de Rayleigh
  real(dp), parameter :: Ra=-(g*beta*(T_C-T_F)*hx**3.0_dp)/(alpha*nu)

  ! Condicion de frontera de pared deslizante
  real(dp), parameter :: Uo=0.0_dp


  ! Coeficientes agrupados para la ecuacion de momentum
  ! Direccion "x"
  real(dp), dimension(nx,ny) :: ap_u,ae_u,aw_u,an_u,as_u,b_u

  ! Direccion "y"
  real(dp), dimension(nx,ny) :: ap_v,ae_v,aw_v,an_v,as_v,b_v

  ! Coeficientes agrupados para la ecuacion de correccion de presion
  real(dp), dimension(nx,ny) :: ap_p,ae_p,aw_p,an_p,as_p,b_p

  ! Velocidades en los centroides de los elementos
  real(dp), dimension(nx,ny) :: u_star,v_star

  ! Variables relacionadas con la presion
  ! (presion supuesta, gradientes en "x" e "y", presion corregida)
  real(dp), dimension(nx,ny) :: Pstar,gPstar_u_vol,gPstar_v_vol,Pprime

  ! Coeficientes agrupados para la ecuacion de energia
  real(dp), dimension(nx,ny) :: ap_T,ae_T,aw_T,an_T,as_T,b_T

  ! Campo de temperatura
  real(dp), dimension(nx,ny) :: T

  ! Conductancia difusiva (ecuacion de momentum)
  real(dp), dimension(nx,ny) :: FluxFe_dif_V,FluxFw_dif_V,FluxFn_dif_V,FluxFs_dif_V

  ! Conductancia difusiva (ecuacion de energia)
  real(dp), dimension(nx,ny) :: FluxFe_dif_T,FluxFw_dif_T,FluxFn_dif_T,FluxFs_dif_T

  ! Flux convectivo (coeficiente ap)
  real(dp), dimension(nx,ny) :: FluxCe_conv,FluxCw_conv,FluxCn_conv,FluxCs_conv

  ! Flux convectivo (demas coeficientes)
  real(dp), dimension(nx,ny) :: FluxFe_conv,FluxFw_conv,FluxFn_conv,FluxFs_conv


  ! Flujo de masa en las caras
  real(dp), dimension(nx,ny) :: me_star,mw_star,mn_star,ms_star


  ! Identificadores para la subrutina del criterio de convergencia
  ! (mi intento de Enum xd)
  integer, parameter :: residual_presion=1
  integer, parameter :: residual_u=2
  integer, parameter :: residual_v=3
  integer, parameter :: residual_me=4
  integer, parameter :: residual_mw=5
  integer, parameter :: residual_mn=6
  integer, parameter :: residual_ms=7
  integer, parameter :: residual_energia=8

  ! Residual actualizado
  real(dp), dimension(8) :: resnew


  ! Error mayor
  real(dp) :: error_mayor=1

  ! Variable de bucle
  integer :: numit=0

  ! Identificador de variables
  integer :: variable
  character(len=20) :: nombre

  ! Variables de la iteracion anterior
  real(dp), dimension(nx,ny) :: Pold,u_old,v_old,T_old
  real(dp), dimension(nx,ny) :: me_old,mw_old,mn_old,ms_old

  ! Velocidades en las caras
  real(dp), dimension(nx,ny) :: ue,uw,vn,vs

  ! Coeficiente d
  real(dp), dimension(nx,ny) :: dE_u,dW_u,dN_v,dS_v,dP_u,dP_v

  ! Velocidades interpoladas en las caras
  real(dp), dimension(nx,ny) :: ue_i,uw_i,vn_i,vs_i

  ! Velocidades en las caras de la iteracion anterior
  real(dp), dimension(nx,ny) :: ue_n,uw_n,vn_n,vs_n

  ! Velocidades interpoladas en las caras de la iteracion anterior
  real(dp), dimension(nx,ny) :: ue_i_n,uw_i_n,vn_i_n,vs_i_n

  ! Magnitud de la velocidad
  real(dp), dimension(nx,ny) :: Umag


end module variables_globales
