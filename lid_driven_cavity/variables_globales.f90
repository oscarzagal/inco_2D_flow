module variables_globales
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  implicit none

  ! Control
  integer :: esquema=1

  ! Variables para bucles
  integer :: i,j,k

  ! Subrelajacion
  ! Se demostro empiricamente que con estos factores de subrelajacion se
  ! consiguien resultados adecuados para numeros de Reynolds de 100 y 400.
  ! Ademas es necesario aumentar el numero de elementos computacionales para
  ! lograr la mejor concordancia con la referencia (pagina 252 PDF Xaman).
  ! Actualizacion: lo mas importante es la malla, para un Re > 100 o 400
  ! es importante que la malla sea muy fina, de lo contrario se obtienen
  ! resultados de baja calidad.
  real(dp), parameter :: lambdaP=0.3_dp
  real(dp), parameter :: lambdau=0.6_dp
  real(dp), parameter :: lambdav=0.6_dp

  ! Variables geometricas
  integer, parameter :: nx=171,ny=171
  ! real(dp), parameter :: hx=0.1_dp,hy=0.1_dp
  real(dp), parameter :: hx=1.0_dp,hy=1.0_dp
  real(dp), dimension(nx) :: x,deltax
  real(dp), dimension(ny) :: y,deltay
  real(dp), dimension(nx,ny) :: vol

  ! Paso de tiempo para el falso transitorio. Se probo que para 0.1 y 0.5
  ! los resultados son iguales, por lo que no hay dependencia del paso de
  ! tiempo falso.
  real(dp), parameter :: deltat=0.5_dp

  ! Interpolacion
  real(dp), dimension(nx,ny) :: ge,gw,gn,gs

  ! Limite de iteraciones
  integer, parameter :: limite=20000

  ! Tolerancia
  real(dp), parameter :: epsilon=1e-5_dp

  ! Propiedades termofisicas
  ! real(dp), parameter :: mu=1.817e-5_dp
  ! real(dp), parameter :: rho=1.2047_dp

  real(dp), parameter :: mu=2.5e-3_dp
  real(dp), parameter :: rho=1.0_dp

  ! Condicion de frontera de pared deslizante
  ! real(dp), parameter :: Uo=1.508e-3_dp
  real(dp), parameter :: Uo=1.0_dp

  ! Numero de Reynolds
  real(dp), parameter :: Re=(rho*hx*Uo)/mu

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

  ! Conductancia difusiva
  real(dp), dimension(nx,ny) :: De,Dw,Dn,Ds

  ! Flujo de masa en las caras
  real(dp), dimension(nx,ny) :: me_star,mw_star,mn_star,ms_star

  ! Identificadores para la subrutina del criterio de convergencia
  integer, parameter :: residual_presion=1
  integer, parameter :: residual_u=2
  integer, parameter :: residual_v=3
  integer, parameter :: residual_me=4
  integer, parameter :: residual_mw=5
  integer, parameter :: residual_mn=6
  integer, parameter :: residual_ms=7

  ! Residual actualizado
  real(dp), dimension(7) :: resnew


  ! Error mayor
  real(dp) :: error_mayor=1

  ! Variable de bucle
  integer :: numit=0

  ! Identificador de variables
  integer :: variable
  character(len=20) :: nombre

  ! Variables de la iteracion anterior
  real(dp), dimension(nx,ny) :: Pold,u_old,v_old
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
