program main
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use variables_globales
  use malla
  use escritura, only:escritura_
  use condiciones_frontera
  use conductancia_difusiva, only:conductancia_difusiva_
  use gradiente_explicito, only:gradiente_explicito_
  use ecuacion_momentum, only:ecuacion_momentum_
  use ecuacion_energia, only:ecuacion_energia_
  use flujo_de_masa, only:actualizar_flujo_de_masa
  use ecuacion_presion, only:ecuacion_presion_
  use correccion, only:corregir_presion_velocidad_flujo
  use criterio_convergencia, only:convergencia
  use reasignacion, only:reasignacion_variables
  use sobrecarga_operadores
  use utilidades

  ! Generacion de la malla computacional
  call mallado(nx,hx,x,deltax)
  call mallado(ny,hy,y,deltay)

  ! Calculo de los volumenes de las celdas
  call volumenes()

  ! Calculo de los factores de ponderacion de la interpolacion por volumen
  call inter_vol()

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

  ! Asignacion de las velocidades en las fronteras
  ! Ecuacion de momentum en "u"
  call calcular_velocidades_fronteras(u_star,Uo)

  ! Ecuacion de momentum en "v"
  call calcular_velocidades_fronteras(v_star,0.0_dp)

  ! Condiciones de frontera para la ecuacion de la energia
  call temperaturas_fronteras()

  T_old=T

  ! Terminos difusivos (solo es necesario calcularlos al principio)
  ! Ecuacion de momentum
  call conductancia_difusiva_(nu,FluxFe_dif_V,FluxFw_dif_V,FluxFn_dif_V, &
  FluxFs_dif_V)

  ! Ecuacion de energia
  call conductancia_difusiva_(alpha,FluxFe_dif_T,FluxFw_dif_T,FluxFn_dif_T &
  ,FluxFs_dif_T)



  ! Inicio del bucle SIMPLE
do while(error_mayor .gt. epsilon)


  ! Campo Pprime iniciado en cero cada iteracion
  Pprime(:,:)=0.0_dp

  ! Resolucion de la ecuacion de momentum
  call ecuacion_momentum_()

  ! Se actualiza el flujo de masa en el dominio
  call actualizar_flujo_de_masa()

  ! Se resuelve la ecuacion de correccion de presion
  call ecuacion_presion_()

    ! write(*,*)"Antes de la correccion"
    ! do j=1,ny
    !   do i=1,nx
    !     write(*,*)"v_star(",i,j,") =",v_star(i,j)
    !     ! write(*,*)"T(",i,j,") =",T(i,j)
    !   end do
    ! end do
    ! write(*,*)" "

  ! Correccion de los campos
  call corregir_presion_velocidad_flujo()

    ! write(*,*)"Despues de la correccion"
    ! do j=1,ny
    !   do i=1,nx
    !     write(*,*)"v_star(",i,j,") =",v_star(i,j)
    !     ! write(*,*)"T(",i,j,") =",T(i,j)
    !   end do
    ! end do
    ! write(*,*)" "

    ! write(*,*)"********************************************"
    ! write(*,*)"T(3,2) =",T(3,2)
    ! write(*,*)"T(1,2) =",T(1,2)
    ! write(*,*)"T(2,3) =",T(2,3)
    ! write(*,*)"T(2,1) =",T(2,1)
    ! write(*,*)"********************************************"

    ! Resolucion de la ecuacion de la energia
    call ecuacion_energia_()

    ! write(*,*)"Coeficiente ap_v"
    ! do j=1,ny
    !   do i=1,nx
    !     write(*,*)"ap_v(",i,j,") =",ap_v(i,j)
    !   end do
    ! end do
    ! write(*,*)" "

    ! write(*,*)"Coeficiente ae_v"
    ! do j=1,ny
    !   do i=1,nx
    !     write(*,*)"ae_v(",i,j,") =",ae_v(i,j)
    !   end do
    ! end do
    ! write(*,*)" "

    ! write(*,*)"Coeficiente aw_v"
    ! do j=1,ny
    !   do i=1,nx
    !     write(*,*)"aw_v(",i,j,") =",aw_v(i,j)
    !   end do
    ! end do
    ! write(*,*)" "

    ! write(*,*)"Coeficiente an_v"
    ! do j=1,ny
    !   do i=1,nx
    !     write(*,*)"an_v(",i,j,") =",an_v(i,j)
    !   end do
    ! end do
    ! write(*,*)" "

    ! write(*,*)"Coeficiente as_v"
    ! do j=1,ny
    !   do i=1,nx
    !     write(*,*)"as_v(",i,j,") =",as_v(i,j)
    !   end do
    ! end do
    ! write(*,*)" "

    ! write(*,*)"Coeficiente b_v"
    ! do j=1,ny
    !   do i=1,nx
    !     write(*,*)"b_v(",i,j,") =",b_v(i,j)
    !   end do
    ! end do
    ! write(*,*)" "

    ! if (numit==2) exit

    ! do j=1,ny
    !   do i=1,nx
    !     write(*,*)"T_old(",i,j,")",T_old(i,j)
    !   end do
    ! end do

    ! do j=1,ny
    !   do i=1,nx
    !     write(*,*)"T(",i,j,")",T(i,j)
    !   end do
    ! end do

  ! write(*,*)"ap_T(2,2) =",ap_T(2,2)
  ! write(*,*)"ae_T(2,2) =",ae_T(2,2)
  ! write(*,*)"aw_T(2,2) =",aw_T(2,2)
  ! write(*,*)"an_T(2,2) =",an_T(2,2)
  ! write(*,*)"as_T(2,2) =",as_T(2,2)
  ! write(*,*)"b_T(2,2) =",b_T(2,2)
  ! write(*,*)"T(2,2) =",T(2,2)
  ! write(*,*)"mw_star(2,2) =",mw_star(2,2)
  ! write(*,*)"alpha =",alpha

  ! if (.true.) exit


  ! Convergencia
  call convergencia(Pstar,Pold,residual_presion)
  call convergencia(u_star,u_old,residual_u)
  call convergencia(v_star,v_old,residual_v)
  call convergencia(me_star,me_old,residual_me)
  call convergencia(mw_star,mw_old,residual_mw)
  call convergencia(mn_star,mn_old,residual_mn)
  call convergencia(ms_star,ms_old,residual_ms)
  call convergencia(T,T_old,residual_energia)

  ! Eleccion del error mayor
  error_mayor=resnew(residual_presion)
  variable=residual_presion
  do i=2,8
    if (resnew(i).gt.error_mayor) then
      error_mayor=resnew(i)
      variable=i
    else
      error_mayor=error_mayor
    end if
  end do

  numit=numit+1

  ! Identificador de variable
  if (variable==1) then
      nombre="Presion"
  else if (variable==2) then
      nombre="u"
  else if (variable==3) then
      nombre="v"
  else if (variable==4) then
      nombre="me"
  else if (variable==5) then
      nombre="mw"
  else if (variable==6) then
      nombre="mn"
  else if (variable==7) then
      nombre="ms"
  else if (variable==8) then
      nombre="Energia"
  end if

  ! do j=1,ny
  !    do i=1,nx
  !       write(*,*)"∑ \dot{m}_{f} =",me_star(i,j)+mw_star(i,j)+mn_star(i,j)+ms_star(i,j)
  !    end do
  ! end do

  write(*,*)numit,error_mayor," variable: ", nombre

  call reasignacion_variables()

  ! Segundo criterio de paro
  if (numit==limite) then
    exit
  end if

end do ! Fin del bucle SIMPLE

  write(*,*)" "
  write(*,*)"deltat falso =",deltat
  write(*,*)"Difusividad termica =",alpha
  write(*,*)"Numero de Prandtl =",nu/alpha
  write(*,*)"Numero de Rayleigh =", Ra
  write(*,*)" "


  ! Magnitud de la velocidad
  do j=1,ny
    do i=1,nx
      Umag(i,j)=(u_star(i,j)**(2.0_dp)+v_star(i,j)**(2.0_dp))**0.5_dp
    end do
  end do

  call escritura_(90,'Pstar.dat','Pstar',Pstar)
  call escritura_(90,'Umag.dat','Umag',Umag)
  call escritura_(10,'me.dat','me',me_star)
  call escritura_(76, 'T.dat', 'T', T)

end program main
