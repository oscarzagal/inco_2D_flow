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

  implicit none

  ! Generacion de la malla computacional
  call mallado(nx,hx,x,deltax)
  call mallado(ny,hy,y,deltay)

  ! Calculo de los volumenes de las celdas
  call volumenes()

  ! Calculo de los factores de ponderacion de la interpolacion por volumen
  call inter_vol()

  ! Generacion de los puntos rotados
  call mallado_rotado()


  ! Inicializacion de campos y variables
  call inicializacion()


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

  ! Correccion de los campos
  call corregir_presion_velocidad_flujo()

    ! Resolucion de la ecuacion de la energia
    call ecuacion_energia_()

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

  write(*,*)numit,error_mayor," variable: ", nombre

  call reasignacion_variables()

  ! Segundo criterio de paro
  if (numit==limite) then
    exit
  end if

end do ! Fin del bucle SIMPLE

  write(*,*)" "
  write(*,*)"Angulo de rotacion =",theta
  write(*,*)"deltat falso =",deltat
  write(*,*)"Difusividad termica =",alpha
  write(*,*)"Numero de Prandtl =",nu/alpha
  write(*,*)"Numero de Rayleigh =", Ra
  write(*,*)" "

  ! do j=1,ny
  !    do i=1,nx
  !       write(*,*)"as_T(",i,j,")",as_T(i,j)
  !    end do
  ! end do


  ! Suma del flujo de masa
  do j=1,ny
     do i=1,nx
        sumMdot(i,j)=me_star(i,j)+mw_star(i,j)+mn_star(i,j)+ms_star(i,j)
     end do
  end do


  ! Magnitud de la velocidad
  do j=1,ny
    do i=1,nx
      Umag(i,j)=(u_star(i,j)**(2.0_dp)+v_star(i,j)**(2.0_dp))**0.5_dp
    end do
  end do

  call escritura_(90,'Pstar.dat','Pstar',Pstar)
  call escritura_(90,'Umag.dat','Umag',Umag)
  call escritura_(10,'sumMdot.dat','sumMdot',sumMdot)
  call escritura_(76, 'T.dat', 'T', T)
  call escritura_(73, 'u.dat', 'u', u_star)
  call escritura_(71, 'v.dat', 'v', v_star)

end program main
