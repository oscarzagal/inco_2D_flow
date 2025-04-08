FC=gfortran
FFLAGS=-O3 -Wall -Wextra -std=f2008 #-Wextra
#FFLAGS=-g -ffpe-trap=zero,invalid,overflow,underflow #debuggear
#FFLAGS=-g #debuggear
 #el orden de compilacion importa
SRC=variables_globales.f90\
		sobrecarga_operadores.f90\
		malla.f90\
		escritura.f90\
		utilidades.f90\
		condiciones_frontera.f90\
		gradiente_explicito.f90\
		gauss_seidel.f90\
		conductancia_difusiva.f90\
		esquema_upwind.f90\
		ecuacion_momentum.f90\
		flujo_de_masa.f90\
		ecuacion_presion.f90\
		ecuacion_energia.f90\
		correccion.f90\
		criterio_convergencia.f90\
		reasignacion.f90\
		main.f90

OBJ=${SRC:.f90=.o}

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

main: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

clean:
	@rm -f *.mod *.o main
