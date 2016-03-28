Program average
implicit none

double precision:: cumul, cumul2, aver, aver2, sigma2, minim, maxim
double precision,dimension(:),allocatable::nada
double precision::dato
integer::i,j,k,columna,num_data
CHARACTER(25)::input

double precision::num_div,delta,delta_2,prueba
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::frecuencias
integer::num_div2

open (13, FILE="entrada.dat")
open (100,FILE="histograma.oup")
read (13,*) input
read (13,*) num_data
read (13,*) columna
read (13,*) delta
open (12, FILE=TRIM(input), STATUS="OLD")


dato=0.0d0
cumul=0.0d0
cumul2=0.0d0
aver2=0.0d0
aver=0.0d0
sigma2=0.0d0
minim=0.0d0
maxim=0.0d0


if ( columna > 1 ) then

   allocate (nada(columna-1))

   read(12,*) nada, maxim
   
   minim=maxim
   
   
   rewind(12)
   
   do i=1,num_data
      
      read(12,*) nada, dato
      
      if ( dato > maxim ) maxim=dato
      if ( dato < minim ) minim=dato
      
      cumul=cumul+dato
      cumul2=cumul2+dato**2
      
      
   end do

end if

if ( columna == 1 ) then


   read(12,*) maxim

   minim=maxim


   rewind(12)

   do i=1,num_data

      read(12,*) dato

      if ( dato > maxim ) maxim=dato
      if ( dato < minim ) minim=dato

      cumul=cumul+dato
      cumul2=cumul2+dato**2


   end do

end if



   
aver=cumul/num_data
aver2=cumul2/num_data
sigma2=aver2-aver**2

print*, "- Datos de la columna", columna," de ",TRIM(input),":"
print*, "   <X> =", aver, "<X^2>-<X>^2 =", sigma2
print*, "   X_minimo =", minim, "X_maximo =", maxim
print*, " "
print*, "- Histograma en histograma.dat"

!histograma:



delta_2=delta/2
num_div=((maxim-minim)/delta)+1
num_div2=num_div

ALLOCATE(frecuencias(num_div2))

frecuencias=0.0d0

REWIND(12)

DO k=1,num_data

   read(12,*) nada, dato

   DO j=1,num_div2
      
      IF (dato < (maxim+delta_2-(j-1)*delta)) THEN
         
         frecuencias(j)=frecuencias(j)+1
         
      ELSE
         
         EXIT
         
      END IF
      
   END DO
END DO


prueba=0

DO j=1,num_div2-1
   
   frecuencias(j)=frecuencias(j)-frecuencias(j+1)
   
   prueba=prueba+frecuencias(j)
END DO

prueba=prueba+frecuencias(num_div2)


!normalizo al area:

prueba=0

DO j=1,num_div2
   prueba=prueba+frecuencias(j)*delta
END DO


frecuencias(:)=frecuencias(:)/prueba


DO i=1,num_div2
   WRITE(100,*) maxim-(i-1)*delta,frecuencias(i)
END DO





end program
