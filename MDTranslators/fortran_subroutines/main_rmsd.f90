program main

use para_rmsd
use pdb_in_out

implicit none

integer::N,i,j
double precision,dimension(:,:),allocatable::struct_ref,struct_2,posnew
double precision,dimension(3,3)::rot
double precision,dimension(3)::center_ref,center_2
double precision::rmsd
double precision,dimension(:,:),allocatable::g

character(len=long_max),dimension(2)::file
character(len=long_max)::argumento1,opcion,opcion2
character(len=long_max)::buffer
integer::num_args

integer::ios
integer,dimension(2)::num_atom

!#############

num_args=IARG()

IF (num_args==1) THEN
   
   CALL getarg(1,buffer)
   argumento1=TRIM(buffer)
   
   IF ((argumento1=="-help").OR.(argumento1=="-h")) THEN
      PRINT*," "
      PRINT*,"rmsd_fit estructura1.pdb estructura2.pdb opcion"
      PRINT*," "
      PRINT*,"El programa fittea la estructura2 sobre la estructura1 minimizando"
      PRINT*,"el rmsd segun las siguientes opciones:"
      PRINT*," "
      PRINT*,"   opcion = ca (los carbonos alpha)"
      PRINT*,"   opcion = backbone (los atomos del tipo CA, N y C)"
      PRINT*,"   opcion = all (todos los atomos)"
      PRINT*," "
      PRINT*,"*: todos los atomos tienen el mismo peso en el analisis del rmsd. "
      PRINT*," "
      STOP
   ELSE
      PRINT*,"ERROR en el numero de argumentos"
      PRINT*,"tip: ./rmsd.exe -help"
      STOP
   END IF
END IF

IF (num_args/=3) THEN
   PRINT*,"ERROR en el numero de argumentos"
   PRINT*,"hint: ./rmsd.exe -help"
   STOP
ELSE

   CALL getarg(1,buffer)
   file(1)=TRIM(buffer)
   
   CALL getarg(2,buffer)
   file(2)=TRIM(buffer)

   CALL getarg(3,buffer)
   opcion=TRIM(buffer)

END IF

IF (opcion/="ca" .and. opcion/="backbone" .and. opcion/="all") THEN
   print*,"Error en la variable opcion = ",TRIM(opcion)
   PRINT*,"hint: ./rmsd.exe -help"
   STOP
END IF


!#########################
num_atom=0
ios=0

! miro el numero de atomos
DO j=1,2
   CALL saco_num_atom(opcion, file(j), num_atom(j), ios)
   IF (ios /= 0) THEN
      PRINT*, 'Error de lectura en ', trim(file(j))
      STOP
   END IF
END DO

IF (num_atom(1)/=num_atom(2)) THEN
     PRINT*,'Error:'
     PRINT*,'El numero de atomos de ',trim(file(1)),' (',num_atom(1),')',' y  ',trim(file(2)),' (',num_atom(2),')',&
          ' no coinciden.'
     STOP
END IF

!leo las coordenadas:
N=num_atom(1)
allocate(struct_ref(N,3), struct_2(N,3),g(N,3))

call leo_pdb(opcion, file(1), struct_ref, N)
call leo_pdb(opcion, file(2), struct_2, N)


call min_rmsd(N,struct_ref,struct_2,rot,center_ref,center_2,rmsd,g)


DEALLOCATE(struct_2)
opcion2='all'
CALL saco_num_atom(opcion2, file(2), N, ios)
ALLOCATE(struct_2(N,3))
call leo_pdb(opcion2, file(2), struct_2, N)

DO i=1,N
   struct_2(i,:)=matmul(transpose(rot(:,:)),struct_2(i,:)-center_2)+center_ref
END DO

print'(A,A,A,A,A)', 'HEADER    ','Ajuste de la estructura ',TRIM(file(2)),' sobre ',TRIM(file(1))
print'(A,A,A,A)', 'HEADER    ','minimizando el RMSD de los atomos: ',TRIM(opcion),'.'
print'(A,A,E,A)', 'HEADER    ','(valor del ajuste RMSD =',rmsd,' )'

call print_pdb (file(2),struct_2)



END program main
