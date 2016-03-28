program main

use para_rmsd
use pdb_in_out
use para_dihedros

implicit none

integer::N,i,j,ii,jj,h
double precision,dimension(:,:),allocatable::struct_ref,struct_2,posnew
double precision,dimension(3,3)::rot
double precision,dimension(3)::center_ref,center_2
double precision::rmsd
double precision,dimension(:,:),allocatable::g

character(len=long_max),dimension(1)::file
character(len=long_max)::argumento1,opcion,opcion2
character(len=long_max)::buffer
integer::num_args

integer::ios
integer,dimension(2)::num_atom,num_backbone

integer::num_angs,Nb,num_phi,num_psi,num_res
integer,dimension(:,:),allocatable::indice,para_ordenar
integer,dimension(:),allocatable::ind_backbone,res,new_ind_backbone,ind_all
real,dimension(:),allocatable::angulo
real,dimension(:),allocatable::phi,psi
CHARACTER(len=10),dimension(:),allocatable::label_backbone,label_all

real,dimension(4,3)::pos_aux
real::val_aux,val_aux2

!#############

num_args=IARG()

IF (num_args==1) THEN
   
   CALL getarg(1,buffer)
   argumento1=TRIM(buffer)
   
   IF ((argumento1=="-help").OR.(argumento1=="-h")) THEN
      PRINT*," "
      PRINT*,"move_dihedral estructura.pdb num_angs atom1 atom2 atom3 atom4 ang atom1 ..."
      PRINT*," "
      PRINT*,"El programa gira el o los angulos dihedros elegidos y saca la nueva estructura"
      PRINT*,"fitteada a la original con el backbone"
      PRINT*," "
      PRINT*," "
      PRINT*,"   estructura.pdb:  nombre del fichero pdb (caracteres)"
      PRINT*,"   num_angs:  cantidad de dihedros a cambiar (entero)"
      PRINT*,"   atom1:  indice del primer atomo que define el dihedro (entero)"
      PRINT*,"   atom2:  indice del segundo atomo que define el dihedro (entero)"
      PRINT*,"   atom?:  ... (entero)"
      PRINT*,"   ang:  valor del angulo dihedro deseado en radianes (real)"
      PRINT*," "
      PRINT*," "
      PRINT*,"*: todos los atomos tienen el mismo peso en el analisis del rmsd. "
      PRINT*," "
      STOP
   ELSE
      PRINT*,"ERROR en el numero de argumentos"
      PRINT*,"tip: ./move_dihedral -help"
      STOP
   END IF
END IF

IF (num_args<7) THEN
   PRINT*,"ERROR en el numero de argumentos"
   PRINT*,"hint: ./move_dihedral -help"
   STOP
ELSE

   CALL getarg(1,buffer)
   file(1)=TRIM(buffer)
   
   CALL getarg(2,buffer)
   READ(buffer,*) num_angs

   IF (num_args/=(2+5*num_angs)) THEN
      PRINT*,"ERROR en el numero de argumentos"
      PRINT*,"hint: ./move_dihedral -help"
      STOP
   END IF

   ALLOCATE(indice(num_angs,4),angulo(num_angs))

   indice=0
   angulo=0.0d0

   ii=2
   DO i=1,num_angs
      DO j=1,4
         ii=ii+1
         CALL getarg(ii,buffer)
         READ(buffer,*) indice(i,j)
      END DO
      ii=ii+1
      CALL getarg(ii,buffer)
      READ(buffer,*) angulo(i)
   END DO


   print'(A,A,A)','HEADER     ',TRIM(file(1)),' con los angulos modificados:'
   DO i=1,num_angs
      print 102,'HEADER    ',indice(i,1),indice(i,2),indice(i,3),indice(i,4),' con valor:',angulo(i)
   END DO



END IF


!#########################
num_atom=0
ios=0

!! miro el numero de atomos
DO j=1,1
   opcion="all"
   CALL saco_num_atom(opcion, file(j), num_atom(j), ios)
   IF (ios /= 0) THEN
      PRINT*, 'Error de lectura en ', trim(file(j))
      STOP
   END IF
   opcion="backbone2"
   CALL saco_num_atom(opcion, file(j), num_backbone(j), ios)
   IF (ios /= 0) THEN
      PRINT*, 'Error de lectura en ', trim(file(j))
      STOP
   END IF
END DO

num_atom(2)=num_atom(1)
num_backbone(2)=num_backbone(1)


!! leo coordenadas de todo
N=num_atom(1)
allocate(struct_ref(N,3), struct_2(N,3),res(N))
struct_ref=0.0d0
struct_2=0.0d0
res=0
opcion="all"
call leo_pdb(opcion, file(1), struct_ref, N)
call residuos_pdb(opcion, file(1), res, N)
IF (res(1)/=1) THEN
   num_res=res(1)-1
   DO i=1,N
      res(i)=res(i)-num_res
   END DO
END IF

num_res=res(N)
!! leo indice de los atomos que hacen el backbone
Nb=num_backbone(1)
allocate(ind_backbone(Nb), label_backbone(Nb), label_all(N), ind_all(N))
ind_backbone=0
opcion="backbone2"
call indices_pdb(opcion, file(1), ind_backbone, label_backbone, Nb)
ind_all=0
opcion="all"
call indices_pdb(opcion, file(1), ind_all, label_all, N)



!! Calculo los dihedros de la estructura original

num_phi=0
num_psi=0

DO i=1,Nb-3
   h=0
   IF (label_backbone(i)=='C') THEN
      h=h+1
      ii=res(ind_backbone(i))
      DO j=1,Nb
         IF (label_backbone(j)=='N'.and.res(ind_backbone(j))==(ii+1)) THEN
            h=h+1
         END IF
         IF (label_backbone(j)=='CA'.and.res(ind_backbone(j))==(ii+1)) THEN
            h=h+1
         END IF         
         IF (label_backbone(j)=='C'.and.res(ind_backbone(j))==(ii+1)) THEN
            h=h+1
         END IF
      END DO
      IF (h==4) THEN
         num_phi=num_phi+1
      END IF
   END IF
   IF (label_backbone(i)=='N') THEN
      h=h+1
      ii=res(ind_backbone(i))
      DO j=1,Nb
         IF (label_backbone(j)=='CA'.and.res(ind_backbone(j))==(ii)) THEN
            h=h+1
         END IF
         IF (label_backbone(j)=='C'.and.res(ind_backbone(j))==(ii)) THEN
            h=h+1
         END IF         
         IF (label_backbone(j)=='N'.and.res(ind_backbone(j))==(ii+1)) THEN
            h=h+1
         END IF
      END DO
      IF (h==4) THEN
         num_psi=num_psi+1
      END IF
   END IF   
END DO


allocate(phi(num_phi),psi(num_psi))
phi=0.0d0
psi=0.0d0

num_phi=0
num_psi=0

DO i=1,Nb-3
   h=0
   IF (label_backbone(i)=='C') THEN
      pos_aux=0.0d0
      val_aux=0.0d0
      pos_aux(1,:)=real(struct_ref(ind_backbone(i),:))
      h=h+1
      ii=res(ind_backbone(i))
      DO j=1,Nb
         IF (label_backbone(j)=='N'.and.res(ind_backbone(j))==(ii+1)) THEN
            pos_aux(2,:)=real(struct_ref(ind_backbone(j),:))
            h=h+1
         END IF
         IF (label_backbone(j)=='CA'.and.res(ind_backbone(j))==(ii+1)) THEN
            pos_aux(3,:)=real(struct_ref(ind_backbone(j),:))
            h=h+1
         END IF         
         IF (label_backbone(j)=='C'.and.res(ind_backbone(j))==(ii+1)) THEN
            pos_aux(4,:)=real(struct_ref(ind_backbone(j),:))
            h=h+1
         END IF
      END DO
      IF (h==4) THEN
         num_phi=num_phi+1
         call calc_dihedro(pos_aux,val_aux)
         phi(num_phi)=val_aux
      END IF
   END IF
   IF (label_backbone(i)=='N') THEN
      pos_aux=0.0d0
      val_aux=0.0d0
      pos_aux(1,:)=real(struct_ref(ind_backbone(i),:))
      h=h+1
      ii=res(ind_backbone(i))
      DO j=1,Nb
         IF (label_backbone(j)=='CA'.and.res(ind_backbone(j))==(ii)) THEN
            pos_aux(2,:)=real(struct_ref(ind_backbone(j),:))
            h=h+1
         END IF
         IF (label_backbone(j)=='C'.and.res(ind_backbone(j))==(ii)) THEN
            pos_aux(3,:)=real(struct_ref(ind_backbone(j),:))
            h=h+1
         END IF         
         IF (label_backbone(j)=='N'.and.res(ind_backbone(j))==(ii+1)) THEN
            pos_aux(4,:)=real(struct_ref(ind_backbone(j),:))
            h=h+1
         END IF
      END DO
      IF (h==4) THEN
         num_psi=num_psi+1
         call calc_dihedro(pos_aux,val_aux)
         psi(num_psi)=val_aux
      END IF
   END IF   
END DO


DO i=1,num_psi
   print 101,'HEADER   ','phi1 ',phi(i),'(',degrees(phi(i)),'°)','psi1 ',psi(i),'(',degrees(psi(i)),'°)'
END DO




!! Buen orden del backbone

ALLOCATE(para_ordenar(num_res,4),new_ind_backbone(Nb))
para_ordenar=0
new_ind_backbone=0
DO i=1,Nb
   h=0
   h=para_ordenar(res(ind_backbone(i)),1)+1
   para_ordenar(res(ind_backbone(i)),1)=h
   para_ordenar(res(ind_backbone(i)),h+1)=i
END DO


h=0

DO j=2,para_ordenar(1,1)+1
   IF (label_backbone(para_ordenar(1,j))=='N') THEN
      h=h+1
      new_ind_backbone(h)=ind_backbone(para_ordenar(1,j))
      exit
   END IF
END DO
DO j=2,para_ordenar(1,1)+1
   IF (label_backbone(para_ordenar(1,j))=='CA') THEN
      h=h+1
      new_ind_backbone(h)=ind_backbone(para_ordenar(1,j))
      exit
   END IF
END DO
DO j=2,para_ordenar(1,1)+1
   IF (label_backbone(para_ordenar(1,j))=='CH3') THEN
      h=h+1
      new_ind_backbone(h)=ind_backbone(para_ordenar(1,j))
      exit
   END IF
END DO
DO j=2,para_ordenar(1,1)+1
   IF (label_backbone(para_ordenar(1,j))=='C') THEN
      h=h+1
      new_ind_backbone(h)=ind_backbone(para_ordenar(1,j))
      exit
   END IF
END DO

DO i=2,num_res-1
   DO j=2,para_ordenar(i,1)+1
      IF (label_backbone(para_ordenar(i,j))=='N') THEN
         h=h+1
         new_ind_backbone(h)=ind_backbone(para_ordenar(i,j))
         exit
      END IF
   END DO
   DO j=2,para_ordenar(i,1)+1
      IF (label_backbone(para_ordenar(i,j))=='CA') THEN
         h=h+1
         new_ind_backbone(h)=ind_backbone(para_ordenar(i,j))
         exit
      END IF
   END DO
   DO j=2,para_ordenar(i,1)+1
      IF (label_backbone(para_ordenar(i,j))=='C') THEN
         h=h+1
         new_ind_backbone(h)=ind_backbone(para_ordenar(i,j))
         exit
      END IF
   END DO
END DO


DO j=2,para_ordenar(num_res,1)+1
   IF (label_backbone(para_ordenar(num_res,j))=='N') THEN
      h=h+1
      new_ind_backbone(h)=ind_backbone(para_ordenar(num_res,j))
      exit
   END IF
END DO
DO j=2,para_ordenar(num_res,1)+1
   IF (label_backbone(para_ordenar(num_res,j))=='CA') THEN
      h=h+1
      new_ind_backbone(h)=ind_backbone(para_ordenar(num_res,j))
      exit
   END IF
END DO
DO j=2,para_ordenar(num_res,1)+1
   IF (label_backbone(para_ordenar(num_res,j))=='C') THEN
      h=h+1
      new_ind_backbone(h)=ind_backbone(para_ordenar(num_res,j))
      exit
   END IF
END DO
DO j=2,para_ordenar(num_res,1)+1
   IF (label_backbone(para_ordenar(num_res,j))=='CH3') THEN
      h=h+1
      new_ind_backbone(h)=ind_backbone(para_ordenar(num_res,j))
      exit
   END IF
END DO      

ind_backbone=new_ind_backbone



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! Reconstruyo la estructura con los nuevos dihedros

struct_2=struct_ref

call cambio_dihedro_d(num_angs,indice,angulo,Nb,ind_backbone,N,struct_2,res,label_all)

allocate(g(N,3))
call min_rmsd(N,struct_ref,struct_2,rot,center_ref,center_2,rmsd,g)


DO i=1,N
   struct_2(i,:)=matmul(transpose(rot(:,:)),struct_2(i,:)-center_2)+center_ref
END DO




print'(A)', 'HEADER    '
print'(A,A,A,A,A)', 'HEADER    ','Nueva configuracion con dihedros elegidos ',TRIM(file(1))
print'(A)', 'HEADER    '
print'(A)', 'HEADER    '
print'(A)', 'HEADER    '

call print_pdb (file(1),struct_2)


101 format (A,A,F8.4,x,A1,F9.4,A3,5x,A,F8.4,x,A1,F9.4,A3)
102 format (A,4I4,A,F8.4)






END program main
