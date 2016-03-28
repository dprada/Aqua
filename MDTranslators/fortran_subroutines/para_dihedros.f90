MODULE para_dihedros

use para_rmsd

IMPLICIT NONE

CONTAINS

SUBROUTINE calc_dihedro (posiciones,angulo)

  IMPLICIT NONE

  REAL, DIMENSION(4,3),INTENT(IN)::posiciones
  REAL::angulo
  REAL,DIMENSION(3)::pos1,pos2,pos3,pos4


  angulo=0.0d0

  pos1=0.0d0
  pos2=0.0d0
  pos3=0.0d0
  pos4=0.0d0

  pos1=posiciones(1,:)
  pos2=posiciones(2,:)
  pos3=posiciones(3,:)
  pos4=posiciones(4,:)

  
  CALL calculo_dihed (pos1,pos2,pos3,pos4,angulo)


END SUBROUTINE calc_dihedro

SUBROUTINE multiplico (a,b,c)

  IMPLICIT NONE
  
  real,dimension(3),intent(in)::a,b
  real,dimension(3),intent(out)::c
  
  c(1)=a(2)*b(3)-a(3)*b(2)
  c(2)=-a(1)*b(3)+a(3)*b(1)
  c(3)=a(1)*b(2)-a(2)*b(1)
  
end SUBROUTINE multiplico

SUBROUTINE multiplico_d (a,b,c)

  IMPLICIT NONE
  
  double precision,dimension(3),intent(in)::a,b
  double precision,dimension(3),intent(out)::c
  
  c(1)=a(2)*b(3)-a(3)*b(2)
  c(2)=-a(1)*b(3)+a(3)*b(1)
  c(3)=a(1)*b(2)-a(2)*b(1)
  
end SUBROUTINE multiplico_d

SUBROUTINE calculo_dihed (atom1,atom2,atom3,atom4,angulo)

  IMPLICIT NONE
  
  real,intent(out)::angulo
  real,dimension(3),intent(in)::atom1,atom2,atom3,atom4
  real,dimension(3)::vec1,vec2,vec3,aux1,aux2,aux3
  real::cosa
  integer::signo
  double precision::pi

  vec1=0.0d0
  vec2=0.0d0
  vec3=0.0d0
  cosa=0.0d0
  signo=0
  angulo=0.0d0
  pi=3.14159265358979
  
  vec1(:)=atom2(:)-atom1(:)
  vec2(:)=atom3(:)-atom2(:)
  vec3(:)=atom4(:)-atom3(:)
  
  aux1=0.0d0
  aux2=0.0d0
  
  CALL multiplico(vec1,vec2,aux1)
  CALL multiplico(vec2,vec3,aux2)
  
  cosa=(aux1(1)*aux2(1)+aux1(2)*aux2(2)+aux1(3)*aux2(3))/(sqrt(dot_product(aux1,aux1))*sqrt(dot_product(aux2,aux2)))

  IF (cosa>=1.0d0) THEN 
     cosa=1.0d0
  END IF
  IF (cosa<=-1.0d0) THEN 
     cosa=-1.0d0
  END IF


  cosa=(acos(cosa))!*(180.0d0/pi)

 ! IF (cosa>180.0d0) THEN
 !    print*,'ERROR EN ANGULOS'
 !    STOP
 ! END IF
  
  aux3=0.0d0
  
  CALL multiplico(aux1,aux2,aux3)
  
  IF ( dot_product(aux3,vec2) <= 0.0d0 ) THEN
     signo=-1
  ELSE
     signo=+1
  END IF
  
  angulo=cosa*signo

END SUBROUTINE calculo_dihed


function degrees(angulo)

  implicit none
  double precision,parameter::pi=3.14159265358979
  real,INTENT(IN)::angulo
  real::degrees

  degrees=real(angulo/pi)*180

end function degrees
  
function radians(angulo)

  implicit none
  double precision,parameter::pi=3.14159265358979
  real,INTENT(IN)::angulo
  real::radians

  radians=real(angulo/180.0d0)*real(pi)

end function radians
  



subroutine cambio_dihedro_d(num_angs,indice,angulo,Nb,ind_backbone,N,struct_2,res,label_all)

IMPLICIT NONE

INTEGER,INTENT(IN)::num_angs,N,Nb
REAL,dimension(num_angs),INTENT(IN)::angulo
INTEGER,DIMENSION(N),INTENT(IN)::res
INTEGER,DIMENSION(num_angs,4)::indice
INTEGER,DIMENSION(Nb),INTENT(IN)::ind_backbone
DOUBLE PRECISION,DIMENSION(N,3),INTENT(INOUT)::struct_2
REAL,DIMENSION(N,3)::estructura,nueva
LOGICAL,DIMENSION(num_angs)::filtro
LOGICAL,DIMENSION(N)::tachado
LOGICAL::switch
REAL::pshi
CHARACTER(len=10),dimension(N),INTENT(IN)::label_all
INTEGER::i,j,contador,num_res,Naux,ii,h
REAL,DIMENSION(3)::pos1,pos2,pos3,pos4,aux,aux2
REAL::cosa
REAL,DIMENSION(Nb-1)::dist
REAL,DIMENSION(Nb-2)::atres

double precision,dimension(3,3)::rot
double precision,dimension(3)::center_ref,center_2
double precision::rmsd
double precision,dimension(:,:),allocatable::g
double precision,dimension(:,:,:),allocatable::rota
double precision,dimension(:,:),allocatable::centro_ref,centro_2
double precision,DIMENSION(3)::daux
double precision,dimension(:,:),allocatable::struct_ref,struct_vieja

filtro=.true.
switch=.false.

estructura=0.0d0
estructura=real(struct_2)
nueva=0.0d0
tachado=.false.
num_res=res(N)


!! Primero pongo el backbone

dist=0.0d0
atres=0.0d0

DO i=1,Nb-1
   aux=0.0d0
   aux(:)=estructura(ind_backbone(i),:)-estructura(ind_backbone(i+1),:)
   
   dist(i)=sqrt(dot_product(aux(:),aux(:)))
END DO

DO i=1,Nb-2
   aux=0.0d0
   aux2=0.0d0
   cosa=0.0d0
   aux(:)=estructura(ind_backbone(i),:)-estructura(ind_backbone(i+1),:)
   aux2(:)=estructura(ind_backbone(i+1),:)-estructura(ind_backbone(i+2),:)
   cosa=dot_product(aux,aux2)/(sqrt(dot_product(aux,aux))*sqrt(dot_product(aux2,aux2)))
   atres(i)=acos(cosa)
END DO



DO i=1,3
   nueva(ind_backbone(i),:)=estructura(ind_backbone(i),:)
   tachado(ind_backbone(i))=.true.
END DO

DO i=4,Nb

   pshi=0.0d0

   switch=.false.
   DO j=1,num_angs
      IF (filtro(j)==.true.) THEN
         IF (ind_backbone(i-3)==indice(j,1)) THEN
            IF (ind_backbone(i-2)==indice(j,2)) THEN
               IF (ind_backbone(i-1)==indice(j,3)) THEN
                  IF (ind_backbone(i)==indice(j,4)) THEN
                     filtro(j)=.false.
                     switch=.true.
                     pshi=angulo(j)
                     exit
                  ENDIF
               END IF
            END IF
         END IF
      END IF
   END DO

   IF (switch==.false.) THEN

      pos1=estructura(ind_backbone(i-3),:)
      pos2=estructura(ind_backbone(i-2),:)
      pos3=estructura(ind_backbone(i-1),:)
      pos4=estructura(ind_backbone(i),:)
      CALL calculo_dihed (pos1,pos2,pos3,pos4,pshi)
      

   END IF

   pos1=nueva(ind_backbone(i-3),:)
   pos2=nueva(ind_backbone(i-2),:)
   pos3=nueva(ind_backbone(i-1),:)
   pos4=0.0d0

   call remonto(pos1,pos2,pos3,pshi,dist(i-1),atres(i-2),pos4)
   CALL calculo_dihed (pos1,pos2,pos3,pos4,pshi)


   nueva(ind_backbone(i),:)=pos4
   tachado(ind_backbone(i))=.true.

END DO

!! coloco los O o H de los grupos peptidicos:

Naux=4
allocate(struct_ref(4,3),struct_vieja(4,3),g(Naux,3))
struct_ref=0.0d0
struct_vieja=0.0d0

DO i=1,N
   IF (label_all(i)=='O') THEN
      ii=res(i)
      DO j=1,N
         IF (label_all(j)=='H'.and.res(j)==(ii+1)) THEN
            tachado(i)=.true.
            tachado(j)=.true.

            rot=0.0d0
            center_ref=0.0d0
            center_2=0.0d0

            DO h=1,N
               IF ((label_all(h)=='CA'.or.label_all(h)=='CH3').and.(res(h)==ii)) THEN
                  struct_ref(1,:)=dble(nueva(h,:))
                  struct_vieja(1,:)=dble(estructura(h,:))
               END IF
               IF ((label_all(h)=='C').and.(res(h)==ii)) THEN
                  struct_ref(2,:)=dble(nueva(h,:))
                  struct_vieja(2,:)=dble(estructura(h,:))
               END IF
               IF ((label_all(h)=='N').and.(res(h)==(ii+1))) THEN
                  struct_ref(3,:)=dble(nueva(h,:))
                  struct_vieja(3,:)=dble(estructura(h,:))
               END IF
               IF ((label_all(h)=='CA'.or.label_all(h)=='CH3').and.(res(h)==(ii+1))) THEN
                  struct_ref(4,:)=dble(nueva(h,:))
                  struct_vieja(4,:)=dble(estructura(h,:))
               END IF
            END DO

            
            call min_rmsd(Naux,struct_ref,struct_vieja,rot,center_ref,center_2,rmsd,g)


            IF (rmsd>0.001) THEN
               print*,'ERROR 00 EN SUBRUTINA cambio_dihedro_d DEL MODULO para_dihedros'
               stop
            END IF

            daux=0.0d0
            daux=dble(estructura(i,:))
            daux(:)=matmul(transpose(rot(:,:)),daux(:)-center_2)+center_ref
            nueva(i,:)=real(daux(:))

            daux=0.0d0
            daux=dble(estructura(j,:))
            daux(:)=matmul(transpose(rot(:,:)),daux(:)-center_2)+center_ref
            nueva(j,:)=real(daux(:))


         END IF
      END DO
   END IF
END DO

deallocate(struct_ref,struct_vieja,g)


!! coloco el resto

!! coloco lo del primer residuo (fije su backbone)

DO i=1,N
   IF (res(i)==1.and.tachado(i)==.false.) THEN
      nueva(i,:)=estructura(i,:)
   END IF
END DO


allocate(struct_ref(3,3),struct_vieja(3,3),g(3,3))
allocate(rota(num_res-1,3,3),centro_ref(num_res-1,3),centro_2(num_res-1,3))
rota=0.0d0
centro_ref=0.0d0
centro_2=0.0d0
rot=0.0d0
center_ref=0.0d0
center_2=0.0d0



Naux=3
DO i=2,num_res-1
   contador=0
   DO j=1,Nb
      IF (res(ind_backbone(j))==i) THEN
         contador=contador+1
         struct_ref(contador,:)=dble(nueva(ind_backbone(j),:))
         struct_vieja(contador,:)=dble(estructura(ind_backbone(j),:))
      END IF
   END DO


   call min_rmsd(Naux,struct_ref,struct_vieja,rot,center_ref,center_2,rmsd,g)

   rota(i-1,:,:)=rot(:,:)
   centro_ref(i-1,:)=center_ref
   centro_2(i-1,:)=center_2

   IF (rmsd>0.001) THEN
      
      print*,'ERROR 01 EN SUBRUTINA cambio_dihedro_d DEL MODULO para_dihedros'
      stop
      
   END IF
END DO

h=0
DO i=1,Nb
   IF (res(ind_backbone(i))==num_res) THEN
      h=h+1
   END IF
END DO

IF (h==3) THEN
   contador=0
   DO j=1,Nb
      IF (res(ind_backbone(j))==num_res) THEN
         contador=contador+1
         struct_ref(contador,:)=dble(nueva(ind_backbone(j),:))
         struct_vieja(contador,:)=dble(estructura(ind_backbone(j),:))
      END IF
   END DO


   call min_rmsd(Naux,struct_ref,struct_vieja,rot,center_ref,center_2,rmsd,g)

   rota(num_res-1,:,:)=rot(:,:)
   centro_ref(num_res-1,:)=center_ref
   centro_2(num_res-1,:)=center_2

   IF (rmsd>0.001) THEN
      
      print*,'ERROR 02 EN SUBRUTINA cambio_dihedro_d DEL MODULO para_dihedros'
      stop
      
   END IF
END IF
IF (h==2) THEN
   DO j=1,Nb
      IF (res(ind_backbone(j))==(num_res-1)) THEN
         contador=1
         struct_ref(contador,:)=dble(nueva(ind_backbone(j),:))
         struct_vieja(contador,:)=dble(estructura(ind_backbone(j),:)) !!es redundante, pero solo queda el ultimo
      END IF
   END DO
   DO j=1,Nb
      IF (res(ind_backbone(j))==num_res) THEN
         contador=contador+1
         struct_ref(contador,:)=dble(nueva(ind_backbone(j),:))
         struct_vieja(contador,:)=dble(estructura(ind_backbone(j),:))
      END IF
   END DO

   call min_rmsd(Naux,struct_ref,struct_vieja,rot,center_ref,center_2,rmsd,g)

   rota(num_res-1,:,:)=rot(:,:)
   centro_ref(num_res-1,:)=center_ref
   centro_2(num_res-1,:)=center_2

   IF (rmsd>0.001) THEN
      
      print*,'ERROR 02 EN SUBRUTINA cambio_dihedro_d DEL MODULO para_dihedros'
      stop
      
   END IF
END IF



DO i=1,N
   IF (tachado(i)==.false.) THEN
      
      IF (res(i)==1) THEN
         nueva(i,:)=estructura(i,:)
         tachado(i)=.true.
      END IF

      IF ((res(i)/=1)) THEN

         daux=0.0d0
         daux=dble(estructura(i,:))
         rot=0.0d0
         center_2=0.0d0
         center_ref=0.0d0
         rot=rota(res(i)-1,:,:)
         center_ref=centro_ref(res(i)-1,:)
         center_2=centro_2(res(i)-1,:)

         daux(:)=matmul(transpose(rot(:,:)),daux(:)-center_2)+center_ref
         nueva(i,:)=real(daux(:))
         
         tachado(i)=.true.
      END IF

   END IF
END DO
         
struct_2=dble(nueva)





END subroutine cambio_dihedro_d



SUBROUTINE remonto(atom1,atom2,atom3,dihedro,distancia,angulo,atom4)

implicit none

real,dimension(3),intent(in)::atom1,atom2,atom3
real,dimension(3),intent(out)::atom4
real::dihedro,distancia,angulo
double precision::mod
double precision,dimension(3)::vec1,vec2,aux,aux2,aux3

atom4=0.0d0
vec1=0.0d0
vec2=0.0d0
aux=0.0d0
aux2=0.0d0
aux3=0.0d0
mod=0.d0


vec1(:)=dble(atom2(:))-dble(atom1(:))
vec2(:)=dble(atom3(:))-dble(atom2(:))

mod=dsqrt(dot_product(vec2(:),vec2(:)))
aux(:)=vec2(:)/mod



aux(:)=dble(distancia)*dcos(dble(angulo))*aux(:)

CALL multiplico_d(vec1,vec2,aux2)
mod=0.0d0
mod=dsqrt(dot_product(aux2(:),aux2(:)))
aux2(:)=aux2(:)/mod

CALL multiplico_d(aux2,vec2,aux3)

aux2(:)=dble(distancia)*dsin(dble(angulo))*dsin(dble(dihedro))*aux2(:)

mod=0.0d0
mod=dsqrt(dot_product(aux3(:),aux3(:)))
aux3(:)=aux3(:)/mod


aux3(:)=dble(distancia)*dsin(dble(angulo))*dcos(dble(dihedro))*aux3(:)

aux(:)=aux(:)+aux2(:)+aux3(:)

mod=0.0d0
mod=dsqrt(dot_product(aux(:),aux(:)))


atom4(:)=atom3(:)+real(aux(:))


END SUBROUTINE remonto




END MODULE para_dihedros
