MODULE GLOB

INTEGER::num_frames,num_parts,dim_mss
INTEGER,DIMENSION(:),ALLOCATABLE::lim_inf_mss,lim_sup_mss
INTEGER,DIMENSION(:,:,:),ALLOCATABLE::traj_mss
INTEGER,DIMENSION(:,:),ALLOCATABLE::labels

CONTAINS

SUBROUTINE INIT_TRAJ_MSS_2_NET ()

  IF (ALLOCATED(traj_mss)) DEALLOCATE(traj_mss)

  ALLOCATE(traj_mss(num_frames,num_parts,dim_mss))

  !traj_mss=0
  !print*,'SIZE', sizeof(traj_mss)/1073741824.0

END SUBROUTINE INIT_TRAJ_MSS_2_NET

SUBROUTINE TRAJ_MSS_2_TRAJ_NODES ()

  !! You have to have right to write on the directory
  
  INTEGER,DIMENSION(:,:),ALLOCATABLE:: traj_nodes
  LOGICAL,DIMENSION(:,:),ALLOCATABLE::occup
  INTEGER,DIMENSION(:,:),ALLOCATABLE::indice
  INTEGER,DIMENSION(:),ALLOCATABLE::deantes
  INTEGER:: nnn,i,j,box,counted,a,b
  CHARACTER*20::f1

  ALLOCATE(traj_nodes(nw,num_frames))

  DO box=1,dim_mss

     IF (box==1) THEN
        counted=1
        traj_nodes=1
     ELSE
        CALL SYSTEM ('mv prov_trad_aux.aux prov_trad_aux_old.aux')
     END IF
     
     ALLOCATE(occup(counted,lim_inf_mss(box):lim_sup_mss(box)),indice(counted,lim_inf_mss(box):lim_sup_mss(box)))
     occup=.false.
   
     DO i=1,num_frames
        DO ii=1,num_parts
           b=traj_nodes(ii,i)
           a=traj_mss(ii,i,box)
           occup(b,a)=.true.
        END DO
     END DO
   
     nnn=0
     WRITE(f1,*) box   !! Aqui habia un '(I)' donde el asterisco
     f1="(I,"//TRIM(ADJUSTL(f1))//"I3)"

     OPEN(21,FILE="prov_trad_aux.aux",status="REPLACE",ACTION="WRITE")
     IF (box==1) THEN
        DO i=1,counted
           DO j=lim_inf_mss(box),lim_sup_mss(box)
              IF (occup(i,j).eqv..true.) THEN
                 nnn=nnn+1
                 indice(i,j)=nnn
                 WRITE(21,f1) nnn,j
              END IF
           END DO
        END DO
     ELSE
        OPEN(61,FILE="prov_trad_aux_old.aux",status="OLD",ACTION="READ")
        ALLOCATE(deantes(box-1))         
        DO i=1,counted
           READ(61,*) a,deantes(:)
           DO j=lim_inf_mss(box),lim_sup_mss(box)
              IF (occup(i,j).eqv..true.) THEN
                 nnn=nnn+1
                 WRITE(21,f1) nnn,deantes(:),j
                 indice(i,j)=nnn
              END IF
           END DO
        END DO
        DEALLOCATE(deantes)
        CLOSE(61)
     END IF
     CLOSE(21)

     DO i=1,num_frames
        DO ii=1,num_parts
           b=traj_nodes(ii,i)
           a=traj_mss(ii,i,box)
           traj_nodes(ii,i)=indice(b,a)
        END DO
     END DO
     DEALLOCATE(indice,occup)
     counted=nnn
     
     IF (box/=1) CALL SYSTEM('rm prov_trad_aux_old.aux')
     
  END DO

  DEALLOCATE(traj_mss)
  ALLOCATE(labels(counted,dim_mss))
  
  OPEN(61,FILE="prov_trad_aux.aux",status="OLD",ACTION="READ")

  DO i=1,counted
     READ(61) a,labels(i,:)
  END DO

  CLOSE(61)
  CALL SYSTEM('rm prov_trad_aux_old.aux')


END SUBROUTINE TRAJ_MSS_2_TRAJ_NODES



SUBROUTINE EVOLUTION_STEP(T_start,T_ind,T_tau,vect_in,N_nodes,Ktot,vect_out)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::N_nodes,Ktot
  INTEGER,DIMENSION(Ktot),INTENT(IN)::T_ind
  DOUBLE PRECISION,DIMENSION(Ktot),INTENT(IN)::T_tau
  INTEGER,DIMENSION(N_nodes+1),INTENT(IN)::T_start
  
  DOUBLE PRECISION,DIMENSION(N_nodes),INTENT(IN)::vect_in
  DOUBLE PRECISION,DIMENSION(N_nodes),INTENT(OUT)::vect_out
  DOUBLE PRECISION,DIMENSION(N_nodes)::Pe

  INTEGER::i,j,ii,jj
  
  Pe=0
  vect_out=0.0d0

  DO i=1,N_nodes
     DO j=T_start(i)+1,T_start(i+1)
        Pe(i)=Pe(i)+T_tau(j)
     END DO
  END DO

  DO i=1,N_nodes
     DO j=T_start(i)+1,T_start(i+1)
        jj=T_ind(j)
        vect_out(jj)=vect_out(jj)+(T_tau(j)/Pe(i))*vect_in(i)
     END DO
  END DO

END SUBROUTINE EVOLUTION_STEP

SUBROUTINE BROWNIAN_RUN (with_loops,T_start,T_ind,T_tau,iseed,node_alpha,length,N_nodes,Ktot,vect_out)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::N_nodes,Ktot,length,with_loops
  INTEGER,DIMENSION(Ktot),INTENT(IN)::T_ind
  DOUBLE PRECISION,DIMENSION(Ktot),INTENT(IN)::T_tau
  INTEGER,DIMENSION(N_nodes+1),INTENT(IN)::T_start
  INTEGER,DIMENSION(4),INTENT(IN)::iseed
  INTEGER,INTENT(IN)::node_alpha
  INTEGER,DIMENSION(length+1),INTENT(OUT)::vect_out

  DOUBLE PRECISION,DIMENSION(N_nodes)::Pe
  INTEGER,DIMENSION(4)::hseed
  INTEGER::place
  DOUBLE PRECISION::dice,bandera

  INTEGER::i,j

  hseed=iseed

  IF (with_loops==1) THEN

     Pe=0
     DO i=1,N_nodes
        DO j=T_start(i)+1,T_start(i+1)
           Pe(i)=Pe(i)+T_tau(j)
        END DO
     END DO
     
     vect_out=0
     place=node_alpha+1
     vect_out(1)=0
     
     DO j=2,length+1
        
        CALL dlarnv(1,hseed,1,dice)
        dice=dice*Pe(place)
        
        bandera=0.0d0
        DO i=T_start(place)+1,T_start(place+1)
           bandera=bandera+T_tau(i)
           IF (dice<=bandera) THEN
              place=T_ind(i)
              exit
           END IF
        END DO
        
        vect_out(j)=place-1
        
     END DO

  ELSE

     Pe=0
     DO i=1,N_nodes
        DO j=T_start(i)+1,T_start(i+1)
           IF (T_ind(j)/=i) THEN
              Pe(i)=Pe(i)+T_tau(j)
           END IF
        END DO
     END DO
     
     vect_out=0
     place=node_alpha+1
     vect_out(1)=0
     
     DO j=2,length+1
        
        CALL dlarnv(1,hseed,1,dice)
        dice=dice*Pe(place)
        
        bandera=0.0d0
        DO i=T_start(place)+1,T_start(place+1)
           IF (T_ind(i)/=place) THEN
              bandera=bandera+T_tau(i)
              IF (dice<=bandera) THEN
                 place=T_ind(i)
                 exit
              END IF
           END IF
        END DO
        
        vect_out(j)=place-1
        
     END DO

  END IF


END SUBROUTINE BROWNIAN_RUN

SUBROUTINE BROWNIAN_RUN_FPT (T_start,T_ind,T_tau,iseed,node_alpha,node_omega,N_nodes,Ktot,steps)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::N_nodes,Ktot
  INTEGER,DIMENSION(Ktot),INTENT(IN)::T_ind
  DOUBLE PRECISION,DIMENSION(Ktot),INTENT(IN)::T_tau
  INTEGER,DIMENSION(N_nodes+1),INTENT(IN)::T_start
  INTEGER,DIMENSION(4),INTENT(IN)::iseed
  INTEGER,INTENT(IN)::node_alpha,node_omega
  INTEGER,INTENT(OUT)::steps

  DOUBLE PRECISION,DIMENSION(N_nodes)::Pe
  INTEGER,DIMENSION(4)::hseed
  INTEGER::place,omega
  DOUBLE PRECISION::dice,bandera

  INTEGER::i,j

  hseed=iseed

  Pe=0
  DO i=1,N_nodes
     DO j=T_start(i)+1,T_start(i+1)
        Pe(i)=Pe(i)+T_tau(j)
     END DO
  END DO

  place=node_alpha+1
  omega=node_omega+1
  
  steps=0
  DO WHILE (place/=omega)
     steps=steps+1
     CALL dlarnv(1,hseed,1,dice)
     dice=dice*Pe(place)
     
     bandera=0.0d0
     DO i=T_start(place)+1,T_start(place+1)
        bandera=bandera+T_tau(i)
        IF (dice<=bandera) THEN
           place=T_ind(i)
           exit
        END IF
     END DO
     
  END DO

END SUBROUTINE BROWNIAN_RUN_FPT


SUBROUTINE DETAILED_BALANCE_DISTANCE(db_dist,p,T_start,T_ind,T_tau,N_nodes,Ktot)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::N_nodes,Ktot
  INTEGER,DIMENSION(Ktot),INTENT(IN)::T_ind
  DOUBLE PRECISION,DIMENSION(Ktot),INTENT(IN)::T_tau
  INTEGER,DIMENSION(N_nodes+1),INTENT(IN)::T_start
  DOUBLE PRECISION,INTENT(IN)::p
  DOUBLE PRECISION,INTENT(OUT)::db_dist


  INTEGER::i,j,ii,jj
  DOUBLE PRECISION::overp,aux_val
  DOUBLE PRECISION::We_total,i_weight,i_trans,j_weight,j_trans

  overp=1.0d0/p
  db_dist=0.0d0

  We_total=0

  DO i=1,N_nodes
     DO j=T_start(i)+1,T_start(i+1)
        We_total=We_total+T_tau(j)
     END DO
  END DO

  DO i=1,N_nodes
     DO ii=T_start(i)+1,T_start(i+1)
        j=T_ind(ii)
        i_trans=T_tau(ii)
        j_trans=0
        DO jj=T_start(j)+1,T_start(j+1)
           IF (i==T_ind(jj)) THEN
              j_trans=T_tau(jj)
              exit
           END IF
        END DO
        db_dist=db_dist+(abs((i_trans-j_trans))/(We_total))**p
     END DO
  END DO

  db_dist=(0.50d0*db_dist)**overp


END SUBROUTINE DETAILED_BALANCE_DISTANCE

SUBROUTINE SYMMETRIZE_NET(newKmax,TT_tau,TT_ind,TT_start,Pe,newKtot,T_ind,T_tau,T_start,N_nodes,Ktot)

  IMPLICIT NONE

  TYPE iarray_pointer
     INTEGER,DIMENSION(:),POINTER::p1
  END TYPE iarray_pointer
  TYPE darray_pointer
     DOUBLE PRECISION,DIMENSION(:),POINTER::d1
  END TYPE darray_pointer


  INTEGER,INTENT(IN)::N_nodes,Ktot,newKtot
  INTEGER,DIMENSION(Ktot),INTENT(IN)::T_ind
  DOUBLE PRECISION,DIMENSION(Ktot),INTENT(IN)::T_tau
  INTEGER,DIMENSION(N_nodes+1),INTENT(IN)::T_start

  DOUBLE PRECISION,DIMENSION(N_nodes),INTENT(OUT)::Pe
  INTEGER,DIMENSION(N_nodes+1),INTENT(OUT)::TT_start
  INTEGER,DIMENSION(newKtot),INTENT(OUT)::TT_ind
  DOUBLE PRECISION,DIMENSION(newKtot),INTENT(OUT)::TT_tau
  INTEGER,INTENT(OUT)::newKmax

  LOGICAL::interruptor
  INTEGER::i,j,l,h,g,gg,destino
  DOUBLE PRECISION::aux,aux2
  integer,dimension(:),allocatable::salidas,auxx1
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::SL,auxx2
  TYPE(iarray_pointer),DIMENSION(:),POINTER::F_ind
  TYPE(darray_pointer),DIMENSION(:),POINTER::flux
  
  Pe=0.0d0
  TT_start=0
  TT_ind=0
  TT_tau=0.0d0
  newKmax=0


  ALLOCATE(F_ind(N_nodes),flux(N_nodes),salidas(N_nodes),SL(N_nodes))
  salidas=0
  SL=0.0d0

  DO i=1,N_nodes

     DO j=T_start(i)+1,T_start(i+1)

        IF (T_ind(j)==i) THEN
           SL(i)=T_tau(j)*2
        ELSE
           
           destino=T_ind(j)

           gg=salidas(i)
           IF (gg==0) THEN
              ALLOCATE(F_ind(i)%p1(1),flux(i)%d1(1))
              F_ind(i)%p1(1)=destino
              flux(i)%d1(1)=T_tau(j)
              salidas(i)=1
           ELSE
              interruptor=.false.
              DO h=1,gg
                 IF (F_ind(i)%p1(h)==destino) THEN
                    flux(i)%d1(h)=flux(i)%d1(h)+T_tau(j)
                    interruptor=.true.
                    exit
                 END IF
              END DO
              IF (interruptor.eqv..false.) THEN
                 ALLOCATE(auxx1(gg+1),auxx2(gg+1))
                 auxx1(1:gg)=F_ind(i)%p1(:)
                 auxx2(1:gg)=flux(i)%d1(:)
                 auxx1(gg+1)=destino
                 auxx2(gg+1)=T_tau(j)
                 salidas(i)=gg+1
                 DEALLOCATE(F_ind(i)%p1,flux(i)%d1)
                 ALLOCATE(F_ind(i)%p1(gg+1),flux(i)%d1(gg+1))
                 F_ind(i)%p1(:)=auxx1(:)
                 flux(i)%d1(:)=auxx2(:)
                 DEALLOCATE(auxx1,auxx2)
              END IF
           END IF

           gg=salidas(destino)   
           IF (gg==0) THEN
              ALLOCATE(F_ind(destino)%p1(1),flux(destino)%d1(1))
              F_ind(destino)%p1(1)=i
              flux(destino)%d1(1)=T_tau(j)
              salidas(destino)=1
           ELSE
              interruptor=.false.
              DO h=1,gg
                 IF (F_ind(destino)%p1(h)==i) THEN
                    flux(destino)%d1(h)=flux(destino)%d1(h)+T_tau(j)
                    interruptor=.true.
                    exit
                 END IF
              END DO
              IF (interruptor.eqv..false.) THEN
                 ALLOCATE(auxx1(gg+1),auxx2(gg+1))
                 auxx1(1:gg)=F_ind(destino)%p1(:)
                 auxx2(1:gg)=flux(destino)%d1(:)
                 auxx1(gg+1)=i
                 auxx2(gg+1)=T_tau(j)
                 salidas(destino)=gg+1
                 DEALLOCATE(F_ind(destino)%p1,flux(destino)%d1)
                 ALLOCATE(F_ind(destino)%p1(gg+1),flux(destino)%d1(gg+1))
                 F_ind(destino)%p1(:)=auxx1(:)
                 flux(destino)%d1(:)=auxx2(:)
                 DEALLOCATE(auxx1,auxx2)
              END IF
           END IF
           
        END IF

     END DO

  END DO

  g=0
  DO i=1,N_nodes
     IF (SL(i)>0) THEN
        salidas(i)=salidas(i)+1
     END IF
  END DO
  g=SUM(salidas(:),DIM=1)



  TT_ind=0
  TT_tau=0
  TT_start=0
  Pe=0

  g=0
  DO i=1,N_nodes
     TT_start(i)=g
     gg=g
     g=g+salidas(i)
     salidas(i)=0
     aux=SL(i)
     IF (aux>0.0d0) THEN
        TT_ind(gg+1)=i
        TT_tau(gg+1)=aux
        Pe(i)=aux
        salidas(i)=1
     END IF
  END DO
  TT_start(N_nodes+1)=g

  DO i=1,N_nodes
     IF ((TT_start(i+1)-TT_start(i)-salidas(i))>0) THEN
        gg=size(F_ind(i)%p1(:),DIM=1)
        aux2=0.0d0
        DO j=1,gg
           g=salidas(i)+1
           salidas(i)=g
           g=g+TT_start(i)
           aux=flux(i)%d1(j)
           TT_ind(g)=F_ind(i)%p1(j)
           TT_tau(g)=aux
           aux2=aux2+aux
        END DO
        Pe(i)=Pe(i)+aux2
     END IF
  END DO

  newKmax=MAXVAL(salidas(:),DIM=1)


  DEALLOCATE(F_ind,flux,salidas,SL)

  Pe=Pe/2.0d0
  TT_tau=TT_tau/2.0d0


END SUBROUTINE SYMMETRIZE_NET



SUBROUTINE GRAD (N_sets,comunidades,T_ind,T_tau,T_start,N_nodes,Ktot)   !!!! Para revisar!!!! T_tau debe ser double precision


  IMPLICIT NONE
  

  INTEGER,INTENT(IN)::N_nodes,Ktot
  INTEGER,DIMENSION(Ktot),INTENT(IN)::T_ind,T_tau
  INTEGER,DIMENSION(N_nodes+1),INTENT(IN)::T_start


  INTEGER,INTENT(OUT)::N_sets
  INTEGER,DIMENSION(N_nodes),INTENT(OUT)::comunidades

  INTEGER,DIMENSION(:),ALLOCATABLE::Pe


  INTEGER:: i,j,jj,g,h,dim

  !! Para minimos:
  INTEGER::dim2,inicio,fin,peso,candidato
  INTEGER,DIMENSION(:),ALLOCATABLE::tope_lista
  INTEGER,DIMENSION(:,:),ALLOCATABLE::lista
  LOGICAL,DIMENSION(:),ALLOCATABLE::label_glob_min,label_loc_min,filtro
  logical::inter

  !! Para las basins:
  INTEGER::hacia,desde
  logical,dimension(:),allocatable::label
  integer,dimension(:),allocatable::vect_aux1

  ALLOCATE(Pe(N_nodes))
  Pe=0

  DO i=1,N_nodes
     DO j=T_start(i)+1,T_start(i+1)
        Pe(i)=Pe(i)+T_tau(j)
     END DO
  END DO


  dim=1
  dim2=dim

  !! Ordeno primero los contactos

  ALLOCATE(lista(N_nodes,dim2),tope_lista(N_nodes))

  lista=0
  tope_lista=0

  DO i=1,N_nodes

     h=0
     g=T_start(i+1)-T_start(i)

     ALLOCATE(filtro(g))
     filtro=.true.

     inicio=T_start(i)+1
     fin=T_start(i+1)
     
     DO j=inicio,fin
        IF (T_ind(j)==i) THEN
           filtro(j-inicio+1)=.false.
           EXIT
        END IF
     END DO

     DO WHILE (h<dim2)
        
        g=inicio-1+MAXLOC(T_tau(inicio:fin),DIM=1,MASK=filtro)
        candidato=g
        peso=Pe(T_ind(candidato))
        DO j=inicio,fin
           
           IF ((T_tau(j)==T_tau(candidato)).and.(filtro(j-inicio+1).eqv..true.)) THEN
              IF (Pe(T_ind(j))>peso) THEN
                 peso=Pe(T_ind(j))
                 candidato=j
              END IF
           END IF
        END DO
        
        h=h+1
        lista(i,h)=T_ind(candidato)
        filtro(candidato-inicio+1)=.false.
        
        IF (COUNT(filtro)==0) THEN
           exit
        END IF
        
     END DO

     DEALLOCATE(filtro)
     tope_lista(i)=h

  END DO


  ALLOCATE(label_glob_min(N_nodes),label_loc_min(N_nodes))
  label_glob_min=.false.
  label_loc_min=.false.


  !! Minimos globales:

  label_glob_min=.false.

  DO i=1,N_nodes

     peso=Pe(i)
     inter=.true.
     DO j=T_start(i)+1,T_start(i+1)
        g=T_ind(j)
        IF (peso<Pe(g)) THEN
           inter=.false.
           exit
        END IF
     END DO
     
     IF (inter.eqv..true.) label_glob_min(i)=.true.
     
  END DO
  

  !! Minimos locales:

  label_loc_min=.false.
  
  DO i=1,N_nodes
     
     peso=Pe(i)
     
     inter=.true.
     DO j=1,tope_lista(i)
        
        IF (peso<=Pe(lista(i,j))) THEN
           inter=.false.
           exit
        END IF
        
     END DO
     
     IF (inter.eqv..true.) label_loc_min(i)=.true.
     
  END DO


  !! Defino las basins

  ALLOCATE(label(N_nodes))
  label=.false.
  label=label_loc_min
  
  allocate(vect_aux1(N_nodes))
  vect_aux1=0
    
  DO i=1,N_nodes
     IF (label_loc_min(i).eqv..true.) THEN
        comunidades(i)=i
     END IF
  END DO


  DO i=1,N_nodes

     desde=i
     h=0
     
     DO WHILE (label(desde).eqv..false.)

        h=h+1
        vect_aux1(h)=desde
        label(desde)=.true.
        
        DO j=1,tope_lista(desde)
           IF (Pe(lista(desde,j))>=Pe(desde)) THEN
              hacia=lista(desde,j)
              inter=.true.
              exit
           END IF
        END DO

        IF (inter.eqv..false.) THEN
           
           exit
        END IF
        
        desde=hacia
        
     END DO
     

     
     IF (comunidades(desde)==0) THEN
        j=desde
     ELSE
        j=comunidades(desde)
     END IF
     
     DO jj=1,h
        comunidades(vect_aux1(jj))=j
     END DO
     
  END DO


  label=.false.
  DO i=1,N_nodes
     label(comunidades(i))=.true.
  END DO
  h=0

  DO i=1,N_nodes
     IF (label(i).eqv..true.) THEN
        h=h+1
     END IF
  END DO
  

  comunidades=comunidades-1

  N_sets=h


  DEALLOCATE(Pe)
  DEALLOCATE(tope_lista,lista,label_glob_min,label_loc_min)
  DEALLOCATE(label,vect_aux1)

END SUBROUTINE GRAD


SUBROUTINE BUILD_NET_BIN (f_traj_nodes,f_net,nw,frames)   !!!! Para revisar el double precision de T_tau

IMPLICIT NONE

TYPE array_pointer
   INTEGER,DIMENSION(:),POINTER::p1
END TYPE array_pointer

CHARACTER(80),INTENT(IN)::f_traj_nodes,f_net
INTEGER,INTENT(IN)::nw,frames

INTEGER::i,j,ii,g,h,kk,l,hh
INTEGER,DIMENSION(:,:),ALLOCATABLE::tray


!from the shared variables:

INTEGER::Ktot,Kmax
INTEGER,DIMENSION(:),ALLOCATABLE::T_start,T_ind,T_tau,Pe


!para la red
INTEGER::N_nodes,contador,x,x2,index,K_out_i
INTEGER,DIMENSION(:),ALLOCATABLE::SL,W,K_out
TYPE(array_pointer),DIMENSION(:),POINTER::C,WK_out
INTEGER,DIMENSION(:),POINTER::aux_puntero,aux_puntero2
LOGICAL::switch
INTEGER,DIMENSION(:),ALLOCATABLE::auxiliar

ALLOCATE(tray(nw,frames))
tray=0

OPEN(unit=21,FILE=TRIM(f_traj_nodes),STATUS='old',action='READ',FORM='UNFORMATTED')

N_nodes=0
DO i=1,frames
   READ (21) ii, tray(:,i)
   j=MAXVAL(tray(:,i))
   IF (j>N_nodes) N_nodes=j
END DO


ALLOCATE(C(N_nodes),K_out(N_nodes))
ALLOCATE(WK_out(N_nodes),SL(N_nodes),W(N_nodes))
ALLOCATE(aux_puntero(25000),aux_puntero2(25000))

DO i=1,N_nodes
   ALLOCATE(C(i)%p1(1),WK_out(i)%p1(1))
   C(i)%p1(:)=0
   WK_out(i)%p1(:)=0
END DO

SL=0
W=0
K_out=0
kk=0
x=0
x2=0
aux_puntero(:)=0
aux_puntero2(:)=0
contador=0

DO index=1,nw

   x=tray(index,1)
   contador=contador+1
   x2=tray(index,2)
   contador=contador+1

   DO i=3,frames

      W(x)=W(x)+1
      
      switch=.true.
      DO g=1,K_out(x)
         IF(C(x)%p1(g)==x2) THEN
            WK_out(x)%p1(g)=WK_out(x)%p1(g)+1
            switch=.false.
            EXIT
         END IF
      END DO

      IF (switch.eqv..true.) THEN
         IF (K_out(x)>0) THEN
            g=K_out(x)
            IF (g>25000) THEN
               print*,'tenemos un problema'
               stop
            END IF
            aux_puntero(1:g)=C(x)%p1(:)
            aux_puntero2(1:g)=WK_out(x)%p1(:)
            DEALLOCATE(C(x)%p1,WK_out(x)%p1)
            ALLOCATE(C(x)%p1(g+1),WK_out(x)%p1(g+1))
            C(x)%p1(1:g)=aux_puntero(:)
            WK_out(x)%p1(1:g)=aux_puntero2(:)
            C(x)%p1(g+1)=x2
            WK_out(x)%p1(g+1)=1
            K_out(x)=K_out(x)+1
         ELSE
            WK_out(x)%p1(1)=1
            C(x)%p1(1)=x2
            K_out(x)=1
         END IF
      END IF

      x=x2  !!bin
      
      x2=tray(index,i)
      contador=contador+1

   END DO

   W(x)=W(x)+1
   
   switch=.true.
   
   DO g=1,K_out(x)
      IF(C(x)%p1(g)==x2) THEN
         WK_out(x)%p1(g)=WK_out(x)%p1(g)+1
         switch=.false.
         EXIT
      END IF
   END DO
   
   IF (switch.eqv..true.) THEN
      IF (K_out(x)>0) THEN
         g=K_out(x)
         IF (g>25000) THEN
            print*,'tenemos un problema'
            stop
         END IF
         aux_puntero(1:g)=C(x)%p1(:)
         aux_puntero2(1:g)=WK_out(x)%p1(:)
         DEALLOCATE(C(x)%p1,WK_out(x)%p1)
         ALLOCATE(C(x)%p1(g+1),WK_out(x)%p1(g+1))
         C(x)%p1(1:g)=aux_puntero(:)
         WK_out(x)%p1(1:g)=aux_puntero2(:)
         C(x)%p1(g+1)=x2
         WK_out(x)%p1(g+1)=1
         K_out(x)=K_out(x)+1
      ELSE
         WK_out(x)%p1(1)=1
         C(x)%p1(1)=x2
         K_out(x)=1
      END IF
   END IF
 
   x=x2
   contador=contador+1

   W(x)=W(x)+1

END DO


!!Saco a un fichero para comparar:

Kmax=maxval(K_out(:))
Ktot=sum(K_out(:))

ALLOCATE(T_start(N_nodes+1),T_ind(Ktot),T_tau(Ktot),Pe(N_nodes))
T_start=0
T_ind=0
T_tau=0
Pe=0


l=0
DO i=1,N_nodes

   h=0
   switch=.false.
   DO g=1,K_out(i)
      IF (C(i)%p1(g)==i) THEN
         switch=.true.
         exit
      END IF
   END DO

   T_start(i)=l
   IF (switch.eqv..true.) THEN
      l=l+1
      T_ind(l)=i
      hh=WK_out(i)%p1(g)
      T_tau(l)=hh
      h=h+hh
      DO j=1,K_out(i)
         IF (j/=g) THEN
            l=l+1
            T_ind(l)=C(i)%p1(j)
            hh=WK_out(i)%p1(j)
            T_tau(l)=hh
            h=h+hh
         END IF
      END DO
   ELSE
      DO j=1,K_out(i)
         l=l+1
         T_ind(l)=C(i)%p1(j)
         hh=WK_out(i)%p1(j)
         h=h+hh
         T_tau(l)=hh
      END DO
   END IF

   Pe(i)=h

END DO
T_start(N_nodes+1)=l

!Sacar a fichero los datos para mi:

open(UNIT=22,file=TRIM(f_net),action='write',status='new')

WRITE (22,*) N_nodes,Kmax,Ktot

DO i=1,N_nodes

   K_out_i=T_start(i+1)-T_start(i)
   ALLOCATE(auxiliar(2*K_out_i))
   auxiliar=0
   
   g=T_start(i)
   
   DO j=1,K_out_i
      auxiliar(2*j-1)=T_ind(g+j)
      auxiliar(2*j)=T_tau(g+j)
   END DO
   
   WRITE (22,*) i,K_out_i,Pe(i),auxiliar(:)
   
   DEALLOCATE(auxiliar)
   
END DO

CLOSE(22)


PRINT*,'# Number of frames analysed:', nw*frames
PRINT*, '# Total weight Links out:',SUM(Pe(:),DIM=1)
PRINT*, '# Total connectivity', Ktot
PRINT*, '# Max. connectivity',Kmax

DEALLOCATE(C,K_out)
DEALLOCATE(WK_out,SL,W)
DEALLOCATE(aux_puntero,aux_puntero2)
DEALLOCATE(tray)



END SUBROUTINE BUILD_NET_BIN



SUBROUTINE BUILD_NET (f_traj_nodes,f_net,nw,frames)   !!!! Para revisar el double precision de T_tau

IMPLICIT NONE

TYPE array_pointer
   INTEGER,DIMENSION(:),POINTER::p1
END TYPE array_pointer

CHARACTER(80),INTENT(IN)::f_traj_nodes,f_net
INTEGER,INTENT(IN)::nw,frames

INTEGER::i,j,ii,g,h,kk,l,hh
INTEGER,DIMENSION(:,:),ALLOCATABLE::tray


!from the shared variables:

INTEGER::Ktot,Kmax
INTEGER,DIMENSION(:),ALLOCATABLE::T_start,T_ind,T_tau,Pe


!para la red
INTEGER::N_nodes,contador,x,x2,index,K_out_i
INTEGER,DIMENSION(:),ALLOCATABLE::SL,W,K_out
TYPE(array_pointer),DIMENSION(:),POINTER::C,WK_out
INTEGER,DIMENSION(:),POINTER::aux_puntero,aux_puntero2
LOGICAL::switch
INTEGER,DIMENSION(:),ALLOCATABLE::auxiliar

ALLOCATE(tray(nw,frames))
tray=0

OPEN(unit=21,FILE=TRIM(f_traj_nodes),STATUS='old',action='READ')

N_nodes=0
DO i=1,frames
   READ (21,*) ii
   DO j=1,nw
      READ (21,*) g
      print*,g
      tray(j,i)=g
      IF (g>N_nodes) N_nodes=g
   END DO
END DO
CLOSE(21)

ALLOCATE(C(N_nodes),K_out(N_nodes))
ALLOCATE(WK_out(N_nodes),SL(N_nodes),W(N_nodes))
ALLOCATE(aux_puntero(25000),aux_puntero2(25000))

DO i=1,N_nodes
   ALLOCATE(C(i)%p1(1),WK_out(i)%p1(1))
   C(i)%p1(:)=0
   WK_out(i)%p1(:)=0
END DO

SL=0
W=0
K_out=0
kk=0
x=0
x2=0
aux_puntero(:)=0
aux_puntero2(:)=0
contador=0

DO index=1,nw

   x=tray(index,1)
   contador=contador+1
   x2=tray(index,2)
   contador=contador+1

   DO i=3,frames

      W(x)=W(x)+1
      
      switch=.true.
      DO g=1,K_out(x)
         IF(C(x)%p1(g)==x2) THEN
            WK_out(x)%p1(g)=WK_out(x)%p1(g)+1
            switch=.false.
            EXIT
         END IF
      END DO

      IF (switch.eqv..true.) THEN
         IF (K_out(x)>0) THEN
            g=K_out(x)
            IF (g>25000) THEN
               print*,'tenemos un problema'
               stop
            END IF
            aux_puntero(1:g)=C(x)%p1(:)
            aux_puntero2(1:g)=WK_out(x)%p1(:)
            DEALLOCATE(C(x)%p1,WK_out(x)%p1)
            ALLOCATE(C(x)%p1(g+1),WK_out(x)%p1(g+1))
            C(x)%p1(1:g)=aux_puntero(:)
            WK_out(x)%p1(1:g)=aux_puntero2(:)
            C(x)%p1(g+1)=x2
            WK_out(x)%p1(g+1)=1
            K_out(x)=K_out(x)+1
         ELSE
            WK_out(x)%p1(1)=1
            C(x)%p1(1)=x2
            K_out(x)=1
         END IF
      END IF

      x=x2  !!bin
      
      x2=tray(index,i)
      contador=contador+1

   END DO

   W(x)=W(x)+1
   
   switch=.true.
   
   DO g=1,K_out(x)
      IF(C(x)%p1(g)==x2) THEN
         WK_out(x)%p1(g)=WK_out(x)%p1(g)+1
         switch=.false.
         EXIT
      END IF
   END DO
   
   IF (switch.eqv..true.) THEN
      IF (K_out(x)>0) THEN
         g=K_out(x)
         IF (g>25000) THEN
            print*,'tenemos un problema'
            stop
         END IF
         aux_puntero(1:g)=C(x)%p1(:)
         aux_puntero2(1:g)=WK_out(x)%p1(:)
         DEALLOCATE(C(x)%p1,WK_out(x)%p1)
         ALLOCATE(C(x)%p1(g+1),WK_out(x)%p1(g+1))
         C(x)%p1(1:g)=aux_puntero(:)
         WK_out(x)%p1(1:g)=aux_puntero2(:)
         C(x)%p1(g+1)=x2
         WK_out(x)%p1(g+1)=1
         K_out(x)=K_out(x)+1
      ELSE
         WK_out(x)%p1(1)=1
         C(x)%p1(1)=x2
         K_out(x)=1
      END IF
   END IF
 
   x=x2
   contador=contador+1

   W(x)=W(x)+1

END DO


!!Saco a un fichero para comparar:

Kmax=maxval(K_out(:))
Ktot=sum(K_out(:))

ALLOCATE(T_start(N_nodes+1),T_ind(Ktot),T_tau(Ktot),Pe(N_nodes))
T_start=0
T_ind=0
T_tau=0
Pe=0


l=0
DO i=1,N_nodes

   h=0
   switch=.false.
   DO g=1,K_out(i)
      IF (C(i)%p1(g)==i) THEN
         switch=.true.
         exit
      END IF
   END DO

   T_start(i)=l
   IF (switch.eqv..true.) THEN
      l=l+1
      T_ind(l)=i
      hh=WK_out(i)%p1(g)
      T_tau(l)=hh
      h=h+hh
      DO j=1,K_out(i)
         IF (j/=g) THEN
            l=l+1
            T_ind(l)=C(i)%p1(j)
            hh=WK_out(i)%p1(j)
            T_tau(l)=hh
            h=h+hh
         END IF
      END DO
   ELSE
      DO j=1,K_out(i)
         l=l+1
         T_ind(l)=C(i)%p1(j)
         hh=WK_out(i)%p1(j)
         h=h+hh
         T_tau(l)=hh
      END DO
   END IF

   Pe(i)=h

END DO
T_start(N_nodes+1)=l

!Sacar a fichero los datos para mi:

open(UNIT=22,file=TRIM(f_net),action='write',status='new')

WRITE (22,*) N_nodes,Kmax,Ktot

DO i=1,N_nodes

   K_out_i=T_start(i+1)-T_start(i)
   ALLOCATE(auxiliar(2*K_out_i))
   auxiliar=0
   
   g=T_start(i)
   
   DO j=1,K_out_i
      auxiliar(2*j-1)=T_ind(g+j)
      auxiliar(2*j)=T_tau(g+j)
   END DO
   
   WRITE (22,*) i,K_out_i,Pe(i),auxiliar(:)
   
   DEALLOCATE(auxiliar)
   
END DO

CLOSE(22)


PRINT*,'# Number of frames analysed:', nw*frames
PRINT*, '# Total weight Links out:',SUM(Pe(:),DIM=1)
PRINT*, '# Total connectivity', Ktot
PRINT*, '# Max. connectivity',Kmax

DEALLOCATE(C,K_out)
DEALLOCATE(WK_out,SL,W)
DEALLOCATE(aux_puntero,aux_puntero2)
DEALLOCATE(tray)



END SUBROUTINE BUILD_NET


SUBROUTINE cfep_pfold (plot,info_1,info_2,A,B,T_ind,T_tau,T_start,length,num_iter,N_nodes,Ktot)

  implicit none

  INTEGER,INTENT(IN)::A,B,N_nodes,Ktot,length,num_iter
  INTEGER,DIMENSION(Ktot),INTENT(IN)::T_ind
  DOUBLE PRECISION,DIMENSION(Ktot),INTENT(IN)::T_tau
  INTEGER,DIMENSION(N_nodes+1),INTENT(IN)::T_start

  REAL,DIMENSION(length,3),INTENT(OUT)::plot
  REAL,DIMENSION(N_nodes,2),INTENT(OUT)::info_2
  INTEGER,DIMENSION(N_nodes),INTENT(OUT)::info_1

  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::Pe
  
  integer::i,j,times,jj,g

  double precision,dimension(:),allocatable::Pf,Pf2
  double precision,dimension(:,:),allocatable::appear
  
  logical,dimension(:),allocatable::filtro,filtro2
  double precision::cut
  double precision::delta_cut
  double precision::Z,Za,Zab


  plot=0.0d0
  info_2=0.0d0
  info_1=0


  ALLOCATE(Pe(N_nodes))
  Pe=0.0d0
  DO i=1,N_nodes
     DO j=T_start(i),T_start(i+1)
        Pe(i)=Pe(i)+T_tau(j)
     END DO
  END DO

  ALLOCATE(Pf(N_nodes),Pf2(N_nodes),appear(N_nodes,2))
  Pf=0.0d0
  Pf2=0.0d0
  appear=0.0d0

  Pf2(A)=1.0d0
  Pf2(B)=0.0d0


  DO times=1,num_iter

     Pf2(A)=1.0d0
     Pf2(B)=0.0d0
     Pf=Pf2
     Pf2=0.0d0
     DO i=1,N_nodes
        DO j=T_start(i)+1,T_start(i+1)
           jj=T_ind(j)
           Pf2(i)=Pf2(i)+(T_tau(j)*Pf(jj))/Pe(i)
        END DO
     END DO

  END DO

  Pf2(A)=1.0d0
  Pf2(B)=0.0d0
  Pf=Pf2
  
  DEALLOCATE(Pf2)
  ALLOCATE(filtro(N_nodes),filtro2(N_nodes))

  filtro2=.false.

  Z=sum(Pe(:),dim=1)

  delta_cut=1.0d0/(length*1.0d0)

  DO g=1,length

     cut=delta_cut*g

     Za=0
     Zab=0
     filtro=.false.

     DO i=1,N_nodes
        IF (Pf(i)<cut) THEN 
           filtro(i)=.true.
        END IF
     END DO

     DO i=1,N_nodes
        IF (filtro(i).eqv..true.) THEN
           Za=Za+Pe(i)
           DO j=T_start(i)+1,T_start(i+1)
              jj=T_ind(j)
              IF (filtro(jj).eqv..false.) THEN
                 Zab=Zab+T_tau(j)
              END IF
           END DO
        END IF
     END DO

     
     !WRITE(154,*) ((Za*1.0d0)/(Z*1.0d0)),-log((Zab*1.0d0)/(Z*1.0d0)),cut
     plot(g,1)=(Za/Z)
     plot(g,2)=-log(Zab/Z)
     plot(g,3)=cut

     DO i=1,N_nodes
        IF (filtro(i).eqv..true.) THEN 
           IF (filtro2(i).eqv..false.) THEN
              filtro2(i)=.true.
              appear(i,1)=Za/Z
              appear(i,2)=-log(Zab/Z)
           END IF
        END IF
     END DO     

  END DO

  filtro=.true.
  DO i=1,N_nodes
     g=minloc(appear(:,1),DIM=1,MASK=filtro(:))
     !WRITE(155,*) g, appear(g,1), appear(g,2)
     info_1(i)=g
     info_2(i,1)=appear(g,1)
     info_2(i,2)=appear(g,2)
     filtro(g)=.false.
  END DO
!  PRINT*,COUNT(filtro)

END SUBROUTINE cfep_pfold



SUBROUTINE cfep_mfpt (plot,info_1,info_2,A,T_ind,T_tau,T_start,length,num_iter,N_nodes,Ktot)

  implicit none

  INTEGER,INTENT(IN)::A,N_nodes,Ktot,length,num_iter
  INTEGER,DIMENSION(Ktot),INTENT(IN)::T_ind
  DOUBLE PRECISION,DIMENSION(Ktot),INTENT(IN)::T_tau
  INTEGER,DIMENSION(N_nodes+1),INTENT(IN)::T_start

  REAL,DIMENSION(length,3),INTENT(OUT)::plot
  REAL,DIMENSION(N_nodes,2),INTENT(OUT)::info_2
  INTEGER,DIMENSION(N_nodes),INTENT(OUT)::info_1

  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::Pe
  
  integer::i,j,times,jj,g

  double precision,dimension(:),allocatable::Pf,Pf2,antes
  double precision,dimension(:,:),allocatable::appear
  
  logical,dimension(:),allocatable::filtro,filtro2
  double precision::cut
  double precision::delta_cut,mm
  double precision::Z,Za,Zab
  logical::interr

  plot=0.0d0
  info_2=0.0d0
  info_1=0


  ALLOCATE(Pe(N_nodes))
  Pe=0.0d0
  DO i=1,N_nodes
     DO j=T_start(i),T_start(i+1)
        Pe(i)=Pe(i)+T_tau(j)
     END DO
  END DO

  ALLOCATE(Pf(N_nodes),Pf2(N_nodes),appear(N_nodes,2),antes(N_nodes))
  Pf=0.0d0
  Pf2=0.0d0
  appear=0.0d0

  Pf2=100000.0d0
  Pf2(A)=0.0d0


  DO times=1,num_iter

     Pf2(A)=0.0d0
     Pf=Pf2
     antes=Pf2
     Pf2=0.0d0
     
     DO i=1,N_nodes
        DO j=T_start(i)+1,T_start(i+1)
           jj=T_ind(j)
           Pf2(i)=Pf2(i)+(T_tau(j)*Pf(jj))/Pe(i)
        END DO
     END DO

     antes=antes-Pf2
     interr=.false.
     DO i=1,N_nodes
        IF (abs(antes(i))>0.0010d0) THEN
           interr=.true.
           exit
        END IF
     END DO
     
     IF (interr.eqv..false.) exit

     Pf2=Pf2+1.0d0
  END DO



  Pf2(A)=0.0d0

  Pf=Pf2

  
  DEALLOCATE(Pf2)
  ALLOCATE(filtro(N_nodes),filtro2(N_nodes))

  filtro2=.false.

  Z=sum(Pe(:),dim=1)

  mm=maxval(Pf(:),DIM=1)

  delta_cut=(mm*1.0d0)/(length*1.0d0)

  DO g=1,length

     cut=delta_cut*g

     Za=0.0d0
     Zab=0.0d0
     filtro=.false.

     DO i=1,N_nodes
        IF (Pf(i)<cut) THEN 
           filtro(i)=.true.
        END IF
     END DO

     DO i=1,N_nodes
        IF (filtro(i).eqv..true.) THEN
           Za=Za+Pe(i)
           DO j=T_start(i)+1,T_start(i+1)
              jj=T_ind(j)
              IF (filtro(jj).eqv..false.) THEN
                 Zab=Zab+T_tau(j)
              END IF
           END DO
        END IF
     END DO

     
     !WRITE(154,*) ((Za*1.0d0)/(Z*1.0d0)),-log((Zab*1.0d0)/(Z*1.0d0)),cut
     plot(g,1)=(Za/Z)
     plot(g,2)=-log(Zab/Z)
     plot(g,3)=cut

     DO i=1,N_nodes
        IF (filtro(i).eqv..true.) THEN 
           IF (filtro2(i).eqv..false.) THEN
              filtro2(i)=.true.
              appear(i,1)=Za/Z
              appear(i,2)=-log(Zab/Z)
           END IF
        END IF
     END DO     

  END DO

  filtro=.true.
  DO i=1,N_nodes
     g=minloc(appear(:,1),DIM=1,MASK=filtro(:))
     !WRITE(155,*) g, appear(g,1), appear(g,2)
     info_1(i)=g
     info_2(i,1)=appear(g,1)
     info_2(i,2)=appear(g,2)
     filtro(g)=.false.
  END DO
!  PRINT*,COUNT(filtro)

END SUBROUTINE cfep_mfpt



SUBROUTINE COMPONENTS(num_comp,componente,T_start,T_ind,T_tau,N_nodes,Ktot)

  IMPLICIT NONE
  
  INTEGER,INTENT(IN)::N_nodes,Ktot
  INTEGER,DIMENSION(Ktot),INTENT(IN)::T_ind
  DOUBLE PRECISION,DIMENSION(Ktot),INTENT(IN)::T_tau
  INTEGER,DIMENSION(N_nodes+1),INTENT(IN)::T_start
  
  
  INTEGER,DIMENSION(N_nodes),INTENT(OUT)::componente
  INTEGER,INTENT(OUT)::num_comp
  
  INTEGER::i,j,jj,g,gg,h,hh
  
  componente=0
  num_comp=0
  
  DO i=1,N_nodes
     
     IF (componente(i)==0) THEN
        num_comp=num_comp+1
        g=num_comp
        componente(i)=g
     ELSE
        g=componente(i)
     END IF
     
     DO j=T_start(i)+1,T_start(i+1)
        jj=T_ind(j)
        IF (componente(jj)==0) THEN
           componente(jj)=g
        ELSE
           IF (componente(jj)/=g) THEN
              IF (componente(jj)<g) THEN
                 gg=componente(jj)
              ELSE
                 gg=g
                 g=componente(jj)
              END IF
              DO h=1,N_nodes
                 IF (componente(h)==g) THEN
                    componente(h)=gg
                 ELSE
                    IF (componente(h)>g) componente(h)=componente(h)-1
                 END IF
              END DO
              g=gg
              num_comp=num_comp-1
           END IF
        END IF
     END DO
  END DO

  componente=componente-1
  
END SUBROUTINE COMPONENTS


SUBROUTINE DIJKSTRA (distancia,node,dim_out,directed,T_start,T_ind,T_tau,N_nodes,Ktot)
 
  IMPLICIT NONE
  
  INTEGER,INTENT(IN)::N_nodes,Ktot,node,directed,dim_out
  INTEGER,DIMENSION(Ktot),INTENT(IN)::T_ind
  DOUBLE PRECISION,DIMENSION(Ktot),INTENT(IN)::T_tau
  INTEGER,DIMENSION(N_nodes+1),INTENT(IN)::T_start
  DOUBLE PRECISION,DIMENSION(N_nodes,dim_out),INTENT(OUT)::distancia
  
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::vect_aux
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro,filtro2
  
  DOUBLE PRECISION::azero
  INTEGER::i,j,jj,g,gg,h,hh
  
  ALLOCATE(filtro(N_nodes),vect_aux(N_nodes))
  
  distancia=0.0d0
  azero=0.0d0

  IF (node>0) THEN
     
     vect_aux=1.0d0/azero
     filtro=.true.
     vect_aux(node)=0.0d0
     
     DO WHILE (COUNT(filtro)>0)
        i=MINLOC(vect_aux(:),DIM=1,MASK=filtro)
        DO j=T_start(i)+1,T_start(i+1)
           g=T_ind(j)
           IF (filtro(g).eqv..true.) THEN
              IF (vect_aux(g)>(vect_aux(i)+T_tau(j))) THEN
                 vect_aux(g)=vect_aux(i)+T_tau(j)
              END IF
           END IF
        END DO
        filtro(i)=.false.
     END DO
     distancia(:,1)=vect_aux(:)

  ELSE

     ALLOCATE(filtro2(N_nodes))
     DO h=1,N_nodes
        vect_aux=1.0d0/azero
        filtro=.true.
        filtro2=.true.
        vect_aux(h)=0.0d0
        IF (directed==0) THEN
           DO i=1,h-1
              vect_aux(i)=distancia(h,i)
              filtro2(i)=.false.
           END DO
        END IF
        DO WHILE (COUNT(filtro)>0)
           i=MINLOC(vect_aux(:),DIM=1,MASK=filtro)
           DO j=T_start(i)+1,T_start(i+1)
              g=T_ind(j)
              IF (filtro2(g).eqv..true.) THEN
                 IF (vect_aux(g)>(vect_aux(i)+T_tau(j))) THEN
                    vect_aux(g)=vect_aux(i)+T_tau(j)
                 END IF
              END IF
           END DO

           filtro(i)=.false.
           filtro2(i)=.false.
        END DO
        distancia(:,h)=vect_aux(:)
     END DO
     DEALLOCATE(filtro2)

  END IF

  DEALLOCATE(filtro,vect_aux)

END SUBROUTINE DIJKSTRA



SUBROUTINE MDS (coordinates,eigenvals,eigenvects,stress,directed,opt,opt_stress,dim,lout,&
     T_start,T_ind,T_tau,distances,N_nodes,Ktot,dim_distances)
 
  IMPLICIT NONE
  
  INTEGER,INTENT(IN)::N_nodes,Ktot,lout,dim,opt,opt_stress,directed,dim_distances
  INTEGER,DIMENSION(Ktot),INTENT(IN)::T_ind
  DOUBLE PRECISION,DIMENSION(Ktot),INTENT(IN)::T_tau
  INTEGER,DIMENSION(N_nodes+1),INTENT(IN)::T_start
  DOUBLE PRECISION,DIMENSION(dim_distances,dim_distances),INTENT(IN)::distances

  DOUBLE PRECISION,DIMENSION(lout),INTENT(OUT)::eigenvals
  DOUBLE PRECISION,DIMENSION(N_nodes,lout),INTENT(OUT)::eigenvects
  DOUBLE PRECISION,DIMENSION(N_nodes,dim),INTENT(OUT)::coordinates
  DOUBLE PRECISION,DIMENSION(lout),INTENT(OUT)::stress

  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::vect_aux,coors_aux
  DOUBLE PRECISION::dd,norm,salida
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::di
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::CC
  INTEGER::num_val,info,Lwork
  INTEGER, DIMENSION(:), ALLOCATABLE::iwork,ifail
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::work

  INTEGER::h,i,j,g,abajo,ii,jj,gg

  eigenvals=0.0d0
  eigenvects=0.0d0
  coordinates=0.0d0
  stress=0.0d0

  ALLOCATE(CC(N_nodes,N_nodes),di(N_nodes))
  di=0.0d0
  dd=0.0d0
  CC=0.0d0

  IF (dim_distances==N_nodes) THEN
     CC=distances
  ELSE
     !! Dijkstra:
     IF (opt==1) THEN
        CALL DIJKSTRA(CC,-1,N_nodes,directed,T_start,T_ind,T_tau,N_nodes,Ktot)
     ELSE
        DO h=1,N_nodes
           DO i=T_start(h)+1,T_start(h+1)
              j=T_ind(i)
              CC(j,h)=T_tau(i)
           END DO
        END DO
     END IF
  END IF
  print*,'ENTRA 1'

  !! Deberia poner aqui una opcion para hacer un guardado opcional distancias=CC
  !! y usarlas luego como opcion al calcular el stress.

  DO i=1,N_nodes
     DO j=1,N_nodes
        CC(j,i)=CC(j,i)**2    ! Necesario en el metodo
        di(i)=di(i)+CC(j,i)
        dd=dd+CC(j,i)
     END DO
     di(i)=di(i)/(N_nodes*1.0d0)
  END DO
  dd=dd/(N_nodes*N_nodes*1.0d0)

  DO i=1,N_nodes
     DO j=1,N_nodes
        CC(j,i)=-0.50d0*(CC(j,i)-di(i)-di(j)+dd)
     END DO
  END DO

  Lwork=8*N_nodes

  ALLOCATE (work(Lwork),iwork(5*N_nodes),ifail(N_nodes))

  eigenvals=0.0d0
  eigenvects=0.0d0
  work=0.0d0
  iwork=0
  ifail=0
  abajo=N_nodes-lout+1


  print*,'ENTRA 2'

  CALL dsyevx ('V','I','U',N_nodes,CC,N_nodes,0,0,abajo,N_nodes,0.0d0,num_val&
       &,eigenvals,eigenvects,N_nodes,work,Lwork,iwork,ifail,info)

  IF (info/=0) THEN
     print*,"#Error with the diagonalization -MDS.f90:dsyevx-"
     print*,"#the array 'work' should has the dimension:",work(1)
  END IF

  DEALLOCATE (work,iwork,ifail,CC)
  DEALLOCATE(di)

  coordinates=0.0d0

  DO j=1,dim
     coordinates(:,j)=sqrt(eigenvals(lout-j+1))*eigenvects(:,lout-j+1)
  END DO

  IF (opt_stress==1) THEN
     ALLOCATE(vect_aux(Ktot),coors_aux(N_nodes))
     vect_aux=0.0d0
     coors_aux=0.0d0
     norm=0.0d0
     DO ii=1,N_nodes
        DO jj=T_start(ii)+1,T_start(ii+1)
           norm=norm+T_tau(jj)**2
        END DO
     END DO
     DO j=1,lout
        if (eigenvals(lout-j+1)>0.0d0) THEN
           salida=0.0d0
           dd=sqrt(eigenvals(lout-j+1))
           coors_aux(:)=dd*eigenvects(:,lout-j+1)
           DO ii=1,N_nodes
              DO jj=T_start(ii)+1,T_start(ii+1)
                 gg=T_ind(jj)
                 vect_aux(jj)=vect_aux(jj)+(coors_aux(ii)-coors_aux(gg))**2
                 salida=salida+(sqrt(vect_aux(jj))-T_tau(jj))**2
              END DO
           END DO
           stress(j)=sqrt(salida/norm)
        END IF
     END DO
     DEALLOCATE(vect_aux,coors_aux)
  END IF
  

END SUBROUTINE MDS


SUBROUTINE MCL (N_sets,comunidades,granularity,epsilon,iterations,T_start,T_ind,T_tau,N_nodes,Ktot)
 
  IMPLICIT NONE
  
  DOUBLE PRECISION, INTENT(IN):: granularity,epsilon
  INTEGER,INTENT(IN)::N_nodes,Ktot,iterations
  INTEGER,DIMENSION(Ktot),INTENT(IN)::T_ind
  DOUBLE PRECISION,DIMENSION(Ktot),INTENT(IN)::T_tau
  INTEGER,DIMENSION(N_nodes+1),INTENT(IN)::T_start

  INTEGER,INTENT(OUT)::N_sets
  INTEGER,DIMENSION(N_nodes),INTENT(OUT)::comunidades

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: Pe
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: matriz, matriz_seg, matriz2

  INTEGER::i,j,ii,g,h
  DOUBLE PRECISION:: aux
  LOGICAL::inter

  ALLOCATE(Pe(N_nodes),matriz(N_nodes,N_nodes),matriz2(N_nodes,N_nodes))
  Pe=0.0d0
  DO i=1,N_nodes
     DO j=T_start(i)+1, T_start(i+1)
        Pe(i)=Pe(i)+T_tau(j)
     END DO
  END DO

  matriz=0.0d0
  matriz2=0.0d0

  DO i=1,N_nodes
     DO j=T_start(i)+1,T_start(i+1)
        matriz(i,T_ind(j))=T_tau(j)/Pe(i)
        
     END DO
  END DO
  
  inter=.true.
  ii=0

  DO WHILE (inter.eqv..true.)
     ii=ii+1

     matriz2=matmul(matriz,matriz)
     
     DO j=1,N_nodes
        DO g=1,N_nodes
           matriz2(j,g)=matriz2(j,g)**granularity
        END DO
     END DO
     
     DO j=1,N_nodes
        aux=sum(matriz2(j,:))
        matriz2(j,:)=matriz2(j,:)/aux
     END DO
     

     IF (iterations<0.1) THEN
        inter=.false.
        DO g=1,N_nodes
           DO h=1,N_nodes
              IF (abs(matriz(g,h)-matriz2(g,h))>epsilon) THEN
                 inter=.true.
                 exit
              END IF
           END DO
           IF (inter.eqv..true.) THEN
              exit
           END IF
        END DO
     END IF

     IF (ii==iterations) THEN
        inter=.false.
     END IF
     
     matriz=matriz2

  END DO

  N_sets=0
  DO g=1,N_nodes
     IF (sum(matriz(:,g))/=0.0d0) THEN
        N_sets=N_sets+1
     END IF
  END DO
  
  h=0
  DO g=1,N_nodes
     IF (sum(matriz(:,g))/=0.0d0) THEN
        h=h+1
        DO j=1,N_nodes
           IF (matriz(j,g)/=0.0d0) THEN
              comunidades(j)=h
           END IF
        END DO
     END IF
  END DO
  
  DEALLOCATE(Pe,matriz,matriz2)

  comunidades=comunidades-1

END SUBROUTINE MCL

END MODULE GLOB



!! f2py --f90flags='-fast' -c -m pyn_fort_net pyn_fort_net.f90
!!f2py --f90flags='-fast' -c -m pyn_fort_net pyn_fort_net.f90 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_def -lpthread
