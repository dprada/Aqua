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

SUBROUTINE WEIGHT_CORE_NEW_K_TOTAL(newKtot,newNnodes,threshold,T_ind,T_tau,T_start,N_nodes,Ktot)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::N_nodes,Ktot
  DOUBLE PRECISION,INTENT(IN)::threshold
  INTEGER,DIMENSION(Ktot),INTENT(IN)::T_ind
  DOUBLE PRECISION,DIMENSION(Ktot),INTENT(IN)::T_tau
  INTEGER,DIMENSION(N_nodes+1),INTENT(IN)::T_start

  INTEGER,INTENT(OUT)::newKtot,newNnodes

  INTEGER::i,j
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro
  DOUBLE PRECISION::Pe

  ALLOCATE(filtro(N_nodes))

  filtro=.TRUE.

  DO i=1,N_nodes
     Pe=0.0d0
     DO j=T_start(i)+1,T_start(i+1)
        Pe=Pe+T_tau(j)
     END DO
     IF (Pe<=threshold) filtro(i)=.FALSE.
  END DO

  newKtot=0
  DO i=1,N_nodes
     IF (filtro(i)==.TRUE.) THEN
        DO j=T_start(i)+1,T_start(i+1)
           IF (filtro(T_ind(j))==.TRUE.) THEN
              newKtot=newKtot+1
           END IF
        END DO
     END IF
  END DO

  newNnodes=COUNT(filtro,DIM=1)

  DEALLOCATE(filtro)

END SUBROUTINE WEIGHT_CORE_NEW_K_TOTAL

SUBROUTINE EXTRACT_NET_NEW_K_TOTAL(newKtot,list_nodes,T_ind,T_tau,T_start,newNnodes,N_nodes,Ktot)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::N_nodes,Ktot,newNnodes
  INTEGER,DIMENSION(Ktot),INTENT(IN)::T_ind
  DOUBLE PRECISION,DIMENSION(Ktot),INTENT(IN)::T_tau
  INTEGER,DIMENSION(N_nodes+1),INTENT(IN)::T_start
  INTEGER,DIMENSION(newNnodes),INTENT(IN)::list_nodes

  INTEGER,INTENT(OUT)::newKtot

  INTEGER::i,j,tt
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro

  ALLOCATE(filtro(N_nodes))

  filtro=.FALSE.

  DO i=1,newNnodes
     filtro(list_nodes(i)+1)=.TRUE.
  END DO

  newKtot=0
  DO i=1,newNnodes
     tt=list_nodes(i)+1
     DO j=T_start(tt)+1,T_start(tt+1)
        IF (filtro(T_ind(j))==.TRUE.) THEN
           newKtot=newKtot+1
        END IF
     END DO
  END DO

  DEALLOCATE(filtro)

END SUBROUTINE EXTRACT_NET_NEW_K_TOTAL

SUBROUTINE EXTRACT_NET(newKmax,TT_tau,TT_ind,TT_start,newKtot,list_nodes,T_ind,T_tau,T_start,N_nodes,Ktot,newNnodes)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::N_nodes,Ktot,newKtot,newNnodes
  INTEGER,DIMENSION(Ktot),INTENT(IN)::T_ind
  DOUBLE PRECISION,DIMENSION(Ktot),INTENT(IN)::T_tau
  INTEGER,DIMENSION(N_nodes+1),INTENT(IN)::T_start
  INTEGER,DIMENSION(newNnodes),INTENT(IN)::list_nodes

  INTEGER,DIMENSION(newNnodes+1),INTENT(OUT)::TT_start
  INTEGER,DIMENSION(newKtot),INTENT(OUT)::TT_ind
  DOUBLE PRECISION,DIMENSION(newKtot),INTENT(OUT)::TT_tau
  INTEGER,INTENT(OUT)::newKmax

  INTEGER::i,j,gg,hh,tt
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro
  INTEGER,DIMENSION(:),ALLOCATABLE::inv_trad

  ALLOCATE(filtro(N_nodes))
  ALLOCATE(inv_trad(N_nodes))

  filtro=.FALSE.
  inv_trad=0

  DO i=1,newNnodes
     filtro(list_nodes(i)+1)=.TRUE.
     inv_trad(list_nodes(i)+1)=i
  END DO

  gg=0
  hh=0
  DO tt=1,newNnodes
     i=list_nodes(tt)+1
     hh=0
     TT_start(inv_trad(i))=gg
     DO j=T_start(i)+1,T_start(i+1)
        IF (filtro(T_ind(j))==.TRUE.) THEN
           gg=gg+1
           hh=hh+1
           TT_tau(gg)=T_tau(j)
           TT_ind(gg)=inv_trad(T_ind(j))
        END IF
     END DO
     IF (newKmax<hh) newKmax=hh
  END DO
  TT_start(newNnodes+1)=gg

  DEALLOCATE(filtro,inv_trad)

END SUBROUTINE EXTRACT_NET


SUBROUTINE WEIGHT_CORE(newKmax,TT_tau,TT_ind,TT_start,trad,newKtot,newNnodes,threshold,T_ind,T_tau,T_start,N_nodes,Ktot)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::N_nodes,Ktot,newKtot,newNnodes
  DOUBLE PRECISION::threshold
  INTEGER,DIMENSION(Ktot),INTENT(IN)::T_ind
  DOUBLE PRECISION,DIMENSION(Ktot),INTENT(IN)::T_tau
  INTEGER,DIMENSION(N_nodes+1),INTENT(IN)::T_start

  INTEGER,DIMENSION(newNnodes+1),INTENT(OUT)::TT_start
  INTEGER,DIMENSION(newKtot),INTENT(OUT)::TT_ind
  INTEGER,DIMENSION(newNnodes),INTENT(OUT)::trad
  DOUBLE PRECISION,DIMENSION(newKtot),INTENT(OUT)::TT_tau
  INTEGER,INTENT(OUT)::newKmax

  INTEGER::i,j,gg,hh
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro
  INTEGER,DIMENSION(:),ALLOCATABLE::inv_trad
  DOUBLE PRECISION::Pe

  ALLOCATE(filtro(N_nodes),inv_trad(N_nodes))
  filtro=.TRUE.
  trad=0
  inv_trad=0
  newKmax=0

  gg=0
  DO i=1,N_nodes
     Pe=0.0d0
     DO j=T_start(i)+1,T_start(i+1)
        Pe=Pe+T_tau(j)
     END DO
     IF (Pe<=threshold) THEN
        filtro(i)=.FALSE.
     ELSE
        gg=gg+1
        trad(gg)=i
        inv_trad(i)=gg
     END IF
  END DO

  gg=0
  hh=0
  DO i=1,N_nodes
     IF (filtro(i)==.TRUE.) THEN
        hh=0
        TT_start(inv_trad(i))=gg
        DO j=T_start(i)+1,T_start(i+1)
           IF (filtro(T_ind(j))==.TRUE.) THEN
              gg=gg+1
              hh=hh+1
              TT_tau(gg)=T_tau(j)
              TT_ind(gg)=inv_trad(T_ind(j))
           END IF
        END DO
        IF (newKmax<hh) newKmax=hh
     END IF
  END DO
  TT_start(newNnodes+1)=gg

  DEALLOCATE(filtro,inv_trad)

END SUBROUTINE WEIGHT_CORE

  


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

SUBROUTINE GRAD_2 (N_sets,comunidades,dim,T_ind,T_tau,T_start,N_nodes,Ktot)   


  IMPLICIT NONE
  

  INTEGER,INTENT(IN)::N_nodes,Ktot,dim
  INTEGER,DIMENSION(Ktot),INTENT(IN)::T_ind
  DOUBLE PRECISION,DIMENSION(Ktot),INTENT(IN)::T_tau
  INTEGER,DIMENSION(N_nodes+1),INTENT(IN)::T_start


  INTEGER,INTENT(OUT)::N_sets
  INTEGER,DIMENSION(N_nodes),INTENT(OUT)::comunidades

  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::Pe


  INTEGER:: i,j,jj,g,h

  !! Para minimos:
  DOUBLE PRECISION::peso
  INTEGER::dim2,inicio,fin,candidato
  INTEGER,DIMENSION(:),ALLOCATABLE::tope_lista
  INTEGER,DIMENSION(:,:),ALLOCATABLE::lista
  LOGICAL,DIMENSION(:),ALLOCATABLE::label_glob_min,label_loc_min,filtro
  logical::inter

  !! Para las basins:
  INTEGER::hacia,desde
  logical,dimension(:),allocatable::label
  integer,dimension(:),allocatable::vect_aux1

  ALLOCATE(Pe(N_nodes))
  Pe=0.0d0

  DO i=1,N_nodes
     DO j=T_start(i)+1,T_start(i+1)
        Pe(i)=Pe(i)+T_tau(j)
     END DO
  END DO

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

END SUBROUTINE GRAD_2


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
  
  integer::AA,BB
  integer::i,j,times,jj,g

  double precision,dimension(:),allocatable::Pf,Pf2
  double precision,dimension(:,:),allocatable::appear
  
  logical,dimension(:),allocatable::filtro,filtro2
  double precision::cut
  double precision::delta_cut
  double precision::Z,Za,Zab

  AA=A+1
  BB=B+1

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

  Pf2(AA)=1.0d0
  Pf2(BB)=0.0d0


  DO times=1,num_iter

     print*,times

     Pf2(AA)=1.0d0
     Pf2(BB)=0.0d0
     Pf=Pf2
     Pf2=0.0d0
     DO i=1,N_nodes
        DO j=T_start(i)+1,T_start(i+1)
           jj=T_ind(j)
           Pf2(i)=Pf2(i)+(T_tau(j)*Pf(jj))/Pe(i)
        END DO
     END DO

  END DO

  Pf2(AA)=1.0d0
  Pf2(BB)=0.0d0
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


SUBROUTINE cfep_pfold2 (opt_bins,plot,node_index,A,B,T_ind,T_tau,T_start,length,num_iter,N_nodes,Ktot)

  implicit none

  INTEGER,INTENT(IN)::opt_bins,A,B,N_nodes,Ktot,length,num_iter
  INTEGER,DIMENSION(Ktot),INTENT(IN)::T_ind
  DOUBLE PRECISION,DIMENSION(Ktot),INTENT(IN)::T_tau
  INTEGER,DIMENSION(N_nodes+1),INTENT(IN)::T_start

  DOUBLE PRECISION,DIMENSION(length,3),INTENT(OUT)::plot
  INTEGER,DIMENSION(N_nodes),INTENT(OUT)::node_index

  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::Pe
  
  integer::AA,BB
  integer::i,j,times,jj,g,gg,ii,kk

  double precision,dimension(:),allocatable::Pf,Pf2

  INTEGER,dimension(:),allocatable::orderpf
  
  logical,dimension(:),allocatable::filtro,filtro2
  double precision::cut
  double precision::delta_cut
  double precision::Z,Za,Zab

  AA=A+1
  BB=B+1

  plot=0.0d0
  node_index=0

  ALLOCATE(Pe(N_nodes))
  Pe=0.0d0
  DO i=1,N_nodes
     DO j=T_start(i)+1,T_start(i+1)
        Pe(i)=Pe(i)+T_tau(j)
     END DO
  END DO

  ALLOCATE(Pf(N_nodes),Pf2(N_nodes))
  Pf=0.0d0
  Pf2=0.0d0

  Pf2(AA)=1.0d0
  Pf2(BB)=0.0d0


  DO times=1,num_iter

     Pf2(AA)=1.0d0
     Pf2(BB)=0.0d0
     Pf=Pf2
     Pf2=0.0d0
     DO i=1,N_nodes
        DO j=T_start(i)+1,T_start(i+1)
           jj=T_ind(j)
           Pf2(i)=Pf2(i)+(T_tau(j)*Pf(jj))/Pe(i)
        END DO
     END DO

  END DO

  Pf2(AA)=1.0d0
  Pf2(BB)=0.0d0
  Pf=Pf2
  
  DEALLOCATE(Pf2)

  IF (opt_bins==1) THEN

     ALLOCATE(filtro(N_nodes),filtro2(N_nodes))
     
     filtro2=.false.
     
     Z=sum(Pe(:),dim=1)
     
     delta_cut=1.0d0/(length*1.0d0)
     
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
                 node_index(i)=g-1
              END IF
           END IF
        END DO
        
     END DO
     
  ELSE

     Z=sum(Pe(:),dim=1)

     ALLOCATE(filtro(N_nodes),orderpf(N_nodes))

     filtro=.TRUE.

     DO i=1,N_nodes
        j=MAXLOC(Pf,DIM=1,MASK=filtro)
        orderpf(i)=j
        filtro(j)=.FALSE.
     END DO

     filtro=.false.
     Za=0.0d0

     DO i=1,N_nodes

        g=orderpf(i)
        filtro(g)=.true.
        Za=Za+Pe(g)
        Zab=0.0d0

        DO j=1,i
           gg=orderpf(j)
           DO jj=T_start(gg)+1,T_start(gg+1)
              kk=T_ind(jj)
              IF (filtro(kk).eqv..false.) THEN
                 Zab=Zab+T_tau(jj)
              END IF
           END DO
        END DO

        plot(i,1)=(Za/Z)
        plot(i,2)=-log(Zab/Z)
        plot(i,3)=Pf(g)
        node_index(g)=i

     END DO

  END IF


END SUBROUTINE cfep_pfold2

SUBROUTINE cfep_pfold3 (opt_bins,plot,node_index,A,B,T_ind,T_tau,T_start,length,num_iter,N_nodes,Ktot)

  implicit none

  TYPE int_pointer
     INTEGER,DIMENSION(:),POINTER::ip
  END TYPE int_pointer
  

  INTEGER,INTENT(IN)::opt_bins,A,B,N_nodes,Ktot,length,num_iter
  INTEGER,DIMENSION(Ktot),INTENT(IN)::T_ind
  DOUBLE PRECISION,DIMENSION(Ktot),INTENT(IN)::T_tau
  INTEGER,DIMENSION(N_nodes+1),INTENT(IN)::T_start

  DOUBLE PRECISION,DIMENSION(length,3),INTENT(OUT)::plot
  INTEGER,DIMENSION(N_nodes),INTENT(OUT)::node_index

  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::Pe
  
  integer::AA,BB
  integer::i,j,times,jj,g,gg,ii,kk,hh,tt

  double precision,dimension(:),allocatable::Pf,Pf2

  INTEGER,dimension(:),allocatable::orderpf
  
  logical,dimension(:),allocatable::filtro,filtro2
  double precision::cut
  double precision::delta_cut
  double precision::Z,Za,Zab,aux,aux_de_pf,aux_a_pf,aux_trans

  !Para ordenar
  integer::dim_buckets,num_occ_buckets
  logical::interr
  TYPE(int_pointer),DIMENSION(:),POINTER::buckets
  INTEGER,DIMENSION(:),ALLOCATABLE::occup_buckets,orden
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::valores
  double precision::lim_top,val_aux

  AA=A+1
  BB=B+1

  plot=0.0d0
  node_index=0
  print*,'entra'
  ALLOCATE(Pe(N_nodes))
  Pe=0.0d0
  DO i=1,N_nodes
     DO j=T_start(i)+1,T_start(i+1)
        Pe(i)=Pe(i)+T_tau(j)
     END DO
  END DO

  ALLOCATE(Pf(N_nodes),Pf2(N_nodes))
  Pf=0.0d0
  Pf2=0.0d0

  Pf2(AA)=1.0d0
  Pf2(BB)=0.0d0


  DO times=1,num_iter
     print*,times
     Pf2(AA)=1.0d0
     Pf2(BB)=0.0d0
     Pf=Pf2
     Pf2=0.0d0
     DO i=1,N_nodes
        DO j=T_start(i)+1,T_start(i+1)
           jj=T_ind(j)
           Pf2(i)=Pf2(i)+(T_tau(j)*Pf(jj))/Pe(i)
        END DO
     END DO

  END DO

  Pf2(AA)=1.0d0
  Pf2(BB)=0.0d0
  Pf=Pf2
  
  DEALLOCATE(Pf2)
  print*,'listo'
  IF (opt_bins==1) THEN

     ALLOCATE(filtro(N_nodes),filtro2(N_nodes))
     
     filtro2=.false.
     
     Z=sum(Pe(:),dim=1)
     
     delta_cut=1.0d0/(length*1.0d0)
     
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
                 node_index(i)=g-1
              END IF
           END IF
        END DO
        
     END DO
     
  ELSE

     Z=sum(Pe(:),dim=1)
     print*,'entra'

     ALLOCATE(orderpf(N_nodes))
     CALL sort_by_buckets (orderpf,0.0d0,1.0d0,100,2500,Pf,N_nodes)

!!$     print*,'ahi va viejo'
!!$
!!$     DO i=1,N_nodes
!!$
!!$        g=orderpf(i)
!!$        cut=Pf(g)
!!$
!!$        Za=0.0d0
!!$        Zab=0.0d0
!!$
!!$        DO j=1,N_nodes
!!$           IF (Pf(j)>cut) THEN
!!$              Za=Za+Pe(j)
!!$           END IF
!!$
!!$           DO jj=T_start(j)+1,T_start(j+1)
!!$              kk=T_ind(jj)
!!$              IF (((Pf(j)>cut).and.(Pf(kk)<=cut)).or.((Pf(j)<=cut).and.(Pf(kk)>cut))) THEN
!!$                 Zab=Zab+T_tau(jj)
!!$              END IF
!!$           END DO
!!$        END DO
!!$
!!$        plot(i,1)=(Za/Z)
!!$        plot(i,2)=-300.0d0*0.0020d0*log(Zab/Z)
!!$        plot(i,3)=Pf(g)
!!$        node_index(i)=g-1
!!$
!!$     END DO

     print*,'ahi va nuevo'

     aux=0.0d0
     plot(N_nodes,1)=aux
     DO i=N_nodes-1,1,-1
        g=orderpf(i+1)
        aux=aux+Pe(g)
        IF (Pf(orderpf(i))<Pf(g)) THEN
           plot(i,1)=aux
        ELSE
           plot(i,1)=plot(i+1,1)
        END IF
     END DO

     DO i=1,N_nodes
        print*,i
        g=orderpf(i)
        aux_de_pf=Pf(g)
        DO j=T_start(g)+1,T_start(g+1)
           kk=T_ind(j)
           aux_a_pf=Pf(kk)
           aux_trans=T_tau(j)
           IF (aux_a_pf>aux_de_pf) THEN
              DO jj=i,1,-1
                 gg=orderpf(jj)
                 IF (Pf(gg)<aux_de_pf) THEN
                    EXIT
                 ELSE
                    plot(jj,2)=plot(jj,2)+aux_trans
                 END IF
              END DO
              DO jj=i+1,N_nodes
                 gg=orderpf(jj)
                 IF (Pf(gg)>=aux_a_pf) THEN
                    EXIT
                 ELSE
                    plot(jj,2)=plot(jj,2)+aux_trans
                 END IF
              END DO
           ELSE IF (aux_a_pf<aux_de_pf) THEN
              !DO jj=i+1,N_nodes
              !   gg=orderpf(jj)
              !   IF (Pf(gg)>aux_de_pf) THEN
              !      EXIT
              !   ELSE
              !      plot(jj,2)=plot(jj,2)+aux_trans
              !   END IF
              !END DO
              DO jj=i,1,-1
                 gg=orderpf(jj)
                 IF (Pf(gg)<aux_de_pf) THEN
                    IF (Pf(gg)<aux_a_pf) THEN
                       EXIT
                    ELSE
                       plot(jj,2)=plot(jj,2)+aux_trans
                    END IF
                 END IF
              END DO
           END IF
        END DO
        plot(i,3)=aux_de_pf
        node_index(i)=g-1
     END DO

     DO i=1,N_nodes
        plot(i,1)=plot(i,1)/Z
        plot(i,2)=-300.0d0*0.0020d0*log(plot(i,2)/Z)
     END DO

  END IF


END SUBROUTINE cfep_pfold3


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
     DO j=T_start(i)+1,T_start(i+1)
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


SUBROUTINE DENDO_TIME (num_steps,N_basins,pertenece_a,T_ind,T_tau,T_start,N_nodes,Ktot)   !!!! Para revisar!!!! T_tau debe ser double precision


  IMPLICIT NONE
  

  INTEGER,INTENT(IN)::num_steps,N_nodes,Ktot,N_basins
  INTEGER,DIMENSION(Ktot),INTENT(IN)::T_ind
  DOUBLE PRECISION,DIMENSION(Ktot),INTENT(IN)::T_tau
  INTEGER,DIMENSION(N_nodes+1),INTENT(IN)::T_start
  INTEGER,DIMENSION(N_nodes),INTENT(IN)::pertenece_a

  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::Pe,poblacion,columna,columna2
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::posicion
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro,filtro2
  INTEGER,DIMENSION(:),ALLOCATABLE::representante,salto,vertical,etiqueta

  INTEGER::N,i,j,g,h,ii,jj,gg,hh,veces,dim2,dim
  INTEGER::factor,candidato
  DOUBLE PRECISION::aux
  LOGICAL::inter,bandera,bandera2

  dim=1

  ALLOCATE(representante(N_basins),poblacion(N_basins),Pe(N_nodes))
  representante=0
  poblacion=0.0d0
  Pe=0.0d0

  DO ii=1,N_nodes
     DO jj=T_start(ii)+1,T_start(ii+1)
        Pe(ii)=Pe(ii)+T_tau(jj)
     END DO
  END DO

  DO ii=1,N_nodes
     gg=pertenece_a(ii)+1
     IF (poblacion(gg)<Pe(ii)) THEN
        representante(gg)=ii
        poblacion(gg)=Pe(ii)
     END IF
  END DO

  N=N_nodes
  dim2=dim

  ALLOCATE(columna(N),filtro(N),columna2(N),filtro2(N))
  ALLOCATE(salto(N_basins),vertical(N_basins))
  salto=0
  vertical=0

  DO gg=1,N_basins

     filtro=.false.
     columna=0.0d0
     filtro2=.false.
     columna2=0.0d0

     g=representante(gg)

     DO i=T_start(g)+1,T_start(g+1)
        filtro(T_ind(i))=.true.
        columna(T_ind(i))=T_tau(i)
     END DO

     DO veces=1,num_steps

        filtro2=.false.
        columna2=0.0d0
   
        DO i=1,N
           IF (filtro(i)==.true.) THEN
              
              DO j=T_start(i)+1,T_start(i+1)
                 columna2(T_ind(j))=columna2(T_ind(j))+T_tau(j)*columna(i)
                 filtro2(T_ind(j))=.true.
              END DO
              
           END IF
        END DO

        filtro=.false.
        columna=0.0d0
        filtro=filtro2
        columna=columna2

        aux=0.0d0
        aux=sum(columna(:),DIM=1,MASK=filtro)
        
        !!print*,veces+1,COUNT(filtro(:),DIM=1),sum(columna(:),DIM=1,MASK=filtro)
        
        columna=columna/aux

        !! Detecto basin
        
        h=0
        inter=.false.
        
        DO WHILE (h<=dim2)
           
           candidato=maxloc(columna(:),DIM=1,MASK=filtro2)
           
           IF (Pe(candidato)>Pe(g)) THEN
              inter=.true.
              IF (pertenece_a(candidato)==pertenece_a(g)) THEN
                 print*,'tenemos un problema con representante',g,candidato
                 stop
              END IF
           ELSE
              filtro2(candidato)=.false.
              h=h+1
           END IF
           
           IF (inter==.true.) THEN
              exit
           END IF
           
        END DO
        
        IF (inter==.true.) THEN
           exit
        END IF
        
        
     END DO

     salto(gg)=veces+1
     vertical(gg)=pertenece_a(candidato)+1
     IF (salto(gg)==(num_steps+2)) salto(gg)=num_steps

  END DO
!!CLOSE(333)
  

  DEALLOCATE(filtro)
  ALLOCATE(etiqueta(N_basins),posicion(N_basins,3))
  ALLOCATE(filtro(N_basins))
  
  etiqueta=0
  filtro=.false.
  posicion=0.0d0
  
  DO i=1,N_basins
     etiqueta(i)=i
  END DO

  DO i=1,N_basins
     IF (salto(i)==num_steps) THEN
        vertical(i)=i
        posicion(i,2)=salto(i)
        posicion(i,3)=1.0d0
        filtro(i)=.true.
     END IF
  END DO

  DO i=1,N_basins
     IF (filtro(i)==.false.) THEN
        inter=.true.
        DO WHILE (inter==.true.)
           g=vertical(i)
           IF (salto(i)>salto(g)) THEN
              !            print*,'ya teniamos un problema'
              vertical(i)=vertical(g)
           ELSE
              inter=.false.
           END IF
        END DO
     END IF
  END DO

  factor=1

  DO i=num_steps,1,-1
     
     bandera=.true.
     DO WHILE (bandera==.true.)
        bandera=.false.
        DO j=1,N_basins
           IF ((salto(j)==i).and.(filtro(j)==.false.)) THEN
              bandera2=.false.
              IF (filtro(vertical(j))==.false.) THEN
                 print*,'Teniamos un problema',salto(j),salto(vertical(j))
                 IF (salto(j)==salto(vertical(j))) THEN
                    bandera=.true.
                    bandera2=.true.
                 END IF
              END IF
              
              IF (bandera2==.false.) THEN
                 factor=factor*(-1)
                 
                 !!WRITE(101,*)'#',factor
                 
                 posicion(j,2)=i
                 posicion(j,3)=posicion(vertical(j),3)+factor
                 !!WRITE(101,*)'#',posicion(j,2)
                 DO jj=1,N_basins
                    IF (factor>0) THEN
                       IF ((filtro(jj)==.true.).and.(posicion(jj,3)>=posicion(j,3))) THEN
                          posicion(jj,3)=posicion(jj,3)+factor
                       END IF
                    END IF
                    IF (factor<0) THEN
                       IF ((filtro(jj)==.true.).and.(posicion(jj,3)<=posicion(j,3))) THEN
                          posicion(jj,3)=posicion(jj,3)+factor
                       END IF
                    END IF
                 END DO
                 
                 filtro(j)=.true.
                 
              END IF
           END IF
        END DO
     END DO
     
  END DO
  
  WRITE(666,*)'####'
  DO i=1,N_basins
     
     WRITE(666,*)'#','altura',posicion(i,3),'minimo',etiqueta(i)
     WRITE(666,*)' '
     WRITE(666,*)posicion(i,1),posicion(i,3)
     WRITE(666,*)posicion(i,2),posicion(i,3)
     WRITE(666,*)' '
     
     IF ((vertical(i)/=0).and.(filtro(i)==.true.)) THEN
        WRITE(666,*)posicion(i,2),posicion(i,3)
        WRITE(666,*)posicion(i,2),posicion(vertical(i),3)
        WRITE(666,*)' '
     END IF
     
  END DO

  DEALLOCATE(representante,poblacion,Pe)
  DEALLOCATE(columna,filtro,columna2,filtro2)
  DEALLOCATE(salto,vertical)
  DEALLOCATE(etiqueta,posicion)

  
END SUBROUTINE DENDO_TIME


SUBROUTINE DENDO_BY_NODES (N_basins,pertenece_a,T_ind,T_tau,T_start,N_nodes,Ktot)

  IMPLICIT NONE

  TYPE int_pointer
     INTEGER,DIMENSION(:),POINTER::ip
  END TYPE int_pointer
  
  TYPE double_pointer
     DOUBLE PRECISION,DIMENSION(:),POINTER::dp
  END TYPE double_pointer

  INTEGER,INTENT(IN)::N_nodes,Ktot,N_basins
  INTEGER,DIMENSION(Ktot),INTENT(IN)::T_ind
  DOUBLE PRECISION,DIMENSION(Ktot),INTENT(IN)::T_tau
  INTEGER,DIMENSION(N_nodes+1),INTENT(IN)::T_start
  INTEGER,DIMENSION(N_nodes),INTENT(IN)::pertenece_a

  INTEGER::gg,hh,ii,jj,kk,ll
  INTEGER,DIMENSION(:),ALLOCATABLE::aux_vect_int

  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::Pe,Pe_basin
  INTEGER,DIMENSION(:),ALLOCATABLE::Popul_basin
  TYPE(int_pointer),DIMENSION(:),POINTER::basin

  

  ALLOCATE(Pe(N_nodes),Pe_basin(N_basins),Popul_basin(N_basins))

  Pe=0.0d0
  Pe_basin=0.0d0
  Popul_basin=0

  DO ii=1,N_nodes
     DO jj=T_start(ii)+1,T_start(ii+1)
        Pe(ii)=Pe(ii)+T_tau(jj)
     END DO
     kk=pertenece_a(ii)
     Pe_basin(kk)=Pe_basin(kk)+Pe(ii)
     Popul_basin(kk)=Popul_basin(kk)+1
  END DO

  ALLOCATE(basin(N_basins))
  DO ii=1,N_basins
     ALLOCATE(basin(ii)%ip(Popul_basin(ii)))
     basin(ii)%ip(:)=0
  END DO

  ALLOCATE(aux_vect_int(N_basins))
  aux_vect_int=0
  DO ii=1,N_nodes
     kk=pertenece_a(ii)
     gg=aux_vect_int(kk)+1
     basin(kk)%ip(gg)=ii
     aux_vect_int(kk)=gg
  END DO

  !! Tiros con el resto de basins

  
END SUBROUTINE DENDO_BY_NODES



SUBROUTINE DENDO_BOTTOM_UP (N_basins,pertenece_a,T_ind,T_tau,T_start,N_nodes,Ktot)   !!!! Para revisar!!!! T_tau debe ser double precision


  IMPLICIT NONE

  TYPE array_pointer
     INTEGER,DIMENSION(:),POINTER::p1
  END TYPE array_pointer
  TYPE(array_pointer),DIMENSION(:),POINTER::cluster
  
  TYPE doble_pointer
     DOUBLE PRECISION,DIMENSION(:),POINTER::dp2
  END TYPE doble_pointer


  INTEGER,INTENT(IN)::N_nodes,Ktot,N_basins
  INTEGER,DIMENSION(Ktot),INTENT(IN)::T_ind
  DOUBLE PRECISION,DIMENSION(Ktot),INTENT(IN)::T_tau
  INTEGER,DIMENSION(N_nodes+1),INTENT(IN)::T_start
  INTEGER,DIMENSION(N_nodes),INTENT(IN)::pertenece_a

  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::Pe,pob_repres
  INTEGER,DIMENSION(:),ALLOCATABLE::representante,poblacion

  INTEGER::N,i,j,g,h,ii,jj,gg,hh,veces,dim2,dim
  integer::ggg,k,kk,iji,dim_supra,contador
  logical,dimension(:),allocatable::filtro,suprafiltro

  integer,dimension(:),allocatable:: auxiliar

  !! Para las basins:

  TYPE(array_pointer),DIMENSION(:),POINTER::basin,listados
  integer,dimension(:),allocatable::unen,lista

  !! Para el dendograma:
  integer,dimension(:),allocatable:: compi,compi_orig,etiqueta
  double precision, dimension(:,:),allocatable:: ejey
  TYPE(doble_pointer),DIMENSION(:),POINTER::dendo,dendo2
  double precision, dimension(:),allocatable:: origen, ejex,comosi,comosi2
  integer::cant_cluster
  double precision::aux_1,aux_2
  double precision::referencia,lim_inf,yavalio
  double precision,dimension(:),allocatable::barrera
  integer::desde,hacia,nada3,nada4,limite
  logical::inter,inter2,inter3

  !! Para transformar el dendograma:
  double precision,dimension(:),allocatable::newejex,newhorizontal,horizontal,aux
  double precision,dimension(:,:),allocatable::newejey
  integer,dimension(:,:),allocatable::ramas
  logical,dimension(:),allocatable::filt_rama,filtro3,pisado
  logical,dimension(:,:),allocatable::filtro2
  double precision,dimension(:,:),allocatable::tree
  integer::candidato,cuantos
  integer::minimo,compruebo
  integer::num_steps

  integer,dimension(:,:),allocatable::suprabasins,nexo,ordeno_dendo
  integer,dimension(:),allocatable::ordeno,pertenece_cluster,poblacion_supra,supracompi,orden_total
  double precision,dimension(:),allocatable::suprabarrera
  logical,dimension(:),allocatable::ford,forda,fordb,selecciono
  integer::dim_bars,dim_total
  double precision::bandera

  print*,'Que pasaaaa'

  ALLOCATE(representante(N_basins),poblacion(N_basins),Pe(N_nodes),pob_repres(N_basins))
  representante=0
  poblacion=0
  Pe=0.0d0
  pob_repres=0.0d0

  DO ii=1,N_nodes
     DO jj=T_start(ii)+1,T_start(ii+1)
        Pe(ii)=Pe(ii)+T_tau(jj)
     END DO
  END DO

  DO ii=1,N_nodes
     gg=pertenece_a(ii)+1
     IF (pob_repres(gg)<Pe(ii)) THEN
        representante(gg)=ii
        pob_repres(gg)=Pe(ii)
     END IF
     poblacion(gg)=poblacion(gg)+1
  END DO

  DEALLOCATE(pob_repres)


  ALLOCATE(basin(N_basins))
  DO i=1,N_basins
     ALLOCATE(basin(i)%p1(poblacion(i)))
     basin(i)%p1(:)=0
  END DO

  poblacion=0

  DO i=1,N_nodes
     gg=pertenece_a(i)+1
     ii=poblacion(gg)+1
     poblacion(gg)=ii
     basin(gg)%p1(ii)=i
  END DO

  print*,'AQUI 1'
!!!!! Barro para el dendo
  limite=10000
  ALLOCATE(dendo(N_basins),dendo2(N_basins),comosi(limite),comosi2(limite),lista(limite))
  ALLOCATE(origen(N_basins),listados(N_basins))
  ALLOCATE(filtro(N_basins),suprafiltro(N_basins),etiqueta(N_basins))
  ALLOCATE(ejex(N_basins),ejey(N_basins,2),barrera(N_basins))
  ALLOCATE(compi(N_basins),compi_orig(N_basins))
  ALLOCATE(pertenece_cluster(N_basins))

  etiqueta=0
  !dendo=0.0d0
  comosi=0.0d0
  comosi2=0.0d0
  origen=0.0d0
  ejey=0.0d0
  ejex=-1.0d0
  compi=1
  compi_orig=0
  barrera=0.0d0

  DO i=1,N_basins
     etiqueta(i)=i
  END DO

  origen=0.0d0
  contador=0
  DO i=1,N_basins
     
     aux_1=0.0d0
     aux_2=0.0d0
     contador=0
     comosi=0.0d0
     comosi2=0.0d0
     lista=0
     
     DO ii=1,poblacion(i)
        
        g=basin(i)%p1(ii)
        
        
        aux_1=Pe(g)
        IF (origen(i)<aux_1) THEN
           origen(i)=aux_1
        END IF
        
        DO j=T_start(g)+1,T_start(g+1)
           
           IF ((pertenece_a(T_ind(j))+1)/=i) THEN
              
              aux_2=Pe(T_ind(j))
              gg=pertenece_a(T_ind(j))+1
              
              inter=.false.
              DO kk=1,contador
                 IF (lista(kk)==gg) THEN
                    yavalio=comosi(kk)
                    inter=.true.
                 END IF
              END DO
              IF (inter==.false.) yavalio=0.0d0
              
              IF ((yavalio<aux_1).and.(yavalio<aux_2)) THEN
                 inter=.false.
                 DO kk=1,contador
                    IF (lista(kk)==gg) THEN
                       inter=.true.
                       exit
                    END IF
                 END DO
                 IF (inter==.false.) THEN                  
                    contador=contador+1
                    IF (contador>limite) THEN
                       print*,'problema, subele el limite'
                       stop
                    END IF
                    lista(contador)=gg
                    kk=contador
                 END IF
                 
                 IF (aux_2>=aux_1) THEN
                    comosi(kk)=aux_1
                    comosi2(kk)=aux_2
                    !dendo(gg,i)=aux_1
                    !dendo2(gg,i)=aux_2
                 ELSE
                    comosi(kk)=aux_2
                    comosi2(kk)=aux_1
                    !dendo(gg,i)=aux_2
                    !dendo2(gg,i)=aux_1
                 END IF
              END IF
           END IF
        END DO
     END DO
     
     ALLOCATE(dendo(i)%dp2(contador),dendo2(i)%dp2(contador),listados(i)%p1(contador))
     DO j=1,contador
        dendo(i)%dp2(j)=comosi(j)
        dendo2(i)%dp2(j)=comosi2(j)
        listados(i)%p1(j)=lista(j)
     END DO
     
  END DO

  print*,'AQUI 2'

  DO i=1,N_basins
     DO j=1,size(dendo(i)%dp2(:))
        g=listados(i)%p1(j)
        inter=.false.
        DO ii=1,size(dendo(g)%dp2(:))
           h=listados(g)%p1(ii)
           IF (h==i) THEN
              inter=.true.
              exit
           END IF
        END DO
        IF (inter==.false.) THEN
           kk=size(dendo(g)%dp2(:))
           comosi=0.0d0
           lista=0
           kk=kk+1
           IF (kk>limite) THEN
              print*,'problema, subele el limite'
              stop
           END IF
           DO jj=1,kk-1
              comosi(jj)=dendo(g)%dp2(jj)
              lista(jj)=listados(g)%p1(jj)
           END DO
           comosi(kk)=0.0d0
           lista(kk)=i
           DEALLOCATE(dendo(g)%dp2, listados(g)%p1)
           ALLOCATE(dendo(g)%dp2(kk),listados(g)%p1(kk))
           dendo(g)%dp2(:)=comosi(:)
           listados(g)%p1(:)=lista(:)
        END IF
     END DO
  END DO

  DO i=1,N_basins
     DO j=1,size(dendo(i)%dp2(:))
        dendo(i)%dp2(j)=-log(dendo(i)%dp2(j))
     END DO
     origen(i)=-log(origen(i))
  END DO
  

  print*,'AQUI 3'

  DO i=1,N_basins
     DO j=1,size(dendo(i)%dp2(:))
        g=listados(i)%p1(j)
        IF (g>i) THEN
           
           DO ii=1,size(dendo(g)%dp2(:))
              h=listados(g)%p1(ii)
              IF (h==i) THEN
                 exit
              END IF
           END DO
           
           IF (dendo(i)%dp2(j)>dendo(g)%dp2(ii)) THEN
              dendo(i)%dp2(j)=dendo(g)%dp2(ii)
           END IF
           dendo(g)%dp2(ii)=dendo(i)%dp2(j)
           
        END IF
     END DO
  END DO

  DO i=1,N_basins
     aux_1=1000.0d0
     DO j=1,size(dendo(i)%dp2(:))
        g=listados(i)%p1(j)
        IF ((i/=g).and.(dendo(i)%dp2(j)<aux_1)) THEN
           aux_1=dendo(i)%dp2(j)
           compi_orig(i)=g
        END IF
     END DO
  END DO

  !! Pongo para coordenadas
  !! Tengo el dendo y los origenes
  
  compi=0
  ejex=0.0d0
  ejey=0.0d0
  
  DO i=1,N_basins
     ejey(i,1)=origen(i)
  END DO
  
  !! ordeno basins y barreras
  
  ALLOCATE(ordeno(N_basins))
  
  filtro=.true.
  DO i=1,N_basins
     iji=MINLOC(origen,DIM=1,MASK=filtro)
     filtro(iji)=.false.
     ordeno(i)=iji
  END DO
  
  DO i=1,N_basins
     print*,ordeno(i),origen(ordeno(i))
  END DO
  
  ALLOCATE(filtro2(N_basins,N_basins))
  filtro2=.false.

  gg=0
  DO i=1,N_basins
     DO j=1,size(dendo(i)%dp2(:))
        g=listados(i)%p1(j)
        IF (dendo(i)%dp2(j)<=1000.0) then
           filtro2(g,i)=.true.
           gg=gg+1
        END IF
     END DO
  END DO
  
  
  gg=gg/2
  print*,'dim dendo',gg,N_basins+gg
  
  ALLOCATE(ordeno_dendo(N_basins+gg,2))
  
  gg=0
  inter=.true.
  DO WHILE (inter==.true.)
     gg=gg+1
     print*,gg
     inter=.false.
     aux_1=1000.0d0
     DO i=1,N_basins
        DO ii=1,size(dendo(i)%dp2(:))
           g=listados(i)%p1(ii)
           IF (g>i) THEN
              IF (filtro2(i,g)==.true.) THEN
                 IF (dendo(i)%dp2(ii)<aux_1) THEN
                    ordeno_dendo(gg,1)=i
                    ordeno_dendo(gg,2)=g
                    aux_1=dendo(i)%dp2(ii)
                    inter=.true.
                 END IF
              END IF
           END IF
        END DO
     END DO
     IF (inter==.true.) THEN
        filtro2(ordeno_dendo(gg,1),ordeno_dendo(gg,2))=.false.
     END IF
  END DO
  

  gg=gg-1
  dim_bars=gg
  dim_total=gg+N_basins
  print*,'dim_total',dim_total
  deallocate(filtro2)
  ALLOCATE(orden_total(dim_total),selecciono(dim_total))
  selecciono=.false.
  print*,N_basins
  ii=1
  jj=1
  DO i=1,dim_total
     IF (ii<=N_basins) THEN
        DO kk=1,size(dendo(ordeno_dendo(jj,2))%dp2(:))
           IF (listados(ordeno_dendo(jj,2))%p1(kk)==ordeno_dendo(jj,1)) THEN
              inter=.true.
              exit
           END IF
        END DO
        IF (inter==.true.) THEN
           IF (origen(ordeno(ii))<=dendo(ordeno_dendo(jj,2))%dp2(kk)) THEN
              orden_total(i)=ii
              ii=ii+1
           ELSE
              selecciono(i)=.true.
              orden_total(i)=jj
              jj=jj+1
           END IF
        ELSE
           PRINT*,'OTRO PROBLEMA'
           STOP
        END IF
        
     ELSE
        orden_total(i)=jj
        jj=jj+1
        selecciono(i)=.true.
     END IF
  END DO
   
  DO i=1,dim_total
     g=orden_total(i)
     IF (selecciono(i)==.false.) THEN
        print*,'basin',origen(ordeno(g)),ordeno(g)
     END IF
     IF (selecciono(i)==.true.) THEN
        inter=.false.
        DO j=1,size(dendo(ordeno_dendo(g,2))%dp2(:))
           h=listados(ordeno_dendo(g,2))%p1(j)
           IF (h==ordeno_dendo(g,1)) THEN
              inter=.true.
              exit
           END IF
        END DO
        IF (inter==.true.) THEN
           print*,'barrera',dendo(ordeno_dendo(g,2))%dp2(j),ordeno_dendo(g,1),ordeno_dendo(g,2)
        ELSE
           print*,'problemica'
           stop
        END IF
     END IF
  END DO


  bandera=0.0d0
  ALLOCATE(cluster(N_basins))
  !DO i=1,N_basins
  !   ALLOCATE(cluster(i)%p1(N_basins+1))
  !   cluster(i)%p1(:)=0
  !END DO
  cant_cluster=0  
  pertenece_cluster=0
  
  DO i=1,dim_total
     g=orden_total(i)
     IF (selecciono(i)==.false.) THEN
        cant_cluster=cant_cluster+1
        ALLOCATE(cluster(cant_cluster)%p1(1))
        cluster(cant_cluster)%p1(1)=ordeno(g)
        pertenece_cluster(ordeno(g))=cant_cluster
     ELSE
        desde=pertenece_cluster(ordeno_dendo(g,1))
        hacia=pertenece_cluster(ordeno_dendo(g,2))
        IF (hacia>desde) THEN
           desde=pertenece_cluster(ordeno_dendo(g,2))
           hacia=pertenece_cluster(ordeno_dendo(g,1))
        END IF
        
        IF (desde/=hacia) THEN
           compi(cluster(desde)%p1(1))=cluster(hacia)%p1(1)
           inter=.false.
           DO kk=1,size(dendo(ordeno_dendo(g,2))%dp2(:))
              h=listados(ordeno_dendo(g,2))%p1(kk)
              IF (h==ordeno_dendo(g,1)) THEN
                 inter=.true.
                 exit
              END IF
           END DO
           IF (inter==.true.) THEN
              ejey(cluster(desde)%p1(1),2)=dendo(ordeno_dendo(g,2))%dp2(kk)
           ELSE
              print*,'probleeeema'
              exit
           END IF
           ii=size(cluster(hacia)%p1(:))
           jj=size(cluster(desde)%p1(:))
           allocate(auxiliar(ii+jj))
           auxiliar=0
           auxiliar(1:ii)=cluster(hacia)%p1(:)
           auxiliar(ii+1:ii+jj)=cluster(desde)%p1(:)
           deallocate(cluster(hacia)%p1)
           allocate(cluster(hacia)%p1(ii+jj))
           cluster(hacia)%p1(:)=auxiliar(:)
           deallocate(auxiliar)
           DO j=(desde+1),cant_cluster
              deallocate(cluster(j-1)%p1)
              iji=size(cluster(j)%p1(:))
              allocate(cluster(j-1)%p1(iji))
              cluster(j-1)%p1(:)=cluster(j)%p1(:)
           END DO
           deallocate(cluster(cant_cluster)%p1)
           cant_cluster=cant_cluster-1
           DO j=1,cant_cluster
              DO iji=1,size(cluster(j)%p1(:))
                 pertenece_cluster(cluster(j)%p1(iji))=j
              END DO
           END DO
        END IF
     END IF
  END DO
  

   
  gg=0
  DO i=1,cant_cluster
     print*,'@@@@'
     DO j=1,SIZE(cluster(i)%p1(:),DIM=1)
        g=cluster(i)%p1(j)
        gg=gg+1
        ejex(g)=dble(gg)
        print'(I,F14.2,2E)',g,ejex(g),ejey(g,1),ejey(g,2)
     END DO
  END DO

  
  gg=cluster(1)%p1(1)
  ejey(gg,2)=20.0d0         
      
  
  !DO i=1,N_basins
  !   DEALLOCATE(cluster(i)%p1)
  !END DO
  !DEALLOCATE(cluster)
  
  OPEN(700,FILE="dendo1_lineas3.oup",STATUS="REPLACE",ACTION="WRITE")
  OPEN(800,FILE="dendo1_puntos3.oup",STATUS="REPLACE",ACTION="WRITE")

  DO i=1,N_basins
     WRITE(700,*) '#'
     WRITE(700,*) ' '
     WRITE(700,*) ejex(i),ejey(i,1)
     WRITE(700,*) ejex(i),ejey(i,2)
     WRITE(700,*) '#'
     WRITE(700,*) ' '
     WRITE(700,*) ejex(i),ejey(i,2)
     IF (compi(i)==0) THEN
        WRITE(700,*) ejex(i),ejey(i,2)
     ELSE
        WRITE(700,*) ejex(compi(i)),ejey(i,2)
     END IF
  END DO

  DO i=1,N_basins
     WRITE(800,*) ejex(i),ejey(i,1)
  END DO
  
  
  WRITE(800,*) ' '
  DO i=1,N_basins
     IF (i/=gg) THEN
        WRITE(800,'(A,x,I3,x,A,x,I,x,A,x,F14.2)') '#  basin:',i,'con representante:',representante(i),'ejex:',ejex(i)
        WRITE(800,'(A,x,I,x,A,x,F8.4)') '#  con el compi:',compi_orig(i),'con ejex',ejex(compi_orig(i))
        WRITE(800,*) ' '
     ELSE
        WRITE(800,'(A,x,I3,x,A,x,I,x,A,x,F14.2)') '#  basin:',i,'con representante:',representante(i),'ejex:',ejex(i)
        WRITE(800,'(A,x,I,x,A,x,F14.2)') '#  con el compi:',compi_orig(i),'con ejex',0.0d0
        WRITE(800,*) ' '      
     END IF
  END DO
  




!!! Transformo el dendograma:
  
  num_steps=N_basins-1
  
  ALLOCATE(horizontal(N_basins))
  ALLOCATE(newejex(N_basins),newhorizontal(N_basins),newejey(N_basins,2))
  ALLOCATE(filtro2(N_basins,N_basins),filtro3(N_basins),pisado(N_basins),aux(N_basins))
  allocate(ramas(num_steps,2),filt_rama(N_basins))
  aux=0.0d0
  newejex=0.0d0
  newhorizontal=0.0d0
  horizontal=0.0d0
  newejey=0.0d0
  ramas=0.0d0
  filt_rama=.false.
  
  DO i=1,N_basins
     newejey(i,1)=ejey(i,1)
     newejey(i,2)=ejey(i,2)
     newejex(i)=ejex(i)
     IF (i/=gg) THEN
        horizontal(i)=ejex(compi(i))
     ELSE
        horizontal(i)=ejex(gg)
     END IF
  END DO

  OPEN(90,FILE="dendo2_lineas3.oup",STATUS="REPLACE",ACTION="WRITE")
  OPEN(91,FILE="dendo2_puntos3.oup",STATUS="REPLACE",ACTION="WRITE")
  
  
  DO i=1,N_basins
     WRITE(91,*) newejex(i),newejey(i,1)
  END DO
  
  WRITE(91,*) ' '
  DO i=1,N_basins
     WRITE(91,'(A,x,I,x,A,x,I,x,A,x,F14.2)') '#  basin:',i,'con representante:',representante(i),'ejex:',newejex(i)
     WRITE(91,*) ' '
  END DO
  
  
  minimo=minloc(ejey(:,1),DIM=1)
  
  k=0
  filtro2=.false.
  filtro3=.false.
  pisado=.false.
  
  DO i=1,N_basins
     aux(i)=ejex(i)
  END DO
  
  DO j=1,N_basins
     
     DO i=1,N_basins
        IF (i/=minimo) THEN
           IF ((horizontal(i)>(ejex(j)-0.10d0)).and.(horizontal(i)<(ejex(j)+0.10d0))) THEN
              filtro2(j,i)=.true.
              print*,'@@',ejex(j),ejex(i),j,i
           END IF
        END IF
     END DO
     
     inter3=.true.
     
     IF (count(filtro2(j,:),DIM=1)/=0) THEN
        
        candidato=minloc(ejey(:,2),DIM=1,MASK=filtro2(j,:))
        
!!!!
        DO ii=1,N_basins
           IF ((ii/=minimo).and.(ii/=candidato)) THEN
              IF ((horizontal(ii)>(ejex(candidato)-0.10d0)).and.(horizontal(ii)<(ejex(candidato)+0.10d0))) THEN
                 filtro3(j)=.true.
                 inter3=.false.
              END IF
           END IF
        END DO
        
        IF (inter3==.true.) THEN
           
           !      print*,ejex(j),ejex(candidato)
           
           k=k+1
           
           ramas(k,1)=j
           ramas(k,2)=candidato
           
           newejex(j)=ejex(j)
           newejex(candidato)=ejex(candidato)
           newejey(j,2)=ejey(candidato,2)
           newejey(candidato,2)=ejey(candidato,2)
           
           aux(j)=((newejex(candidato)-newejex(j))/2.0d0)+newejex(j)
           
           filtro2(j,candidato)=.false.
           
           WRITE(90,*) '#'
           WRITE(90,*) ' '
           WRITE(90,*) newejex(j),newejey(j,1)
           WRITE(90,*) newejex(j),newejey(j,2)
           WRITE(90,*) '#'
           WRITE(90,*) ' '
           WRITE(90,*) newejex(ramas(k,1)),newejey(ramas(k,1),2)
           WRITE(90,*) newejex(ramas(k,2)),newejey(ramas(k,2),2)
           WRITE(90,*) '#'
           WRITE(90,*) ' '
           WRITE(90,*) newejex(candidato),newejey(candidato,2)
           WRITE(90,*) newejex(candidato),newejey(candidato,1)
           
           newejey(j,2)=newejey(candidato,2)
           newejey(j,1)=newejey(j,2) !!
           
           pisado(j)=.true.
           pisado(candidato)=.true.
           
        END IF
     END IF
     
  END DO



  print*,'***'
  
  compruebo=0
  
  DO i=1,N_basins
     
     compruebo=compruebo+count(filtro2(:,i),DIM=1)
     
  END DO
  
  
  DO jj=1,N_basins
     IF (filtro3(jj)==.true.) THEN
        print*,'habia en',jj
     END IF
  END DO
  
  !DO jj=1,1
  DO WHILE (compruebo/=0)
     
     DO j=1,N_basins
        
        IF ((count(filtro2(j,:),DIM=1)/=0)) THEN
           
           candidato=minloc(ejey(:,2),DIM=1,MASK=filtro2(j,:))
           
           IF (count(filtro2(candidato,:),DIM=1)==0) THEN !!
              
              IF (filtro3(j)==.true.) THEN
                 
                 print*,'******',j,candidato,ejex(j),ejex(candidato)!,eje
                 
                 print*, '#'
                 print*, ' '
                 print*, aux(j),newejey(j,1)
                 print*, aux(j),ejey(candidato,2)    
                 
                 print*, '#'
                 print*, ' '
                 print*, aux(candidato),ejey(candidato,2)
                 print*, aux(candidato),newejey(candidato,1)    
                 
                 
                 print*, '#'
                 print*, ' '
                 print*, aux(j),ejey(candidato,2)
                 print*, aux(candidato),ejey(candidato,2)    
                 
                 
                 WRITE(90,*) '#'
                 WRITE(90,*) ' '
                 WRITE(90,*) aux(j),newejey(j,1)
                 WRITE(90,*) aux(j),ejey(candidato,2)    
                 
                 WRITE(90,*) '#'
                 WRITE(90,*) ' '
                 WRITE(90,*) aux(candidato),ejey(candidato,2)
                 WRITE(90,*) aux(candidato),newejey(candidato,1)    
                 
                 
                 WRITE(90,*) '#'
                 WRITE(90,*) ' '
                 WRITE(90,*) aux(j),ejey(candidato,2)
                 WRITE(90,*) aux(candidato),ejey(candidato,2)    
                 
                 newejey(j,2)=ejey(candidato,2)
                 newejey(j,1)=newejey(j,2) !!
                 
                 aux_1=aux(j)
                 aux(j)=((aux(candidato)-aux_1)/2.0d0)+aux_1
                 
                 filtro2(j,candidato)=.false.
                 pisado(j)=.true.
                 
              ELSE
                 
                 print*,'***',j,candidato,ejex(j),ejex(candidato)!,eje
                 
                 WRITE(90,*) '#'
                 WRITE(90,*) ' '
                 WRITE(90,*) aux(j),newejey(j,2)
                 WRITE(90,*) aux(j),ejey(candidato,2)    
                 
                 WRITE(90,*) '#'
                 WRITE(90,*) ' '
                 WRITE(90,*) aux(candidato),ejey(candidato,2)
                 WRITE(90,*) aux(candidato),newejey(candidato,1)    
                 
                 
                 WRITE(90,*) '#'
                 WRITE(90,*) ' '
                 WRITE(90,*) aux(j),ejey(candidato,2)
                 WRITE(90,*) aux(candidato),ejey(candidato,2)    
                 
                 newejey(j,2)=ejey(candidato,2)
                 newejey(j,1)=newejey(j,2) !!
                 
                 aux_1=aux(j)
                 aux(j)=((aux(candidato)-aux_1)/2.0d0)+aux_1
                 
                 filtro2(j,candidato)=.false.
                 pisado(j)=.true.
                 
              END IF !!
              
           END IF
           
        end if
        
     END DO
     
     compruebo=0
     DO i=1,N_basins
        
        compruebo=compruebo+count(filtro2(:,i),DIM=1)
        
     END DO
     print*,'compruebo',compruebo
     
     
  END DO
  !stop
  
  WRITE(90,*) '#'
  WRITE(90,*) ' '
  WRITE(90,*) aux(minimo),newejey(minimo,2)
  WRITE(90,*) aux(minimo),newejey(minimo,2)+5.0d0    


END SUBROUTINE DENDO_BOTTOM_UP

SUBROUTINE SORT_BY_BUCKETS (order,lim_inf,lim_sup,dim_boxes,max_popul,valores,N_vals)

  IMPLICIT NONE

  TYPE int_pointer
     INTEGER,DIMENSION(:),POINTER::ip
  END TYPE int_pointer

  INTEGER,INTENT(IN)::N_vals,max_popul,dim_boxes
  DOUBLE PRECISION,INTENT(IN)::lim_inf,lim_sup
  DOUBLE PRECISION,DIMENSION(N_vals),INTENT(IN)::valores
  INTEGER,DIMENSION(N_vals),INTENT(OUT)::order

  INTEGER::ii,jj,kk,hh,gg,tt,entra,iii
  INTEGER::num_buckets,num_occ_buckets,prox,topprox,num_occ_boxes
  INTEGER,DIMENSION(:),ALLOCATABLE::occ_buckets,reord,occ_boxes
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::lims_buckets,aux_lims,vect_aux
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::vals_aux
  DOUBLE PRECISION::linf,lsup
  TYPE(int_pointer),DIMENSION(:),POINTER::buckets
  LOGICAL::interr
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro

  num_occ_buckets=1
  num_buckets=1
  ALLOCATE(occ_buckets(1),buckets(1),lims_buckets(2,1))
  occ_buckets(1)=N_vals
  ALLOCATE(buckets(1)%ip(N_vals))
  DO ii=1,N_vals
     buckets(1)%ip(ii)=ii
  END DO
  lims_buckets(1,1)=lim_inf
  lims_buckets(2,1)=lim_sup
  
  ALLOCATE(occ_boxes(dim_boxes+1),reord(dim_boxes+1))
  ALLOCATE(aux_lims(2,10))
  topprox=10

  interr=.FALSE.
  IF (N_vals>max_popul) interr=.TRUE.
  DO WHILE (interr==.TRUE.)
     interr=.FALSE.
     prox=0
     gg=0
     entra=0
     DO ii=1,num_buckets
        IF (occ_buckets(ii)>max_popul) THEN
           num_occ_boxes=0
           occ_boxes=0
           entra=entra+1
           linf=lims_buckets(1,entra)
           lsup=lims_buckets(2,entra)
           IF (linf<lsup) THEN
              DO jj=1,occ_buckets(ii)
                 hh=buckets(ii)%ip(jj)
                 tt=int(((valores(hh)-linf)/(lsup-linf))*dim_boxes)+1
                 order(hh)=tt
                 occ_boxes(tt)=occ_boxes(tt)+1
              END DO
              kk=0
              DO jj=1,dim_boxes+1
                 IF (occ_boxes(jj)>0) THEN
                    num_occ_boxes=num_occ_boxes+1
                    kk=kk+1
                    reord(jj)=kk
                    IF (occ_boxes(jj)>max_popul) THEN
                       interr=.TRUE.
                       prox=prox+1
                       IF (prox>topprox) THEN
                          ALLOCATE(vect_aux(2,topprox))
                          vect_aux(:,:)=aux_lims(:,:)
                          DEALLOCATE(aux_lims)
                          ALLOCATE(aux_lims(2,topprox+10))
                          aux_lims(:,1:topprox)=vect_aux(:,:)
                          DEALLOCATE(vect_aux)
                          topprox=topprox+10
                       END IF
                       aux_lims(1,prox)=(jj-1)*((lsup-linf)/dim_boxes)+linf
                       aux_lims(2,prox)=jj*((lsup-linf)/dim_boxes)+linf
                    END IF
                 END IF
              END DO
              DO jj=1,occ_buckets(ii)
                 hh=buckets(ii)%ip(jj)
                 tt=order(hh)
                 order(hh)=reord(tt)+gg
              END DO
              gg=gg+num_occ_boxes
           ELSE
              gg=gg+1
              DO jj=1,occ_buckets(ii)
                 hh=buckets(ii)%ip(jj)
                 order(hh)=gg
              END DO
              prox=prox+1
              IF (prox>topprox) THEN
                 ALLOCATE(vect_aux(2,topprox))
                 vect_aux(:,:)=aux_lims(:,:)
                 DEALLOCATE(aux_lims)
                 ALLOCATE(aux_lims(2,topprox+10))
                 aux_lims(:,1:topprox)=vect_aux(:,:)
                 DEALLOCATE(vect_aux)
                 topprox=topprox+10
              END IF
              aux_lims(1,prox)=linf
              aux_lims(2,prox)=lsup
           END IF
        ELSE
           gg=gg+1
           DO jj=1,occ_buckets(ii)
              hh=buckets(ii)%ip(jj)
              order(hh)=gg
           END DO
        END IF
        DEALLOCATE(buckets(ii)%ip)
     END DO

     DEALLOCATE(lims_buckets)
     ALLOCATE(lims_buckets(2,prox))
     lims_buckets(:,:)=aux_lims(:,1:prox)

     DEALLOCATE(buckets,occ_buckets)
     num_buckets=gg
     ALLOCATE(buckets(num_buckets),occ_buckets(num_buckets))
     occ_buckets=0
     DO ii=1,N_vals
        tt=order(ii)
        occ_buckets(tt)=occ_buckets(tt)+1
     END DO
     DO ii=1,num_buckets
        ALLOCATE(buckets(ii)%ip(occ_buckets(ii)))
        occ_buckets(ii)=0
     END DO
     DO ii=1,N_vals
        tt=order(ii)
        gg=occ_buckets(tt)+1
        buckets(tt)%ip(gg)=ii
        occ_buckets(tt)=gg
     END DO

  END DO

  DEALLOCATE(aux_lims,lims_buckets)
  DEALLOCATE(occ_boxes,reord)

  gg=0
  DO ii=1,num_buckets
     tt=occ_buckets(ii)
     ALLOCATE(vals_aux(tt),filtro(tt))
     filtro=.TRUE.
     DO jj=1,tt
        hh=buckets(ii)%ip(jj)
        vals_aux(jj)=valores(hh)
     END DO
     DO jj=1,tt
        gg=gg+1
        hh=MINLOC(vals_aux,DIM=1,MASK=filtro)
        filtro(hh)=.FALSE.
        hh=buckets(ii)%ip(hh)
        order(gg)=hh
     END DO
     DEALLOCATE(vals_aux,filtro)
     DEALLOCATE(buckets(ii)%ip)
  END DO
  DEALLOCATE(buckets)
  DEALLOCATE(occ_buckets)

END SUBROUTINE SORT_BY_BUCKETS


END MODULE GLOB



!! f2py --f90flags='-fast' -c -m pyn_fort_net pyn_fort_net.f90
!!f2py --f90flags='-fast' -c -m pyn_fort_net pyn_fort_net.f90 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_def -lpthread
