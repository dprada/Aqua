MODULE GLOB
!!! traj: frame,num_part,dim

  INTEGER,DIMENSION(:),ALLOCATABLE::T_ind,T_start
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::T_tau
  INTEGER,DIMENSION(:,:),ALLOCATABLE::labels
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::labels_daux
  
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::distrib,distrib_x
  
  
  TYPE array_pointer
     INTEGER,DIMENSION(:),POINTER::p1
  END TYPE array_pointer
  
  TYPE darray_pointer
     DOUBLE PRECISION,DIMENSION(:),POINTER::d1
  END TYPE darray_pointer
  

CONTAINS

  SUBROUTINE free_memory_ts ()
    
    IF (ALLOCATED(T_ind))   DEALLOCATE(T_ind)
    IF (ALLOCATED(T_tau))   DEALLOCATE(T_tau)
    IF (ALLOCATED(T_start)) DEALLOCATE(T_start)
    IF (ALLOCATED(labels))  DEALLOCATE(labels)
    IF (ALLOCATED(labels_daux))  DEALLOCATE(labels_daux)
    
  END SUBROUTINE free_memory_ts
  
  SUBROUTINE traj2net(len_str,traj_full,ranges,num_frames,num_parts,dimensions,tray)
    
    IMPLICIT NONE
    INTEGER,INTENT(IN)::len_str,num_parts,num_frames,dimensions
    INTEGER,DIMENSION(num_frames,num_parts,dimensions),INTENT(IN)::traj_full
    INTEGER,DIMENSION(dimensions,2),INTENT(IN)::ranges
    INTEGER,DIMENSION(num_frames,num_parts),INTENT(OUT)::tray
    
    INTEGER::N_nodes,Ktot
    INTEGER::index
    INTEGER::ii,jj,gg,kk,ll,hh,hhh,aa,bb
    INTEGER::llevamos,cajon,contador
    INTEGER::x,x2
    
    INTEGER,DIMENSION(:,:),ALLOCATABLE::indice
    LOGICAL,DIMENSION(:,:),ALLOCATABLE::ocupado
    INTEGER,DIMENSION(:),ALLOCATABLE::deantes,aux
    CHARACTER*40::f1,f2,f3
    
    INTEGER,DIMENSION(:),ALLOCATABLE::SL,W,K_out
    TYPE(array_pointer),DIMENSION(:),POINTER::C,WK_out
    INTEGER::Kmax
    INTEGER,DIMENSION(:),POINTER::aux_puntero,aux_puntero2
    LOGICAL::switch
    
    CALL free_memory_ts ()
    
    
    IF ((dimensions>999999).or.(len_str>999999)) THEN
       print*, 'ERROR in fortran traj2net'
       stop
    END IF
    
    
    DO cajon=1,dimensions
       
       IF (cajon==1) THEN
          llevamos=1
          tray=1
       ELSE
          CALL SYSTEM ('mv trad_aux.aux trad_aux_old.aux')
       END IF
       
       ALLOCATE(ocupado(llevamos,ranges(cajon,1):ranges(cajon,2)),indice(llevamos,ranges(cajon,1):ranges(cajon,2)))
       ocupado=.false.
       
       DO ii=1,num_parts
          DO jj=1,num_frames
             bb=tray(jj,ii)
             aa=traj_full(jj,ii,cajon)
             ocupado(bb,aa)=.true.
          END DO
       END DO
       
       contador=0
       
       WRITE(f1,'(I6)') cajon
       WRITE(f2,'(I6)') len_str
       f3="(I,"//TRIM(ADJUSTL(f1))//"I"//TRIM(ADJUSTL(f2))//")"
       
       OPEN(21,FILE="trad_aux.aux",status="REPLACE",ACTION="WRITE")
       IF (cajon==1) THEN
          DO ii=1,llevamos
             DO jj=ranges(1,1),ranges(1,2)
                IF (ocupado(ii,jj).eqv..true.) THEN
                   contador=contador+1
                   indice(ii,jj)=contador
                   WRITE(21,f3) contador,jj
                END IF
             END DO
          END DO
       ELSE
          OPEN(61,FILE="trad_aux_old.aux",status="OLD",ACTION="READ")
          ALLOCATE(deantes(cajon-1))         
          DO ii=1,llevamos
             READ(61,*) aa,deantes(:)
             DO jj=ranges(cajon,1),ranges(cajon,2)
                IF (ocupado(ii,jj).eqv..true.) THEN
                   contador=contador+1
                   WRITE(21,f3) contador,deantes(:),jj
                   indice(ii,jj)=contador
                END IF
             END DO
          END DO
          DEALLOCATE(deantes)
          CLOSE(61)
       END IF
       
       CLOSE(21)
       
       DEALLOCATE(ocupado)
       
       DO ii=1,num_parts
          DO jj=1,num_frames
             bb=tray(jj,ii)
             aa=traj_full(jj,ii,cajon)
             tray(jj,ii)=indice(bb,aa)
          END DO
       END DO
       
       DEALLOCATE(indice)
       llevamos=contador
       
       IF (cajon/=1) CALL SYSTEM('rm trad_aux_old.aux')
       
       !print*,'>>',cajon,llevamos
    END DO
    
    N_nodes=llevamos  
    ALLOCATE(labels(N_nodes,dimensions),deantes(dimensions))
    OPEN(21,FILE="trad_aux.aux",status="OLD",ACTION="READ")
    DO ii=1,N_nodes
       READ(21,*) llevamos,deantes(:)
       labels(ii,:)=deantes(:)
    END DO
    CLOSE(21)
    DEALLOCATE(deantes)
    CALL SYSTEM('rm trad_aux.aux')
    ALLOCATE(C(N_nodes),K_out(N_nodes))
    ALLOCATE(WK_out(N_nodes),SL(N_nodes),W(N_nodes))
    ALLOCATE(aux_puntero(25000),aux_puntero2(25000))
    
    DO ii=1,N_nodes
       ALLOCATE(C(ii)%p1(1),WK_out(ii)%p1(1))
       C(ii)%p1(:)=0
       WK_out(ii)%p1(:)=0
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
    
    DO index=1,num_parts
       
       x=tray(1,index)
       contador=contador+1
       x2=tray(2,index)
       contador=contador+1
       
       DO ii=3,num_frames
          
          W(x)=W(x)+1
          
          switch=.true.
          DO gg=1,K_out(x)
             IF(C(x)%p1(gg)==x2) THEN
                WK_out(x)%p1(gg)=WK_out(x)%p1(gg)+1
                switch=.false.
                EXIT
             END IF
          END DO
          
          IF (switch.eqv..true.) THEN
             IF (K_out(x)>0) THEN
                gg=K_out(x)
                IF (gg>25000) THEN
                   print*,'tenemos un problema'
                   stop
                END IF
                aux_puntero(1:gg)=C(x)%p1(:)
                aux_puntero2(1:gg)=WK_out(x)%p1(:)
                DEALLOCATE(C(x)%p1,WK_out(x)%p1)
                ALLOCATE(C(x)%p1(gg+1),WK_out(x)%p1(gg+1))
                C(x)%p1(1:gg)=aux_puntero(:)
                WK_out(x)%p1(1:gg)=aux_puntero2(:)
                C(x)%p1(gg+1)=x2
                WK_out(x)%p1(gg+1)=1
                K_out(x)=K_out(x)+1
             ELSE
                WK_out(x)%p1(1)=1
                C(x)%p1(1)=x2
                K_out(x)=1
             END IF
          END IF
          
          x=x2  !!bin
          
          x2=tray(ii,index)
          contador=contador+1
          
       END DO
       
       W(x)=W(x)+1
       
       switch=.true.
       
       DO gg=1,K_out(x)
          IF(C(x)%p1(gg)==x2) THEN
             WK_out(x)%p1(gg)=WK_out(x)%p1(gg)+1
             switch=.false.
             EXIT
          END IF
       END DO
       
       IF (switch.eqv..true.) THEN
          IF (K_out(x)>0) THEN
             gg=K_out(x)
             IF (gg>25000) THEN
                print*,'tenemos un problema'
                stop
             END IF
             aux_puntero(1:gg)=C(x)%p1(:)
             aux_puntero2(1:gg)=WK_out(x)%p1(:)
             DEALLOCATE(C(x)%p1,WK_out(x)%p1)
             ALLOCATE(C(x)%p1(gg+1),WK_out(x)%p1(gg+1))
             C(x)%p1(1:gg)=aux_puntero(:)
             WK_out(x)%p1(1:gg)=aux_puntero2(:)
             C(x)%p1(gg+1)=x2
             WK_out(x)%p1(gg+1)=1
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
    
    Kmax=maxval(K_out(:))
    Ktot=sum(K_out(:))
    
    
    ALLOCATE(T_start(N_nodes+1),T_ind(Ktot),T_tau(Ktot))
    T_start=0
    T_ind=0
    T_tau=0
    
    ll=0
    DO ii=1,N_nodes
       
       hh=0
       switch=.false.
       DO gg=1,K_out(ii)
          IF (C(ii)%p1(gg)==ii) THEN
             switch=.true.
             exit
          END IF
       END DO
       
       T_start(ii)=ll
       IF (switch.eqv..true.) THEN
          ll=ll+1
          T_ind(ll)=ii
          hhh=WK_out(ii)%p1(gg)
          T_tau(ll)=hhh
          hh=hh+hhh
          DO jj=1,K_out(ii)
             IF (jj/=gg) THEN
                ll=ll+1
                T_ind(ll)=C(ii)%p1(jj)
                hhh=WK_out(ii)%p1(jj)
                T_tau(ll)=hhh
                hh=hh+hhh
             END IF
          END DO
       ELSE
          DO jj=1,K_out(ii)
             ll=ll+1
             T_ind(ll)=C(ii)%p1(jj)
             hhh=WK_out(ii)%p1(jj)
             hh=hh+hhh
             T_tau(ll)=hhh
          END DO
       END IF
       
    END DO
    
    T_start(N_nodes+1)=ll
    
    tray=tray-1
    DEALLOCATE(C,K_out)
    DEALLOCATE(WK_out,SL,W)
    DEALLOCATE(aux_puntero,aux_puntero2)
    
  end subroutine traj2net
  
  SUBROUTINE trajbinning2net(len_str,traj_full,ranges,bins,num_frames,num_parts,dimensions,tray)
    
    IMPLICIT NONE
    INTEGER,INTENT(IN)::len_str,num_parts,num_frames,dimensions
    DOUBLE PRECISION,DIMENSION(num_frames,num_parts,dimensions),INTENT(IN)::traj_full
    DOUBLE PRECISION,DIMENSION(dimensions,2),INTENT(IN)::ranges
    INTEGER,DIMENSION(dimensions),INTENT(IN)::bins
    INTEGER,DIMENSION(num_frames,num_parts),INTENT(OUT)::tray
    
    DOUBLE PRECISION,DIMENSION(dimensions)::delta_l

    INTEGER::N_nodes,Ktot
    INTEGER::index
    INTEGER::ii,jj,gg,kk,ll,hh,hhh,aa,bb
    INTEGER::llevamos,cajon,contador
    INTEGER::x,x2
    
    INTEGER,DIMENSION(:,:),ALLOCATABLE::indice
    LOGICAL,DIMENSION(:,:),ALLOCATABLE::ocupado
    INTEGER,DIMENSION(:),ALLOCATABLE::deantes,aux
    CHARACTER*40::f1,f2,f3
    
    INTEGER,DIMENSION(:),ALLOCATABLE::SL,W,K_out
    TYPE(array_pointer),DIMENSION(:),POINTER::C,WK_out
    INTEGER::Kmax
    INTEGER,DIMENSION(:),POINTER::aux_puntero,aux_puntero2
    LOGICAL::switch
    
    delta_l=0.0d0
    DO ii=1,dimensions
       delta_l(ii)=(ranges(ii,2)-ranges(ii,1))/(1.0d0*bins(ii))
    END DO

    CALL free_memory_ts ()
    
    
    IF ((dimensions>999999).or.(len_str>999999)) THEN
       print*, 'ERROR in fortran traj2net'
       stop
    END IF
    
    
    DO cajon=1,dimensions
       
       IF (cajon==1) THEN
          llevamos=1
          tray=1
       ELSE
          CALL SYSTEM ('mv trad_aux.aux trad_aux_old.aux')
       END IF
       
       ALLOCATE(ocupado(llevamos,bins(cajon)),indice(llevamos,bins(cajon)))
       ocupado=.false.
       
       DO ii=1,num_parts
          DO jj=1,num_frames
             bb=tray(jj,ii)
             aa=traj_full(jj,ii,cajon)
             ocupado(bb,aa)=.true.
          END DO
       END DO
       
       contador=0
       
       WRITE(f1,'(I6)') cajon
       WRITE(f2,'(I6)') len_str
       f3="(I,"//TRIM(ADJUSTL(f1))//"I"//TRIM(ADJUSTL(f2))//")"
       
       OPEN(21,FILE="trad_aux.aux",status="REPLACE",ACTION="WRITE")
       IF (cajon==1) THEN
          DO ii=1,llevamos
             DO jj=ranges(1,1),ranges(1,2)
                IF (ocupado(ii,jj).eqv..true.) THEN
                   contador=contador+1
                   indice(ii,jj)=contador
                   WRITE(21,f3) contador,jj
                END IF
             END DO
          END DO
       ELSE
          OPEN(61,FILE="trad_aux_old.aux",status="OLD",ACTION="READ")
          ALLOCATE(deantes(cajon-1))         
          DO ii=1,llevamos
             READ(61,*) aa,deantes(:)
             DO jj=ranges(cajon,1),ranges(cajon,2)
                IF (ocupado(ii,jj).eqv..true.) THEN
                   contador=contador+1
                   WRITE(21,f3) contador,deantes(:),jj
                   indice(ii,jj)=contador
                END IF
             END DO
          END DO
          DEALLOCATE(deantes)
          CLOSE(61)
       END IF
       
       CLOSE(21)
       
       DEALLOCATE(ocupado)
       
       DO ii=1,num_parts
          DO jj=1,num_frames
             bb=tray(jj,ii)
             aa=traj_full(jj,ii,cajon)
             tray(jj,ii)=indice(bb,aa)
          END DO
       END DO
       
       DEALLOCATE(indice)
       llevamos=contador
       
       IF (cajon/=1) CALL SYSTEM('rm trad_aux_old.aux')
       
       !print*,'>>',cajon,llevamos
    END DO
    
    N_nodes=llevamos  
    ALLOCATE(labels(N_nodes,dimensions),deantes(dimensions))
    OPEN(21,FILE="trad_aux.aux",status="OLD",ACTION="READ")
    DO ii=1,N_nodes
       READ(21,*) llevamos,deantes(:)
       labels(ii,:)=deantes(:)
    END DO
    CLOSE(21)
    DEALLOCATE(deantes)
    CALL SYSTEM('rm trad_aux.aux')
    ALLOCATE(C(N_nodes),K_out(N_nodes))
    ALLOCATE(WK_out(N_nodes),SL(N_nodes),W(N_nodes))
    ALLOCATE(aux_puntero(25000),aux_puntero2(25000))
    
    DO ii=1,N_nodes
       ALLOCATE(C(ii)%p1(1),WK_out(ii)%p1(1))
       C(ii)%p1(:)=0
       WK_out(ii)%p1(:)=0
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
    
    DO index=1,num_parts
       
       x=tray(1,index)
       contador=contador+1
       x2=tray(2,index)
       contador=contador+1
       
       DO ii=3,num_frames
          
          W(x)=W(x)+1
          
          switch=.true.
          DO gg=1,K_out(x)
             IF(C(x)%p1(gg)==x2) THEN
                WK_out(x)%p1(gg)=WK_out(x)%p1(gg)+1
                switch=.false.
                EXIT
             END IF
          END DO
          
          IF (switch.eqv..true.) THEN
             IF (K_out(x)>0) THEN
                gg=K_out(x)
                IF (gg>25000) THEN
                   print*,'tenemos un problema'
                   stop
                END IF
                aux_puntero(1:gg)=C(x)%p1(:)
                aux_puntero2(1:gg)=WK_out(x)%p1(:)
                DEALLOCATE(C(x)%p1,WK_out(x)%p1)
                ALLOCATE(C(x)%p1(gg+1),WK_out(x)%p1(gg+1))
                C(x)%p1(1:gg)=aux_puntero(:)
                WK_out(x)%p1(1:gg)=aux_puntero2(:)
                C(x)%p1(gg+1)=x2
                WK_out(x)%p1(gg+1)=1
                K_out(x)=K_out(x)+1
             ELSE
                WK_out(x)%p1(1)=1
                C(x)%p1(1)=x2
                K_out(x)=1
             END IF
          END IF
          
          x=x2  !!bin
          
          x2=tray(ii,index)
          contador=contador+1
          
       END DO
       
       W(x)=W(x)+1
       
       switch=.true.
       
       DO gg=1,K_out(x)
          IF(C(x)%p1(gg)==x2) THEN
             WK_out(x)%p1(gg)=WK_out(x)%p1(gg)+1
             switch=.false.
             EXIT
          END IF
       END DO
       
       IF (switch.eqv..true.) THEN
          IF (K_out(x)>0) THEN
             gg=K_out(x)
             IF (gg>25000) THEN
                print*,'tenemos un problema'
                stop
             END IF
             aux_puntero(1:gg)=C(x)%p1(:)
             aux_puntero2(1:gg)=WK_out(x)%p1(:)
             DEALLOCATE(C(x)%p1,WK_out(x)%p1)
             ALLOCATE(C(x)%p1(gg+1),WK_out(x)%p1(gg+1))
             C(x)%p1(1:gg)=aux_puntero(:)
             WK_out(x)%p1(1:gg)=aux_puntero2(:)
             C(x)%p1(gg+1)=x2
             WK_out(x)%p1(gg+1)=1
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
    
    Kmax=maxval(K_out(:))
    Ktot=sum(K_out(:))
    
    
    ALLOCATE(T_start(N_nodes+1),T_ind(Ktot),T_tau(Ktot))
    T_start=0
    T_ind=0
    T_tau=0
    
    ll=0
    DO ii=1,N_nodes
       
       hh=0
       switch=.false.
       DO gg=1,K_out(ii)
          IF (C(ii)%p1(gg)==ii) THEN
             switch=.true.
             exit
          END IF
       END DO
       
       T_start(ii)=ll
       IF (switch.eqv..true.) THEN
          ll=ll+1
          T_ind(ll)=ii
          hhh=WK_out(ii)%p1(gg)
          T_tau(ll)=hhh
          hh=hh+hhh
          DO jj=1,K_out(ii)
             IF (jj/=gg) THEN
                ll=ll+1
                T_ind(ll)=C(ii)%p1(jj)
                hhh=WK_out(ii)%p1(jj)
                T_tau(ll)=hhh
                hh=hh+hhh
             END IF
          END DO
       ELSE
          DO jj=1,K_out(ii)
             ll=ll+1
             T_ind(ll)=C(ii)%p1(jj)
             hhh=WK_out(ii)%p1(jj)
             hh=hh+hhh
             T_tau(ll)=hhh
          END DO
       END IF
       
    END DO
    
    T_start(N_nodes+1)=ll
    
    tray=tray-1
    DEALLOCATE(C,K_out)
    DEALLOCATE(WK_out,SL,W)
    DEALLOCATE(aux_puntero,aux_puntero2)
    
  end subroutine trajbinning2net
  

  SUBROUTINE trajnodes2trajclusters (aux_list,trajnodes,num_nodes,frames,particles,trajout)
    
    IMPLICIT NONE
    INTEGER,INTENT(IN)::frames,particles,num_nodes
    INTEGER,DIMENSION(frames,particles,1),INTENT(IN)::trajnodes
    INTEGER,DIMENSION(num_nodes),INTENT(IN)::aux_list
    INTEGER,DIMENSION(frames,particles,1),INTENT(OUT)::trajout
    
    INTEGER::ii,jj,kk,gg,ll
    
    DO ii=1,frames
       DO jj=1,particles
          ll=trajnodes(ii,jj,1)
          gg=aux_list(ll+1)
          trajout(ii,jj,1)=gg
       END DO
    END DO
    
  END SUBROUTINE trajnodes2trajclusters
  
  
  SUBROUTINE PCA (num_eigenvs,traj,frames,num_parts,dims,eigenvals,eigenvects)
    
    IMPLICIT NONE
    
    INTEGER,INTENT(IN)::frames,num_parts,dims,num_eigenvs
    DOUBLE PRECISION,DIMENSION(frames,num_parts,dims),INTENT(IN)::traj
    DOUBLE PRECISION,DIMENSION(num_eigenvs),INTENT(OUT)::eigenvals
    DOUBLE PRECISION,DIMENSION(num_eigenvs,num_parts,dims),INTENT(OUT)::eigenvects
    
    INTEGER::ii,jj,kk,gg,dim_mat
    INTEGER::xx,yy,jj1,kk1,jj2,kk2
    DOUBLE PRECISION::aux_coor,aver_xx
    DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::C,CC
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::vals,aver
    INTEGER,DIMENSION(:,:),ALLOCATABLE::aux_ind
    DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::vects
    
    !Para diagonalizar:
    INTEGER::num_val,info,eigenlimit,lwork
    INTEGER, DIMENSION(:), ALLOCATABLE::iwork,ifail
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::work
    
    
    eigenvals=0.0d0
    eigenvects=0.0d0
    
    dim_mat=num_parts*dims

    ALLOCATE(C(dim_mat,dim_mat),aver(dim_mat),aux_ind(num_parts,dims))
    
    C=0.0d0
    aver=0.0d0
    aux_ind=0
    
    
    gg=0
    DO kk1=1,dims
       DO jj1=1,num_parts
          gg=gg+1
          aux_ind(jj1,kk1)=gg
       END DO
    END DO
    
    DO ii=1,frames
       
       DO jj1=1,num_parts
          DO kk1=1,dims
             xx=aux_ind(jj1,kk1)
             aux_coor=traj(ii,jj1,kk1)
             
             aver(xx)=aver(xx)+aux_coor
             
             DO jj2=jj1,num_parts
                DO kk2=kk1,dims
                   yy=aux_ind(jj2,kk2)
                   
                   C(xx,yy)=C(xx,yy)+aux_coor*traj(ii,jj2,kk2)
                   
                END DO
             END DO
          END DO
       END DO
       
    END DO
    
    aver=aver/(frames*1.0d0)
    
    DO kk1=1,dims
       DO jj1=1,num_parts
          xx=aux_ind(jj1,kk1)
          aver_xx=aver(xx)
          DO jj2=jj1,num_parts
             DO kk2=kk1,dims
                yy=aux_ind(jj2,kk2)
                C(xx,yy)=C(xx,yy)/(frames*1.0d0)-aver_xx*aver(yy)
                C(yy,xx)=C(xx,yy)
             END DO
          END DO
       END DO
    END DO
    
    eigenlimit=dim_mat-num_eigenvs+1
    lwork=8*dim_mat
    ALLOCATE(work(lwork),iwork(5*dim_mat),ifail(dim_mat))
    ALLOCATE(vects(num_eigenvs,dim_mat),vals(num_eigenvs))
    vects=0.0d0
    vals=0.0d0
    work=0.0d0
    iwork=0
    ifail=0
    
    CALL dsyevx ('V','I','U',dim_mat,C,dim_mat,0,0,eigenlimit,dim_mat,0.0d0,num_val&
         &,vals,vects,dim_mat,work,8*dim_mat,iwork,ifail,info)
    
    IF (info/=0) THEN
       
       print*,'# Error diagonalising the covariance matrix.'
       print*,'the array "work" should have the dimension:', work(1)
       STOP
       
    END IF
    
    DO ii=1,num_eigenvs
       gg=num_eigenvs-ii+1
       eigenvals(ii)=vals(gg)
    END DO
    
    DO ii=1,num_eigenvs
       gg=num_eigenvs-ii+1
       DO jj1=1,num_parts
          DO kk1=1,dims
             xx=aux_ind(jj1,kk1)
             eigenvects(ii,jj1,kk1)=vects(gg,xx)
          END DO
       END DO
    END DO
    
    DEALLOCATE(aux_ind,aver)
    DEALLOCATE(C,work,iwork,ifail)
    DEALLOCATE(vals,vects)
    
  END SUBROUTINE PCA
  

  subroutine prada1 (ybins,bins,min,max,delta_x,rv_min,rv_max,traj,tw,num_parts,len_traj,traj_out)
    
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ybins,bins,tw,len_traj,num_parts,rv_min,rv_max
    DOUBLE PRECISION,INTENT(IN)::delta_x,min,max
    DOUBLE PRECISION,DIMENSION(len_traj,num_parts,1),INTENT(IN)::traj
    INTEGER,DIMENSION(len_traj-2*tw,num_parts,bins),INTENT(OUT)::traj_out
    
    INTEGER::Ltw,Ltw1
    INTEGER::nn,ii,tt,kk,gg
    DOUBLE PRECISION::delta_y
    LOGICAL::filt_min,filt_max
    
    traj_out=0
    Ltw=(2*tw+1)
    Ltw1=Ltw-1
    delta_y=(1.0d0*ybins)/(1.0d0*Ltw)
    
    filt_min=.FALSE.
    filt_max=.FALSE.
    
    IF (rv_min==1) THEN
       filt_min=.TRUE.
    END IF
    IF (rv_max==1) THEN
       filt_max=.TRUE.
    END IF

    IF ((filt_min.eqv..TRUE.).OR.(filt_max.eqv..TRUE.)) THEN
       DO nn=1,num_parts
          DO ii=1,len_traj
             tt=INT((traj(ii,nn,1)-min)/delta_x)+1
             IF (filt_min.eqv..TRUE.) THEN
                IF (tt<1) THEN
                   tt=1
                ELSE
                   tt=tt+1
                END IF
             END IF
             IF (tt>bins) tt=bins
             DO kk=ii-tw,ii+tw
                gg=kk-tw
                IF ((gg>0).and.(gg<=(len_traj-Ltw1))) THEN
                   traj_out(gg,nn,tt)=traj_out(gg,nn,tt)+1
                END IF
             END DO
          END DO
       END DO
    ELSE
       DO nn=1,num_parts
          DO ii=1,len_traj
             tt=INT((traj(ii,nn,1)-min)/delta_x)+1
             IF (tt>bins) THEN
                tt=bins
             END IF
             DO kk=ii-tw,ii+tw
                gg=kk-tw
                IF ((gg>0).and.(gg<=(len_traj-Ltw1))) THEN
                   traj_out(gg,nn,tt)=traj_out(gg,nn,tt)+1
                END IF
             END DO
          END DO
       END DO
    END IF
    
    DO nn=1,num_parts
       DO ii=1,len_traj-Ltw1
          DO kk=1,bins
             tt=CEILING(delta_y*traj_out(ii,nn,kk))
             traj_out(ii,nn,kk)=tt
          END DO
       END DO
    END DO
    
    
  END SUBROUTINE prada1
  

  subroutine prada11 (ybins,bins,min,max,delta_x,rv_min,rv_max,traj,tw,num_parts,len_traj,traj_out)
    
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ybins,bins,tw,len_traj,num_parts,rv_min,rv_max
    DOUBLE PRECISION,INTENT(IN)::delta_x,min,max
    DOUBLE PRECISION,DIMENSION(len_traj,num_parts,1),INTENT(IN)::traj
    INTEGER,DIMENSION(len_traj-2*tw,num_parts,bins),INTENT(OUT)::traj_out
    
    INTEGER::Ltw,Ltw1
    INTEGER::nn,ii,tt,kk,gg,qq
    DOUBLE PRECISION::delta_y
    LOGICAL::filt_min,filt_max
    INTEGER,DIMENSION(:),ALLOCATABLE::histo_aux

    traj_out=0
    Ltw=(2*tw+1)
    Ltw1=Ltw-1
    delta_y=(1.0d0*ybins)/(1.0d0*Ltw)
    
    filt_min=.FALSE.
    filt_max=.FALSE.
    
    IF (rv_min==1) THEN
       filt_min=.TRUE.
    END IF
    IF (rv_max==1) THEN
       filt_max=.TRUE.
    END IF

    ALLOCATE(histo_aux(bins))

    IF ((filt_min.eqv..TRUE.).OR.(filt_max.eqv..TRUE.)) THEN
       
       DO nn=1,num_parts
          histo_aux=0
          DO ii=1,Ltw
             tt=INT((traj(ii,nn,1)-min)/delta_x)+1
             IF (filt_min.eqv..TRUE.) THEN
                IF (tt<1) THEN
                   tt=1
                ELSE
                   tt=tt+1
                END IF
             END IF
             IF (tt>bins) tt=bins
             histo_aux(tt)=histo_aux(tt)+1
          END DO
          DO ii=1,len_traj-Ltw
             traj_out(tw+ii,nn,:)=histo_aux(:)
             qq=INT((traj(Ltw+ii,nn,1)-min)/delta_x)+1
             IF (filt_min.eqv..TRUE.) THEN
                IF (qq<1) THEN
                   qq=1
                ELSE
                   qq=qq+1
                END IF
             END IF
             IF (qq>bins) qq=bins
             tt=INT((traj(ii,nn,1)-min)/delta_x)+1
             IF (filt_min.eqv..TRUE.) THEN
                IF (tt<1) THEN
                   tt=1
                ELSE
                   tt=tt+1
                END IF
             END IF
             IF (tt>bins) tt=bins
             histo_aux(qq)=histo_aux(qq)+1
             histo_aux(tt)=histo_aux(tt)-1
          END DO
          traj_out(len_traj-Ltw1,nn,:)=histo_aux(:)
       END DO
       
       
    ELSE

       DO nn=1,num_parts
          histo_aux=0
          DO ii=1,Ltw
             tt=INT((traj(ii,nn,1)-min)/delta_x)+1
             IF (tt>bins) tt=bins
             histo_aux(tt)=histo_aux(tt)+1
          END DO
          DO ii=1,len_traj-Ltw
             traj_out(tw+ii,nn,:)=histo_aux(:)
             qq=INT((traj(Ltw+ii,nn,1)-min)/delta_x)+1
             IF (qq>bins) qq=bins
             tt=INT((traj(ii,nn,1)-min)/delta_x)+1
             IF (tt>bins) tt=bins
             histo_aux(qq)=histo_aux(qq)+1
             histo_aux(tt)=histo_aux(tt)-1
          END DO
          traj_out(len_traj-Ltw1,nn,:)=histo_aux(:)
       END DO

    END IF

    DO nn=1,num_parts
       DO ii=1,len_traj-Ltw1
          DO kk=1,bins
             tt=CEILING(delta_y*traj_out(ii,nn,kk))
             traj_out(ii,nn,kk)=tt
          END DO
       END DO
    END DO
    
    
  END SUBROUTINE prada11
  


  subroutine prada2 (ybins,sbins,bins,min,max,delta_x,traj,tw,num_parts,len_traj,traj_out)
    
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ybins,sbins,bins,tw,len_traj,num_parts
    DOUBLE PRECISION,INTENT(IN)::delta_x,min,max
    DOUBLE PRECISION,DIMENSION(len_traj,num_parts,1),INTENT(IN)::traj
    INTEGER,DIMENSION(len_traj-2*tw,num_parts,bins+1),INTENT(OUT)::traj_out
    
    INTEGER::Ltw,Ltw1
    INTEGER::nn,ii,tt,kk,gg
    DOUBLE PRECISION::delta_y,sigma
    
    traj_out=0
    Ltw=(2*tw+1)
    Ltw1=Ltw-1
    delta_y=(1.0d0*ybins)/(1.0d0*Ltw)
    
    
    DO nn=1,num_parts
       DO ii=tw+1,len_traj-tw-1
          DO kk=ii-tw,ii+tw
             tt=CEILING((traj(kk,nn,1)-min)/delta_x)
             IF (tt>bins) tt=bins
             traj_out(ii,nn,tt)=traj_out(ii,nn,tt)+1
          END DO
          sigma=0.0d0
          DO kk=ii-tw,ii+tw-1
             sigma=sigma+(traj(kk+1,nn,1)-traj(kk,nn,1))**2
          END DO
          traj_out(ii,nn,bins+1)=int((sigma/(1.0d0*Ltw))*sbins)
       END DO
    END DO
    
    DO nn=1,num_parts
       DO ii=1,len_traj-Ltw1
          DO kk=1,bins
             tt=CEILING(delta_y*traj_out(ii,nn,kk))
             traj_out(ii,nn,kk)=tt
          END DO
       END DO
    END DO
    
  END SUBROUTINE prada2
  
  
  
  
  
  
  subroutine ganna (opt_range,opt,ibins,imin,imax,idelta_x,rv_min,rv_max,traj,ksi,tw,num_parts,len_traj,traj_out)
    
    IMPLICIT NONE
    INTEGER,INTENT(IN)::opt_range,opt,ibins,tw,len_traj,num_parts,rv_min,rv_max
    DOUBLE PRECISION,INTENT(IN)::idelta_x,imin,imax
    DOUBLE PRECISION,INTENT(IN)::ksi
    DOUBLE PRECISION,DIMENSION(len_traj,num_parts,1),INTENT(IN)::traj
    INTEGER,DIMENSION(len_traj-2*tw,num_parts,1),INTENT(OUT)::traj_out
    
    INTEGER::bins
    DOUBLE PRECISION::delta_x,min,max,sobra
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::cumul
    INTEGER::Ltw,Ltw1,num_rep
    DOUBLE PRECISION::dsm
    DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::cumul_list,aux_list
    LOGICAL,DIMENSION(:),ALLOCATABLE::filter
    DOUBLE PRECISION::dist
    INTEGER::ii,jj,gg,kk,tt,lll,nn
    LOGICAL::switch,filt_min,filt_max
    
    filt_min=.FALSE.
    filt_max=.FALSE.
    
!!! For histogram:
    
    bins=ibins
    max=imax
    min=imin
    delta_x=idelta_x
    
    IF (opt_range==0) THEN
       IF (opt==1) THEN
          bins=CEILING((max-min)/delta_x)
          sobra=(bins*delta_x-(max-min))/2.0d0
          bins=bins+1
          min=min-sobra
          max=max+sobra
       ELSE
          delta_x=(max-min)/(bins*1.0d0)
       END IF
    ELSE
       IF (opt==1) THEN
          bins=CEILING((max-min)/delta_x)
       ELSE
          delta_x=(max-min)/(bins*1.0d0)
       END IF
    END IF
    
    IF (rv_min==1) THEN
       bins=bins+1
       filt_min=.TRUE.
    END IF
    IF (rv_max==1) THEN
       bins=bins+1
       filt_max=.TRUE.
    END IF
    
    !!
    ALLOCATE(cumul(bins))
    cumul=0.0d0
    
    Ltw=(2*tw+1)
    Ltw1=Ltw-1
    
    traj_out=0
    num_rep=0
    dsm=sqrt(2.0d0/(1.0d0*Ltw))*ksi !Kolmogorov-Smirnov
    
    DO nn=1,num_parts
       DO ii=1,len_traj-2*tw
          cumul=0.0d0
          IF ((filt_min.eqv..TRUE.).OR.(filt_max.eqv..TRUE.)) THEN
             DO kk=ii,ii+Ltw1
                tt=CEILING((traj(kk,nn,1)-min)/delta_x)
                IF (filt_min.eqv..TRUE.) THEN
                   IF (tt<1) THEN
                      tt=1
                   ELSE
                      tt=tt+1
                   END IF
                END IF
                IF (tt>bins) tt=bins
                cumul(tt)=cumul(tt)+1.0d0
             END DO
          ELSE
             DO kk=ii,ii+Ltw1
                tt=CEILING((traj(kk,nn,1)-min)/delta_x)
                IF (tt==0) tt=1
                cumul(tt)=cumul(tt)+1.0d0
             END DO
          END IF
          cumul=cumul/(1.0d0*Ltw)
          DO kk=1,bins-1
             cumul(kk+1)=cumul(kk+1)+cumul(kk)
          END DO
          switch=.true.
          DO jj=ii-1,1,-1
             gg=traj_out(jj,nn,1)
             IF (filter(gg).eqv..true.) THEN
                dist=MAXVAL(abs(cumul_list(gg,:)-cumul(:)),DIM=1)
                IF (dist<=dsm) THEN
                   traj_out(ii,nn,1)=gg
                   switch=.false.
                   EXIT
                END IF
                filter(gg)=.false.
                IF (COUNT(filter,DIM=1)==0) THEN
                   EXIT
                END IF
             END IF
          END DO
          IF (switch.eqv..true.) THEN
             DO jj=num_rep,1,-1
                IF (COUNT(filter,DIM=1)==0) THEN
                   EXIT
                END IF
                IF (filter(jj).eqv..true.) THEN
                   dist=MAXVAL(abs(cumul_list(jj,:)-cumul(:)),DIM=1)
                   IF (dist<=dsm) THEN
                      traj_out(ii,nn,1)=jj
                      switch=.false.
                      EXIT
                   END IF
                   filter(jj)=.false.
                END IF
             END DO
          END IF
          
          IF (switch.eqv..true.) THEN
             IF (num_rep>0) THEN
                ALLOCATE(aux_list(num_rep,bins))
                aux_list=cumul_list
                DEALLOCATE(cumul_list,filter)
                num_rep=num_rep+1
                ALLOCATE(cumul_list(num_rep,bins),filter(num_rep))
                filter=.true.
                cumul_list(1:(num_rep-1),:)=aux_list
                DEALLOCATE(aux_list)
                cumul_list(num_rep,:)=cumul
                traj_out(ii,nn,1)=num_rep
             ELSE
                num_rep=1
                ALLOCATE(cumul_list(num_rep,bins),filter(num_rep))
                filter=.true.
                cumul_list(num_rep,:)=cumul
                traj_out(ii,nn,1)=num_rep
             END IF
          ELSE
             filter=.true.
          END IF
       END DO
    END DO
    
    traj_out=traj_out-1
    
    DEALLOCATE(cumul_list,filter,cumul)
    
    
  END subroutine ganna


  subroutine free_distrib ()
    
    IMPLICIT NONE
    
    if (ALLOCATED(distrib)) DEALLOCATE(distrib)
    if (ALLOCATED(distrib_x)) DEALLOCATE(distrib_x)
    
  end subroutine free_distrib
  
  subroutine life_time_dist (opt_norm,traj,state,num_frames,num_parts,dims,num_states,mean)
    
    IMPLICIT NONE
    INTEGER,INTENT(IN):: opt_norm,num_parts,dims,num_frames,num_states
    DOUBLE PRECISION,DIMENSION(num_frames,num_parts,dims),INTENT(IN):: traj
    DOUBLE PRECISION,DIMENSION(num_states),INTENT(IN):: state
    DOUBLE PRECISION,INTENT(OUT):: mean
    
    INTEGER:: ii,jj,kk,ll,gg,contador,kkk,lll
    LOGICAL::inside_old,inside
    INTEGER:: contador_total
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::distrib_aux
    
    gg=100
    contador=0
    contador_total=0
    mean=0.0d0
    
    IF (ALLOCATED(distrib)) DEALLOCATE(distrib)
    IF (ALLOCATED(distrib_x)) DEALLOCATE(distrib_x)
    
    ALLOCATE(distrib(gg))
    distrib=0.0d0
    
    DO kkk=1,num_parts
       DO lll=1,dims
          contador=0
          inside_old=.false.
          DO ii=1,num_frames
             inside=.false.
             DO jj=1,num_states
                IF (traj(ii,kkk,lll)==state(jj)) THEN
                   inside=.true.
                   contador=contador+1
                   EXIT
                END IF
             END DO
             
             IF ((inside_old.eqv..true.).and.(inside.eqv..false.)) THEN
                IF (contador>gg) THEN
                   ALLOCATE(distrib_aux(gg))
                   distrib_aux(:)=distrib(:)
                   DEALLOCATE(distrib)
                   ALLOCATE(distrib(contador))
                   distrib(:gg)=distrib_aux(:)
                   distrib((gg+1):)=0.0d0
                   gg=contador
                   DEALLOCATE(distrib_aux)
                END IF
                distrib(contador)=distrib(contador)+1.0d0
                contador_total=contador_total+contador
                mean=mean+1.0d0
                contador=0
             END IF
             
             inside_old=inside
             
          END DO
       END DO
    END DO
    
    mean=(contador_total*1.0d0)/mean
    
    IF (opt_norm==1) THEN
       distrib(:)=distrib(:)/(contador_total*1.0d0)
    END IF
    
    jj=0
    DO ii=1,gg
       IF (distrib(ii)>0.0d0) THEN
          jj=jj+1
       END IF
    END DO
    
    ALLOCATE(distrib_aux(jj),distrib_x(jj))
    jj=0
    DO ii=1,gg
       IF (distrib(ii)>0.0d0) THEN
          jj=jj+1
          distrib_aux(jj)=distrib(ii)
          distrib_x(jj)=ii
       END IF
    END DO
    DEALLOCATE(distrib)
    ALLOCATE(distrib(jj))
    distrib=distrib_aux
    DEALLOCATE(distrib_aux)
    
  end subroutine life_time_dist
  
  subroutine fpt_dist (opt_norm,opt_from_state,opt_from_segment,opt_to_state,opt_to_segment, &
       from_state,from_segment,to_state,to_segment,traj,num_frames,num_parts,dims,from_num_states,to_num_states,mean)
    
    
    IMPLICIT NONE
    INTEGER,INTENT(IN):: opt_norm,opt_from_state,opt_from_segment,opt_to_state,opt_to_segment
    INTEGER,INTENT(IN):: num_parts,dims,num_frames
    INTEGER,INTENT(IN):: from_num_states,to_num_states
    DOUBLE PRECISION,DIMENSION(num_frames,num_parts,dims),INTENT(IN):: traj
    DOUBLE PRECISION,DIMENSION(from_num_states),INTENT(IN):: from_state
    DOUBLE PRECISION,DIMENSION(to_num_states),INTENT(IN):: to_state
    DOUBLE PRECISION,DIMENSION(2),INTENT(IN)::from_segment,to_segment
    DOUBLE PRECISION,INTENT(OUT):: mean
    
    INTEGER:: ii,jj,kk,ll,gg,contador,kkk,lll
    LOGICAL::entro,inside_to,inside_from
    INTEGER(KIND=8):: contador_total
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::distrib_aux
    DOUBLE PRECISION::pos
    
    gg=100
    contador=0
    contador_total=0
    mean=0.0d0
    
    IF (ALLOCATED(distrib)) DEALLOCATE(distrib)
    IF (ALLOCATED(distrib_x)) DEALLOCATE(distrib_x)
    
    ALLOCATE(distrib(gg))
    distrib=0.0d0

    DO kkk=1,num_parts
       DO lll=1,dims
          entro=.false.
          contador=0
          DO ii=num_frames,1,-1
             pos=traj(ii,kkk,lll)
             
             inside_to=.false.
             IF (opt_to_segment==1) THEN
                IF ((to_segment(1)<pos).and.(pos<to_segment(2))) THEN
                   inside_to=.true.
                   entro=.true.
                END IF
             ELSE
                DO jj=1,to_num_states
                   IF (pos==to_state(jj)) THEN
                      inside_to=.true.
                      entro=.true.
                      EXIT
                   END IF
                END DO
             END IF
             
             inside_from=.false.
             IF (inside_to.eqv..false.) THEN
                IF (opt_from_segment==1) THEN
                   IF ((from_segment(1)<pos).and.(pos<from_segment(2))) THEN
                      inside_from=.true.
                   END IF
                ELSE IF (opt_from_state==1) THEN
                   DO jj=1,from_num_states
                      IF (pos==from_state(jj)) THEN
                         inside_from=.true.
                         EXIT
                      END IF
                   END DO
                ELSE
                   inside_from=.true.
                END IF
             END IF
             
             IF ((inside_to.eqv..false.).and.(entro.eqv..true.)) THEN
                contador=contador+1
             ELSE
                contador=0
             END IF
             
             IF ((inside_from.eqv..true.).and.(entro.eqv..true.)) THEN
                IF (contador>gg) THEN
                   ALLOCATE(distrib_aux(gg))
                   distrib_aux(:)=distrib(:)
                   DEALLOCATE(distrib)
                   ALLOCATE(distrib(contador))
                   distrib(:gg)=distrib_aux(:)
                   distrib((gg+1):)=0.0d0
                   gg=contador
                   DEALLOCATE(distrib_aux)
                END IF
                distrib(contador)=distrib(contador)+1.0d0
                contador_total=contador_total+contador
                mean=mean+1.0d0
             END IF
             
          END DO
       END DO
    END DO
    
    IF (mean>0.0d0) THEN
       mean=(contador_total*1.0d0)/mean
    ELSE
       mean=0.0d0
    END IF
    
    IF (opt_norm==1) THEN
       IF (contador_total==0) THEN
          distrib=0.0d0
       ELSE
          distrib(:)=distrib(:)/(contador_total*1.0d0)
       END IF
    END IF
    
    jj=0
    DO ii=1,gg
       IF (distrib(ii)>0.0d0) THEN
          jj=jj+1
       END IF
    END DO
    
    IF (jj/=0) THEN
       ALLOCATE(distrib_aux(jj),distrib_x(jj))
       jj=0
       DO ii=1,gg
          IF (distrib(ii)>0.0d0) THEN
             jj=jj+1
             distrib_aux(jj)=distrib(ii)
             distrib_x(jj)=ii
          END IF
       END DO
       DEALLOCATE(distrib)
       ALLOCATE(distrib(jj))
       distrib=distrib_aux
       DEALLOCATE(distrib_aux)
    ELSE
       DEALLOCATE(distrib)
       ALLOCATE(distrib(1),distrib_x(1))
       distrib=0.0d0
       distrib_x=0.0d0
    END IF
    
  end subroutine fpt_dist
  
  
  subroutine tt_dist (opt_norm,opt_noreturn,opt_from_state,opt_from_segment,opt_to_state,opt_to_segment, &
       from_state,from_segment,to_state,to_segment,traj,num_frames,num_parts,dims,from_num_states,to_num_states,mean)
    
    
    IMPLICIT NONE
    INTEGER,INTENT(IN):: opt_norm,opt_from_state,opt_from_segment,opt_to_state,opt_to_segment,opt_noreturn
    INTEGER,INTENT(IN):: num_parts,dims,num_frames
    INTEGER,INTENT(IN):: from_num_states,to_num_states
    DOUBLE PRECISION,DIMENSION(num_frames,num_parts,dims),INTENT(IN):: traj
    DOUBLE PRECISION,DIMENSION(from_num_states),INTENT(IN):: from_state
    DOUBLE PRECISION,DIMENSION(to_num_states),INTENT(IN):: to_state
    DOUBLE PRECISION,DIMENSION(2),INTENT(IN)::from_segment,to_segment
    DOUBLE PRECISION,INTENT(OUT):: mean

    INTEGER:: ii,jj,kk,ll,gg,contador,kkk,lll
    LOGICAL::entro,inside_to,inside_from,last_from
    INTEGER(KIND=8):: contador_total
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::distrib_aux
    DOUBLE PRECISION::pos
    
    gg=100
    contador=0
    contador_total=0
    mean=0.0d0
    
    IF (ALLOCATED(distrib)) DEALLOCATE(distrib)
    IF (ALLOCATED(distrib_x)) DEALLOCATE(distrib_x)
    
    ALLOCATE(distrib(gg))
    distrib=0.0d0
    
    DO kkk=1,num_parts
       DO lll=1,dims
          entro=.false.
          last_from=.false.
          contador=0
          DO ii=num_frames,1,-1
             pos=traj(ii,kkk,lll)
             
             inside_to=.false.
             IF (opt_to_segment==1) THEN
                IF ((to_segment(1)<pos).and.(pos<to_segment(2))) THEN
                   inside_to=.true.
                   entro=.true.
                END IF
             ELSE
                DO jj=1,to_num_states
                   IF (pos==to_state(jj)) THEN
                      inside_to=.true.
                      entro=.true.
                      EXIT
                   END IF
                END DO
             END IF
             
             inside_from=.false.
             IF (inside_to.eqv..false.) THEN
                IF (opt_from_segment==1) THEN
                   IF ((from_segment(1)<pos).and.(pos<from_segment(2))) THEN
                      inside_from=.true.
                   END IF
                ELSE IF (opt_from_state==1) THEN
                   DO jj=1,from_num_states
                      IF (pos==from_state(jj)) THEN
                         inside_from=.true.
                         EXIT
                      END IF
                   END DO
                ELSE
                   inside_from=.true.
                END IF
             END IF
             
             IF ((inside_to.eqv..false.).and.(entro.eqv..true.)) THEN
                contador=contador+1
             ELSE
                contador=0
             END IF
             
             IF ((inside_from.eqv..true.).and.(entro.eqv..true.).and.(last_from.eqv..false.)) THEN
                
                IF (contador>gg) THEN
                   ALLOCATE(distrib_aux(gg))
                   distrib_aux(:)=distrib(:)
                   DEALLOCATE(distrib)
                   ALLOCATE(distrib(contador))
                   distrib(:gg)=distrib_aux(:)
                   distrib((gg+1):)=0.0d0
                   gg=contador
                   DEALLOCATE(distrib_aux)
                END IF
                distrib(contador)=distrib(contador)+1.0d0
                contador_total=contador_total+contador
                mean=mean+1.0d0
                
                IF (opt_noreturn==1) THEN
                   entro=.false.
                   contador=0
                END IF
                
             END IF
             
             last_from=inside_from
             
          END DO
       END DO
    END DO
    
    IF (mean>0.0d0) THEN
       mean=(contador_total*1.0d0)/mean
    ELSE
       mean=0.0d0
    END IF
    
    IF (opt_norm==1) THEN
       IF (contador_total==0) THEN
          distrib=0.0d0
       ELSE
          distrib(:)=distrib(:)/(contador_total*1.0d0)
       END IF
    END IF
    
    jj=0
    DO ii=1,gg
       IF (distrib(ii)>0.0d0) THEN
          jj=jj+1
       END IF
    END DO
    
    IF (jj/=0) THEN
       ALLOCATE(distrib_aux(jj),distrib_x(jj))
       jj=0
       DO ii=1,gg
          IF (distrib(ii)>0.0d0) THEN
             jj=jj+1
             distrib_aux(jj)=distrib(ii)
             distrib_x(jj)=ii
          END IF
       END DO
       DEALLOCATE(distrib)
       ALLOCATE(distrib(jj))
       distrib=distrib_aux
       DEALLOCATE(distrib_aux)
    ELSE
       DEALLOCATE(distrib)
       ALLOCATE(distrib(1),distrib_x(1))
       distrib=0.0d0
       distrib_x=0.0d0
    END IF
    
  end subroutine tt_dist


subroutine fcpt_dist (opt_norm,opt_noreturn,opt_states,opt_segments, &
     states,segments,commitment,traj,num_frames,num_parts,dims,num_states,num_segments,num_commits,mean)


  IMPLICIT NONE
  INTEGER,INTENT(IN):: opt_norm,opt_states,opt_segments,opt_noreturn
  INTEGER,INTENT(IN):: num_parts,dims,num_frames
  INTEGER,INTENT(IN):: num_states,num_segments,num_commits
  DOUBLE PRECISION,DIMENSION(num_frames,num_parts,dims),INTENT(IN):: traj
  DOUBLE PRECISION,DIMENSION(num_states),INTENT(IN):: states
  DOUBLE PRECISION,DIMENSION(num_segments,2),INTENT(IN)::segments
  INTEGER,DIMENSION(num_commits),INTENT(IN)::commitment
  DOUBLE PRECISION,INTENT(OUT):: mean

  INTEGER:: ii,jj,kk,ll,gg,contador,kkk,lll,toca
  INTEGER,DIMENSION(num_commits)::visited
  LOGICAL::entro,inside_to,inside_false
  LOGICAL:: hecho,listo,touch
  INTEGER(KIND=8):: contador_total
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::distrib_aux
  DOUBLE PRECISION::pos

  gg=100
  contador=0
  contador_total=0
  mean=0.0d0
  visited=0
  toca=0

  IF (ALLOCATED(distrib)) DEALLOCATE(distrib)
  IF (ALLOCATED(distrib_x)) DEALLOCATE(distrib_x)

  ALLOCATE(distrib(gg))
  distrib=0.0d0

  DO kkk=1,num_parts
     DO lll=1,dims
        entro=.false.
        touch=.false.
        contador=0
        DO ii=num_frames,1,-1

           pos=traj(ii,kkk,lll)

           inside_to=.false.
           listo=.false.
           IF (opt_segments==1) THEN
              IF ((segments(num_segments,1)<pos).and.(pos<segments(num_segments,2))) THEN
                 inside_to=.true.
                 entro=.true.
                 touch=.false.
                 visited=0
                 visited(num_segments)=1
                 toca=num_segments-1
              END IF
           ELSE
              IF (pos==states(num_states)) THEN
                 inside_to=.true.
                 entro=.true.
                 touch=.false.
                 visited=0
                 visited(num_states)=1
                 toca=num_states-1
              END IF
           END IF

           listo=.false.
           IF ((entro.eqv..true.).and.(inside_to.eqv..false.)) THEN
              IF (opt_segments==1) THEN
                 !IF ((from_segment(1)<pos).and.(pos<from_segment(2))) THEN
                 !   !!...
                 !END IF
              ELSE 
                 IF ((commitment(toca+1)==1).and.(pos==states(toca+1))) THEN
                    IF (visited(1)/=1) THEN
                       listo=.true.
                    END IF
                 END IF
                 IF (listo.eqv..false.) THEN
                    IF ((pos==states(toca)).and.(commitment(toca)==1)) THEN
                       visited(toca)=1
                       toca=toca-1
                       IF (toca==0) toca=1
                       listo=.true.
                    END IF
                    IF ((pos/=states(toca)).and.(commitment(toca)==0)) THEN
                       visited(toca)=0
                       toca=toca-1
                       IF (toca==0) toca=1
                       IF (pos==states(toca).and.(commitment(toca)==1)) THEN ! traj=[1,1,2] for states [1,0,2][T,F,T]
                          visited(toca)=1
                          toca=toca-1
                          IF (toca==0) toca=1
                       END IF
                       listo=.true.
                    END IF
                 END IF
                 IF (listo.eqv..false.) THEN
                    IF (pos/=states(toca+1).and.commitment(toca+1)==0) THEN
                       visited(1)=0   ! [1,3,1,2] with [1,0,2][T,F,T]
                       listo=.true.
                    END IF
                 END IF
                 IF (listo.eqv..false.) THEN
                    contador=0
                    visited=0
                    touch=.false.
                    entro=.false.
                 END IF
              END IF
           END IF
           
           IF (pos==states(1)) THEN
              touch=.true.
           END IF
           IF (opt_noreturn==1) THEN
              IF ((touch.eqv..true.).and.(pos/=states(1))) THEN
                 listo=.false.
                 contador=0
                 visited=0
                 touch=.false.
                 entro=.false.
              END IF
           END IF

           IF (listo.eqv..true.) THEN
              contador=contador+1
              IF (visited(1)==1) THEN
                 HECHO=.true.
                 IF (opt_states==1) THEN
                    DO kk=num_states,1,-1
                       IF (visited(kk)/=commitment(kk)) THEN
                          HECHO=.false.
                          EXIT
                       END IF
                    END DO
                 ELSE
                    DO kk=num_segments,1,-1
                       IF (visited(kk)/=commitment(kk)) THEN
                          HECHO=.false.
                          EXIT
                       END IF
                    END DO
                 END IF
                 IF (HECHO.eqv..true.) THEN
                    IF (contador>gg) THEN
                       ALLOCATE(distrib_aux(gg))
                       distrib_aux(:)=distrib(:)
                       DEALLOCATE(distrib)
                       ALLOCATE(distrib(contador))
                       distrib(:gg)=distrib_aux(:)
                       distrib((gg+1):)=0.0d0
                       gg=contador
                       DEALLOCATE(distrib_aux)
                    END IF
                    distrib(contador)=distrib(contador)+1.0d0
                    contador_total=contador_total+contador
                    mean=mean+1.0d0
                 END IF
              END IF
           ELSE
              contador=0
           END IF

        END DO
     END DO
  END DO
  
  IF (mean>0.0d0) THEN
     mean=(contador_total*1.0d0)/mean
  ELSE
     mean=0.0d0
  END IF

  IF (opt_norm==1) THEN
     IF (contador_total==0) THEN
        distrib=0.0d0
     ELSE
        distrib(:)=distrib(:)/(contador_total*1.0d0)
     END IF
  END IF

  jj=0
  DO ii=1,gg
     IF (distrib(ii)>0.0d0) THEN
        jj=jj+1
     END IF
  END DO

  IF (jj/=0) THEN
     ALLOCATE(distrib_aux(jj),distrib_x(jj))
     jj=0
     DO ii=1,gg
        IF (distrib(ii)>0.0d0) THEN
           jj=jj+1
           distrib_aux(jj)=distrib(ii)
           distrib_x(jj)=ii
        END IF
     END DO
     DEALLOCATE(distrib)
     ALLOCATE(distrib(jj))
     distrib=distrib_aux
     DEALLOCATE(distrib_aux)
  ELSE
     DEALLOCATE(distrib)
     ALLOCATE(distrib(1),distrib_x(1))
     distrib=0.0d0
     distrib_x=0.0d0
  END IF


end subroutine fcpt_dist


subroutine ctt_dist (opt_norm,opt_noreturn,opt_states,opt_segments, &
     states,segments,commitment,traj,num_frames,num_parts,dims,num_states,num_segments,num_commits,mean)


  IMPLICIT NONE
  INTEGER,INTENT(IN):: opt_norm,opt_states,opt_segments,opt_noreturn
  INTEGER,INTENT(IN):: num_parts,dims,num_frames
  INTEGER,INTENT(IN):: num_states,num_segments,num_commits
  DOUBLE PRECISION,DIMENSION(num_frames,num_parts,dims),INTENT(IN):: traj
  DOUBLE PRECISION,DIMENSION(num_states),INTENT(IN):: states
  DOUBLE PRECISION,DIMENSION(num_segments,2),INTENT(IN)::segments
  INTEGER,DIMENSION(num_commits),INTENT(IN)::commitment
  DOUBLE PRECISION,INTENT(OUT):: mean

  INTEGER:: ii,jj,kk,ll,gg,contador,kkk,lll,toca
  INTEGER,DIMENSION(num_commits)::visited
  LOGICAL::entro,inside_to,inside_false
  LOGICAL:: hecho,listo,touch,last_in
  INTEGER(KIND=8):: contador_total
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::distrib_aux
  DOUBLE PRECISION::pos

  gg=100
  contador=0
  contador_total=0
  mean=0.0d0
  visited=0
  toca=0

  IF (ALLOCATED(distrib)) DEALLOCATE(distrib)
  IF (ALLOCATED(distrib_x)) DEALLOCATE(distrib_x)

  ALLOCATE(distrib(gg))
  distrib=0.0d0

  DO kkk=1,num_parts
     DO lll=1,dims
        entro=.false.
        touch=.false.
        contador=0
        last_in=.false.
        DO ii=num_frames,1,-1

           pos=traj(ii,kkk,lll)

           inside_to=.false.
           listo=.false.
           IF (opt_segments==1) THEN
              IF ((segments(num_segments,1)<pos).and.(pos<segments(num_segments,2))) THEN
                 inside_to=.true.
                 entro=.true.
                 touch=.false.
                 visited=0
                 visited(num_segments)=1
                 toca=num_segments-1
              END IF
           ELSE
              IF (pos==states(num_states)) THEN
                 inside_to=.true.
                 entro=.true.
                 touch=.false.
                 visited=0
                 visited(num_states)=1
                 toca=num_states-1
              END IF
           END IF

           listo=.false.
           IF ((entro.eqv..true.).and.(inside_to.eqv..false.)) THEN
              IF (opt_segments==1) THEN
                 !IF ((from_segment(1)<pos).and.(pos<from_segment(2))) THEN
                 !   !!...
                 !END IF
              ELSE 
                 IF ((commitment(toca+1)==1).and.(pos==states(toca+1))) THEN
                    IF (visited(1)/=1) THEN
                       listo=.true.
                    END IF
                 END IF
                 IF (listo.eqv..false.) THEN
                    IF ((pos==states(toca)).and.(commitment(toca)==1)) THEN
                       visited(toca)=1
                       toca=toca-1
                       IF (toca==0) toca=1
                       listo=.true.
                    END IF
                    IF ((pos/=states(toca)).and.(commitment(toca)==0)) THEN
                       visited(toca)=0
                       toca=toca-1
                       IF (toca==0) toca=1
                       IF (pos==states(toca).and.(commitment(toca)==1)) THEN ! traj=[1,1,2] for states [1,0,2][T,F,T]
                          visited(toca)=1
                          toca=toca-1
                          IF (toca==0) toca=1
                       END IF
                       listo=.true.
                    END IF
                 END IF
                 IF (listo.eqv..false.) THEN
                    IF (pos/=states(toca+1).and.commitment(toca+1)==0) THEN
                       visited(1)=0   ! [1,3,1,2] with [1,0,2][T,F,T]
                       listo=.true.
                    END IF
                 END IF
                 IF (listo.eqv..false.) THEN
                    contador=0
                    visited=0
                    touch=.false.
                    entro=.false.
                 END IF
              END IF
           END IF
           
           IF (pos==states(1)) THEN
              touch=.true.
           END IF
           IF (opt_noreturn==1) THEN
              IF ((touch.eqv..true.).and.(pos/=states(1))) THEN
                 listo=.false.
                 contador=0
                 visited=0
                 touch=.false.
                 entro=.false.
              END IF
           END IF

           IF (listo.eqv..true.) THEN
              contador=contador+1
              IF ((visited(1)==1).and.(last_in.eqv..false.)) THEN
                 HECHO=.true.
                 IF (opt_states==1) THEN
                    DO kk=num_states,1,-1
                       IF (visited(kk)/=commitment(kk)) THEN
                          HECHO=.false.
                          EXIT
                       END IF
                    END DO
                 ELSE
                    DO kk=num_segments,1,-1
                       IF (visited(kk)/=commitment(kk)) THEN
                          HECHO=.false.
                          EXIT
                       END IF
                    END DO
                 END IF
                 IF (HECHO.eqv..true.) THEN
                    IF (contador>gg) THEN
                       ALLOCATE(distrib_aux(gg))
                       distrib_aux(:)=distrib(:)
                       DEALLOCATE(distrib)
                       ALLOCATE(distrib(contador))
                       distrib(:gg)=distrib_aux(:)
                       distrib((gg+1):)=0.0d0
                       gg=contador
                       DEALLOCATE(distrib_aux)
                    END IF
                    distrib(contador)=distrib(contador)+1.0d0
                    contador_total=contador_total+contador
                    mean=mean+1.0d0
                 END IF
              END IF
           ELSE
              contador=0
           END IF
           
           last_in=touch

        END DO
     END DO
  END DO
  
  IF (mean>0.0d0) THEN
     mean=(contador_total*1.0d0)/mean
  ELSE
     mean=0.0d0
  END IF

  IF (opt_norm==1) THEN
     IF (contador_total==0) THEN
        distrib=0.0d0
     ELSE
        distrib(:)=distrib(:)/(contador_total*1.0d0)
     END IF
  END IF

  jj=0
  DO ii=1,gg
     IF (distrib(ii)>0.0d0) THEN
        jj=jj+1
     END IF
  END DO

  IF (jj/=0) THEN
     ALLOCATE(distrib_aux(jj),distrib_x(jj))
     jj=0
     DO ii=1,gg
        IF (distrib(ii)>0.0d0) THEN
           jj=jj+1
           distrib_aux(jj)=distrib(ii)
           distrib_x(jj)=ii
        END IF
     END DO
     DEALLOCATE(distrib)
     ALLOCATE(distrib(jj))
     distrib=distrib_aux
     DEALLOCATE(distrib_aux)
  ELSE
     DEALLOCATE(distrib)
     ALLOCATE(distrib(1),distrib_x(1))
     distrib=0.0d0
     distrib_x=0.0d0
  END IF


end subroutine ctt_dist

END MODULE GLOB
