!!!#################################
!!!####  COMMON VARIABLES AND
!!!####  SUBROUTINES TO UPLOAD THEM
!!!#################################

MODULE GLOB

INTEGER::num_atoms,num_nodes,num_categories,num_sets
INTEGER,DIMENSION(:),ALLOCATABLE::node_run_ats,atom2node,trad2py_node,trad2py_atom,atomspernode
INTEGER,DIMENSION(:),ALLOCATABLE::node_category,node_hbs_num,node_bs_num
INTEGER,DIMENSION(:),ALLOCATABLE::atom_hbs_num,atom_bs_num
INTEGER,DIMENSION(:),ALLOCATABLE::symm_ats_crits,symm_ats_start,symm_ats

INTEGER,DIMENSION(:),ALLOCATABLE::vecti_aux

INTEGER::T_num_hbs,T_num_bs
INTEGER,DIMENSION(:),ALLOCATABLE::T_hbs_start,T_bs_start,T_hbs_ind,T_bs_ind

INTEGER,DIMENSION(:),ALLOCATABLE::mss_ind_atoms,mss_ind_nodes,mss_symm,mss

CONTAINS


SUBROUTINE load_topol(inode_run_ats,iatom2node,itrad2py_node,itrad2py_atom,&
     isymm_ats_start,isymm_ats_crits,isymm_ats,inode_category,inum_nodes,inum_atoms,isymm_ats_dim)

  IMPLICIT NONE
  INTEGER,INTENT(IN)::inum_nodes,inum_atoms,isymm_ats_dim
  INTEGER,DIMENSION(inum_nodes+1),INTENT(IN)::inode_run_ats
  INTEGER,DIMENSION(inum_atoms),INTENT(IN)::iatom2node
  INTEGER,DIMENSION(inum_nodes),INTENT(IN)::itrad2py_node
  INTEGER,DIMENSION(inum_atoms),INTENT(IN)::itrad2py_atom
  INTEGER,DIMENSION(inum_nodes+1),INTENT(IN)::isymm_ats_start
  INTEGER,DIMENSION(inum_nodes),INTENT(IN)::isymm_ats_crits
  INTEGER,DIMENSION(isymm_ats_dim),INTENT(IN)::isymm_ats
  INTEGER,DIMENSION(inum_nodes),INTENT(IN)::inode_category

  INTEGER::ii

  num_nodes=inum_nodes
  num_atoms=inum_atoms

  IF (ALLOCATED(node_run_ats)) DEALLOCATE(node_run_ats)
  IF (ALLOCATED(atom2node))    DEALLOCATE(atom2node)
  IF (ALLOCATED(trad2py_node)) DEALLOCATE(trad2py_node)
  IF (ALLOCATED(trad2py_atom)) DEALLOCATE(trad2py_atom)
  IF (ALLOCATED(atomspernode)) DEALLOCATE(atomspernode)

  ALLOCATE(node_run_ats(num_nodes+1),atom2node(num_atoms))
  node_run_ats(:)=inode_run_ats(:)
  atom2node(:)=iatom2node(:)

  ALLOCATE(trad2py_node(num_nodes),trad2py_atom(num_atoms))
  trad2py_node(:)=itrad2py_node(:)
  trad2py_atom(:)=itrad2py_atom(:)

  ALLOCATE(atomspernode(num_nodes))

  DO ii=1,num_nodes
     atomspernode(ii)=node_run_ats(ii+1)-node_run_ats(ii)
  END DO

  IF (ALLOCATED(symm_ats)) DEALLOCATE(symm_ats)
  IF (ALLOCATED(symm_ats_crits)) DEALLOCATE(symm_ats_crits)
  IF (ALLOCATED(symm_ats_start)) DEALLOCATE(symm_ats_start)

  ALLOCATE(symm_ats(isymm_ats_dim),symm_ats_crits(num_nodes),symm_ats_start(num_nodes+1))
  symm_ats(:)=isymm_ats(:)
  symm_ats_start(:)=isymm_ats_start(:)
  symm_ats_crits(:)=isymm_ats_crits(:)

  IF (ALLOCATED(node_category)) DEALLOCATE(node_category)
  IF (ALLOCATED(node_hbs_num)) DEALLOCATE(node_hbs_num)
  IF (ALLOCATED(node_bs_num)) DEALLOCATE(node_bs_num)
  IF (ALLOCATED(atom_hbs_num))    DEALLOCATE(atom_hbs_num)
  IF (ALLOCATED(atom_bs_num))     DEALLOCATE(atom_bs_num)

  ALLOCATE(node_category(num_nodes))
  ALLOCATE(node_hbs_num(num_nodes),node_bs_num(num_nodes))
  ALLOCATE(atom_hbs_num(num_atoms),atom_bs_num(num_atoms))

  node_category(:)=inode_category(:)

END SUBROUTINE load_topol

SUBROUTINE load_net(iT_hbs_start,iT_bs_start,iT_hbs_ind,iT_bs_ind,iT_num_hbs,iT_num_bs,inum_atoms)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::inum_atoms,iT_num_hbs,iT_num_bs
  INTEGER,DIMENSION(inum_atoms+1),INTENT(IN)::iT_hbs_start,iT_bs_start
  INTEGER,DIMENSION(iT_num_hbs),INTENT(IN)::iT_hbs_ind
  INTEGER,DIMENSION(iT_num_bs),INTENT(IN)::iT_bs_ind

  INTEGER::ii,jj,kk

  IF (ALLOCATED(T_hbs_start))  DEALLOCATE(T_hbs_start)
  IF (ALLOCATED(T_bs_start))   DEALLOCATE(T_bs_start)
  IF (ALLOCATED(T_hbs_ind))    DEALLOCATE(T_hbs_ind)
  IF (ALLOCATED(T_bs_ind))     DEALLOCATE(T_bs_ind)

  T_num_hbs=iT_num_hbs
  T_num_bs=iT_num_bs

  ALLOCATE(T_hbs_start(num_atoms+1),T_bs_start(num_atoms+1))
  ALLOCATE(T_hbs_ind(T_num_hbs),T_bs_ind(T_num_bs))

  T_hbs_start(:)=iT_hbs_start(:)
  T_bs_start(:)=iT_bs_start(:)
  T_hbs_ind(:)=iT_hbs_ind(:)
  T_bs_ind(:)=iT_bs_ind(:)

  node_hbs_num(:)=0
  node_bs_num(:)=0
  DO ii=1,num_atoms
     kk=atom2node(ii)
     jj=T_hbs_start(ii+1)-T_hbs_start(ii)
     atom_hbs_num(ii)=jj
     node_hbs_num(kk)=node_hbs_num(kk)+jj
     jj=T_bs_start(ii+1)-T_bs_start(ii)
     atom_bs_num(ii)=jj
     node_bs_num(kk)=node_bs_num(kk)+jj
  END DO

END SUBROUTINE load_net

!!!#################################
!!!####  GENERAL FUNCTIONS TO SORT
!!!#################################

!!$SUBROUTINE SORT_INT_ATS_1SH (num_ats,dim_vecti_aux,order,valores,interruptor,num_crits)
!!$
!!$  IMPLICIT NONE
!!$
!!$  INTEGER,INTENT(IN)::num_ats
!!$  INTEGER,INTENT(INOUT)::dim_vecti_aux
!!$  INTEGER,DIMENSION(num_ats),INTENT(INOUT)::order
!!$  INTEGER,DIMENSION(dim_vecti_aux),INTENT(IN)::valores
!!$  LOGICAL,INTENT(INOUT)::interruptor
!!$  INTEGER,INTENT(INOUT)::num_crits
!!$
!!$  INTEGER::ii,jj,kk,ll,gg,idim,new_num_crits,tope
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::val_aux,ind_aux,order_aux
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::vals,inds
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::new_symm_aux,cajon
!!$  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro
!!$
!!$  new_num_crits=0
!!$
!!$  gg=0
!!$  DO ii=1,num_crits
!!$     gg=gg+1
!!$     idim=vecti_aux(gg)
!!$     ALLOCATE(val_aux(idim),ind_aux(idim),order_aux(idim))
!!$     DO jj=1,idim
!!$        gg=gg+1
!!$        kk=vecti_aux(gg)
!!$        ll=order(kk)
!!$        ind_aux(jj)=kk
!!$        order_aux(jj)=ll
!!$        val_aux(jj)=valores(gg)
!!$     END DO
!!$     ALLOCATE(filtro(idim),vals(idim),inds(idim))
!!$     filtro=.TRUE.
!!$     DO jj=1,idim
!!$        kk=MAXLOC(val_aux,DIM=1,MASK=filtro(:))
!!$        filtro(kk)=.FALSE.
!!$        ll=ind_aux(jj)
!!$        inds(jj)=ll
!!$        vals(jj)=val_aux(kk)
!!$        order(ll)=order_aux(kk)
!!$     END DO
!!$     interruptor=.FALSE.
!!$     DO jj=2,idim
!!$        IF (vals(jj-1)==vals(jj)) THEN
!!$           IF (interruptor.eqv..FALSE.) THEN
!!$              order_aux(1)=inds(jj-1)
!!$              order_aux(2)=inds(jj)
!!$              ll=2
!!$              interruptor=.TRUE.
!!$           ELSE
!!$              ll=ll+1
!!$              order_aux(ll)=inds(jj)
!!$           END IF
!!$        ELSE
!!$           IF (interruptor.eqv..True.) THEN
!!$              interruptor=.FALSE.
!!$              new_num_crits=new_num_crits+1
!!$              IF (new_num_crits==1) THEN
!!$                 tope=ll+1
!!$                 ALLOCATE(new_symm_aux(tope))
!!$                 new_symm_aux(1)=ll
!!$                 new_symm_aux(2:tope)=order_aux(1:ll)
!!$              ELSE
!!$                 ALLOCATE(cajon(tope))
!!$                 cajon(:)=new_symm_aux(:)
!!$                 DEALLOCATE(new_symm_aux)
!!$                 ALLOCATE(new_symm_aux(tope+1+ll))
!!$                 new_symm_aux(1:tope)=cajon(:)
!!$                 tope=tope+1
!!$                 new_symm_aux(tope)=ll
!!$                 new_symm_aux((tope+1):(tope+ll))=order_aux(1:ll)
!!$                 tope=tope+ll
!!$                 DEALLOCATE(cajon)
!!$              END IF
!!$           END IF
!!$        END IF
!!$     END DO
!!$     IF (interruptor.eqv..True.) THEN
!!$        interruptor=.FALSE.
!!$        new_num_crits=new_num_crits+1
!!$        IF (new_num_crits==1) THEN
!!$           tope=ll+1
!!$           ALLOCATE(new_symm_aux(tope))
!!$           new_symm_aux(1)=ll
!!$           new_symm_aux(2:tope)=order_aux(1:ll)
!!$        ELSE
!!$           ALLOCATE(cajon(tope))
!!$           cajon(:)=new_symm_aux(:)
!!$           DEALLOCATE(new_symm_aux)
!!$           ALLOCATE(new_symm_aux(tope+1+ll))
!!$           new_symm_aux(1:tope)=cajon(:)
!!$           tope=tope+1
!!$           new_symm_aux(tope)=ll
!!$           new_symm_aux((tope+1):(tope+ll))=order_aux(1:ll)
!!$           tope=tope+ll
!!$           DEALLOCATE(cajon)
!!$        END IF
!!$     END IF
!!$     DEALLOCATE(val_aux,ind_aux,order_aux)
!!$     DEALLOCATE(filtro,vals,inds)
!!$  END DO
!!$
!!$  num_crits=new_num_crits
!!$  DEALLOCATE(vecti_aux)
!!$  IF (num_crits==0) THEN
!!$     interruptor=.FALSE.
!!$     dim_vecti_aux=0
!!$  ELSE
!!$     interruptor=.TRUE.
!!$     ALLOCATE(vecti_aux(tope))
!!$     dim_vecti_aux=tope
!!$     vecti_aux(:)=new_symm_aux(:)
!!$     DEALLOCATE(new_symm_aux)
!!$  END IF
!!$
!!$END SUBROUTINE SORT_INT_ATS_1SH
!!$
!!$SUBROUTINE SORT_INT_MATRIX_ATS_1SH (num_ats,dim_vecti_aux,dim_matrix,order,valores,interruptor,num_crits)
!!$
!!$  IMPLICIT NONE
!!$
!!$  INTEGER,INTENT(IN)::num_ats,dim_matrix
!!$  INTEGER,INTENT(INOUT)::dim_vecti_aux
!!$  INTEGER,DIMENSION(num_ats),INTENT(INOUT)::order
!!$  INTEGER,DIMENSION(dim_vecti_aux,dim_matrix),INTENT(IN)::valores
!!$  LOGICAL,INTENT(INOUT)::interruptor
!!$  INTEGER,INTENT(INOUT)::num_crits
!!$
!!$  INTEGER::ii,jj,kk,ll,gg,nn,pp,idim,new_num_crits,tope
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::val_aux,ind_aux,order_aux
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::vals,inds
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::new_symm_aux,cajon
!!$  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro
!!$  LOGICAL,DIMENSION(:,:),ALLOCATABLE::suprafiltro
!!$
!!$  ALLOCATE(suprafiltro(dim_vecti_aux,dim_matrix))
!!$  suprafiltro(:,:)=.TRUE.
!!$
!!$  pp=1
!!$  DO WHILE ((num_crits>0).AND.(pp<=dim_matrix))
!!$     
!!$     new_num_crits=0
!!$     
!!$     gg=0
!!$     DO ii=1,num_crits
!!$        gg=gg+1
!!$        idim=vecti_aux(gg)
!!$        ALLOCATE(val_aux(idim),ind_aux(idim),order_aux(idim))
!!$        DO jj=1,idim
!!$           gg=gg+1
!!$           kk=vecti_aux(gg)
!!$           ll=order(kk)
!!$           ind_aux(jj)=kk
!!$           order_aux(jj)=ll
!!$           nn=MAXLOC(valores(gg,:),DIM=1,MASK=suprafiltro(gg,:))
!!$           suprafiltro(gg,nn)=.FALSE.
!!$           val_aux(jj)=valores(gg,nn)
!!$        END DO
!!$        ALLOCATE(filtro(idim),vals(idim),inds(idim))
!!$        filtro=.TRUE.
!!$        DO jj=1,idim
!!$           kk=MAXLOC(val_aux,DIM=1,MASK=filtro(:))
!!$           filtro(kk)=.FALSE.
!!$           ll=ind_aux(jj)
!!$           inds(jj)=ll
!!$           vals(jj)=val_aux(kk)
!!$           order(ll)=order_aux(kk)
!!$        END DO
!!$        interruptor=.FALSE.
!!$        DO jj=2,idim
!!$           IF (vals(jj-1)==vals(jj)) THEN
!!$              IF (interruptor.eqv..FALSE.) THEN
!!$                 order_aux(1)=inds(jj-1)
!!$                 order_aux(2)=inds(jj)
!!$                 ll=2
!!$                 interruptor=.TRUE.
!!$              ELSE
!!$                 ll=ll+1
!!$                 order_aux(ll)=inds(jj)
!!$              END IF
!!$           ELSE
!!$              IF (interruptor.eqv..True.) THEN
!!$                 interruptor=.FALSE.
!!$                 new_num_crits=new_num_crits+1
!!$                 IF (new_num_crits==1) THEN
!!$                    tope=ll+1
!!$                    ALLOCATE(new_symm_aux(tope))
!!$                    new_symm_aux(1)=ll
!!$                    new_symm_aux(2:tope)=order_aux(1:ll)
!!$                 ELSE
!!$                    ALLOCATE(cajon(tope))
!!$                    cajon(:)=new_symm_aux(:)
!!$                    DEALLOCATE(new_symm_aux)
!!$                    ALLOCATE(new_symm_aux(tope+1+ll))
!!$                    new_symm_aux(1:tope)=cajon(:)
!!$                    tope=tope+1
!!$                    new_symm_aux(tope)=ll
!!$                    new_symm_aux((tope+1):(tope+ll))=order_aux(1:ll)
!!$                    tope=tope+ll
!!$                    DEALLOCATE(cajon)
!!$                 END IF
!!$              END IF
!!$           END IF
!!$        END DO
!!$        IF (interruptor.eqv..True.) THEN
!!$           interruptor=.FALSE.
!!$           new_num_crits=new_num_crits+1
!!$           IF (new_num_crits==1) THEN
!!$              tope=ll+1
!!$              ALLOCATE(new_symm_aux(tope))
!!$              new_symm_aux(1)=ll
!!$              new_symm_aux(2:tope)=order_aux(1:ll)
!!$           ELSE
!!$              ALLOCATE(cajon(tope))
!!$              cajon(:)=new_symm_aux(:)
!!$              DEALLOCATE(new_symm_aux)
!!$              ALLOCATE(new_symm_aux(tope+1+ll))
!!$              new_symm_aux(1:tope)=cajon(:)
!!$              tope=tope+1
!!$              new_symm_aux(tope)=ll
!!$              new_symm_aux((tope+1):(tope+ll))=order_aux(1:ll)
!!$              tope=tope+ll
!!$              DEALLOCATE(cajon)
!!$           END IF
!!$        END IF
!!$        DEALLOCATE(val_aux,ind_aux,order_aux)
!!$        DEALLOCATE(filtro,vals,inds)
!!$     END DO
!!$     
!!$     num_crits=new_num_crits
!!$     DEALLOCATE(vecti_aux)
!!$     IF (num_crits==0) THEN
!!$        interruptor=.FALSE.
!!$        dim_vecti_aux=0
!!$     ELSE
!!$        interruptor=.TRUE.
!!$        ALLOCATE(vecti_aux(tope))
!!$        dim_vecti_aux=tope
!!$        vecti_aux(:)=new_symm_aux(:)
!!$        DEALLOCATE(new_symm_aux)
!!$     END IF
!!$
!!$     pp=pp+1
!!$  END DO
!!$
!!$  DEALLOCATE(suprafiltro)
!!$
!!$
!!$END SUBROUTINE SORT_INT_MATRIX_ATS_1SH


!!!#################################
!!!####  BUILD FIRST SHELL
!!!#################################


!!$SUBROUTINE build_shell1st (core)
!!$ 
!!$  IMPLICIT NONE
!!$ 
!!$  INTEGER,INTENT(IN)::core
!!$ 
!!$  INTEGER::nn,nn2,mm,ii,jj,gg,kk,ll
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::order_ats_1sh,symm_ats_1sh,num_ats_bonded_1sh
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::order_ats_bonded_1sh,symm_ats_bonded_1sh
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::order_nodes_1sh,order_nodes_bonded_1sh
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::boxi,boxi2
!!$
!!$  IF (ALLOCATED(mss_ind_atoms)) DEALLOCATE(mss_ind_atoms)
!!$  IF (ALLOCATED(mss_ind_nodes)) DEALLOCATE(mss_ind_nodes)
!!$  IF (ALLOCATED(mss_symm))      DEALLOCATE(mss_symm)
!!$  IF (ALLOCATED(mss))           DEALLOCATE(mss)
!!$
!!$  nn=atomspernode(core)
!!$  nn2=2*nn
!!$ 
!!$  ALLOCATE(order_ats_1sh(nn),order_nodes_1sh(nn),symm_ats_1sh(nn),num_ats_bonded_1sh(nn2))
!!$
!!$  CALL order_ats_1st(core,nn,order_ats_1sh,symm_ats_1sh)
!!$
!!$  gg=0
!!$  DO ii=1,nn
!!$     jj=order_ats_1sh(ii)
!!$     num_ats_bonded_1sh(gg+1)=atom_hbs_num(jj)
!!$     num_ats_bonded_1sh(gg+2)=atom_bs_num(jj)
!!$     gg=gg+2
!!$  END DO
!!$
!!$  mm=SUM(num_ats_bonded_1sh(:),DIM=1)
!!$  
!!$  ALLOCATE(order_ats_bonded_1sh(mm),order_nodes_bonded_1sh(mm),symm_ats_bonded_1sh(mm))
!!$
!!$  gg=0
!!$  mm=0
!!$  symm_ats_bonded_1sh(:)=0
!!$  DO ii=1,nn
!!$     jj=order_ats_1sh(ii)
!!$     kk=num_ats_bonded_1sh(gg+1)
!!$     IF (kk>0) THEN
!!$        ll=T_hbs_start(jj)+1
!!$        IF (kk==1) THEN
!!$           mm=mm+1
!!$           order_ats_bonded_1sh(mm)=T_hbs_ind(ll)
!!$        ELSE
!!$           ALLOCATE(boxi(kk),boxi2(kk))
!!$           boxi(:)=T_hbs_ind((ll):(ll+kk-1))
!!$           CALL order_ats_bonded_1st(core,kk,boxi,boxi2)
!!$           order_ats_bonded_1sh((mm+1):(mm+kk))=boxi(:)
!!$           symm_ats_bonded_1sh((mm+1):(mm+kk))=boxi2(:)
!!$           DEALLOCATE(boxi,boxi2)
!!$           mm=mm+kk
!!$        END IF
!!$     END IF
!!$     kk=num_ats_bonded_1sh(gg+2)
!!$     IF (kk>0) THEN
!!$        ll=T_bs_start(jj)+1
!!$        IF (kk==1) THEN
!!$           mm=mm+1
!!$           order_ats_bonded_1sh(mm)=T_bs_ind(ll)
!!$        ELSE
!!$           ALLOCATE(boxi(kk),boxi2(kk))
!!$           boxi(:)=T_bs_ind((ll):(ll+kk-1))
!!$           CALL order_ats_bonded_1st(core,kk,boxi,boxi2)
!!$           order_ats_bonded_1sh((mm+1):(mm+kk))=boxi(:)
!!$           symm_ats_bonded_1sh((mm+1):(mm+kk))=boxi2(:)
!!$           DEALLOCATE(boxi,boxi2)
!!$           mm=mm+kk
!!$        END IF
!!$     END IF
!!$     gg=gg+2
!!$  END DO
!!$
!!$  ALLOCATE(mss_ind_atoms(1+nn+nn2+mm),mss_symm(1+nn+nn2+mm))
!!$  ALLOCATE(mss_ind_nodes(1+nn+nn2+mm),mss(1+nn+nn2+mm))
!!$
!!$  !Translate
!!$  order_nodes_1sh(:)        =atom2node(order_ats_1sh(:))
!!$  order_nodes_bonded_1sh(:) =atom2node(order_ats_bonded_1sh(:))
!!$  order_ats_1sh(:)          =trad2py_atom(order_ats_1sh(:))
!!$  order_ats_bonded_1sh(:)   =trad2py_atom(order_ats_bonded_1sh(:))
!!$  order_nodes_1sh(:)        =trad2py_node(order_nodes_1sh(:))
!!$  order_nodes_bonded_1sh(:) =trad2py_node(order_nodes_bonded_1sh(:))
!!$
!!$
!!$  !Build
!!$
!!$  mss_symm(:)=0
!!$  mss_ind_atoms(1)=nn
!!$  mss_ind_nodes(1)=nn
!!$  ii=2
!!$  jj=1+nn
!!$  mss_ind_atoms(ii:jj)=order_ats_1sh(:)
!!$  mss_ind_nodes(ii:jj)=order_nodes_1sh(:)
!!$  mss_symm(ii:jj)=symm_ats_1sh(:)
!!$  ii=jj+1
!!$  jj=jj+nn2
!!$  mss_ind_atoms(ii:jj)=num_ats_bonded_1sh(:)
!!$  mss_ind_nodes(ii:jj)=num_ats_bonded_1sh(:)
!!$  ii=jj+1
!!$  jj=jj+mm
!!$  mss_ind_atoms(ii:jj)=order_ats_bonded_1sh(:)
!!$  mss_ind_nodes(ii:jj)=order_nodes_bonded_1sh(:)
!!$  mss_symm(ii:jj)=symm_ats_bonded_1sh(:)
!!$
!!$  DEALLOCATE(order_ats_1sh,order_ats_bonded_1sh)
!!$  DEALLOCATE(order_nodes_1sh,order_nodes_bonded_1sh)
!!$  DEALLOCATE(symm_ats_1sh,symm_ats_bonded_1sh,num_ats_bonded_1sh)
!!$
!!$
!!$
!!$END SUBROUTINE build_shell1st
!!$
!!$
!!$SUBROUTINE order_ats_1st (core,num_ats,order,symm_ats_1sh)
!!$
!!$  IMPLICIT NONE
!!$
!!$  INTEGER,INTENT(IN)::core,num_ats
!!$  INTEGER,DIMENSION(num_ats),INTENT(OUT)::order,symm_ats_1sh
!!$
!!$  INTEGER::ii,jj,kk,gg,ll,num_crits
!!$  LOGICAL::interruptor
!!$
!!$  symm_ats_1sh(:)=0
!!$
!!$  jj=0
!!$  DO ii=node_run_ats(core)+1,node_run_ats(core+1)
!!$     jj=jj+1
!!$     order(jj)=ii
!!$  END DO
!!$ 
!!$  num_crits=symm_ats_crits(core)
!!$
!!$
!!$  IF (num_crits>0) THEN
!!$
!!$     ii=symm_ats_start(core)
!!$     jj=symm_ats_start(core+1)
!!$     gg=jj-ii
!!$     ALLOCATE(vecti_aux(gg))
!!$     vecti_aux(:)=symm_ats((ii+1):jj)
!!$
!!$     interruptor=.true.
!!$     IF (interruptor.eqv..true.) CALL SORTBYNUMHBS_ATS_1SH(core,num_ats,order,interruptor,num_crits,gg)
!!$     IF (interruptor.eqv..true.) CALL SORTBYNUMBS_ATS_1SH(core,num_ats,order,interruptor,num_crits,gg)
!!$
!!$     IF (interruptor.eqv..true.) THEN
!!$        gg=0
!!$        DO ii=1,num_crits
!!$           gg=gg+1
!!$           kk=vecti_aux(gg)
!!$           DO jj=1,kk
!!$              gg=gg+1
!!$              ll=vecti_aux(gg)
!!$              symm_ats_1sh(ll)=ii
!!$           END DO
!!$        END DO
!!$        DEALLOCATE(vecti_aux)
!!$        !print*,'Todavia',core
!!$     END IF
!!$
!!$  END IF
!!$
!!$END SUBROUTINE order_ats_1st
!!$
!!$SUBROUTINE order_ats_bonded_1st (core,num_ats,order,symm_ats_bonded_1sh)
!!$ 
!!$  IMPLICIT NONE
!!$ 
!!$  INTEGER,INTENT(IN)::core,num_ats
!!$  INTEGER,DIMENSION(num_ats),INTENT(INOUT)::order
!!$  INTEGER,DIMENSION(num_ats),INTENT(OUT)::symm_ats_bonded_1sh
!!$ 
!!$  INTEGER::ii,jj,kk,gg,ll,num_crits
!!$  LOGICAL::interruptor
!!$ 
!!$  symm_ats_bonded_1sh(:)=0
!!$
!!$  num_crits=1
!!$  gg=num_ats+1
!!$
!!$  ALLOCATE(vecti_aux(gg))
!!$  vecti_aux=(/num_ats,(ii,ii=1,num_ats,1)/)
!!$
!!$  interruptor=.true.
!!$
!!$  IF (interruptor.eqv..true.) CALL SORTBYCATEGORY_NODES_BONDED_1SH (core,num_ats,order,interruptor,num_crits,gg)
!!$  IF (interruptor.eqv..true.) CALL SORTBYSUPRASETS_NODES_BONDED_1SH (core,num_ats,order,interruptor,num_crits,gg)
!!$
!!$
!!$
!!$  IF (interruptor.eqv..true.) THEN
!!$     gg=0
!!$     DO ii=1,num_crits
!!$        gg=gg+1
!!$        kk=vecti_aux(gg)
!!$        DO jj=1,kk
!!$           gg=gg+1
!!$           ll=vecti_aux(gg)
!!$           symm_ats_bonded_1sh(ll)=ii
!!$        END DO
!!$     END DO
!!$     DEALLOCATE(vecti_aux)
!!$  END IF
!!$
!!$
!!$END SUBROUTINE order_ats_bonded_1st
!!$
!!$
!!$SUBROUTINE SORTBYNUMHBS_ATS_1SH (core,num_ats,order,interruptor,num_crits,dim_vecti_aux)
!!$
!!$  IMPLICIT NONE
!!$
!!$  INTEGER,INTENT(IN)::core,num_ats
!!$  INTEGER,INTENT(INOUT)::dim_vecti_aux
!!$  INTEGER,DIMENSION(num_ats),INTENT(INOUT)::order
!!$  LOGICAL,INTENT(INOUT)::interruptor
!!$  INTEGER,INTENT(INOUT)::num_crits
!!$
!!$  INTEGER::ii,jj,gg,kk,ll,idim
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::valores
!!$
!!$  ALLOCATE(valores(dim_vecti_aux))
!!$
!!$  gg=0
!!$  DO ii=1,num_crits
!!$     gg=gg+1
!!$     idim=vecti_aux(gg)
!!$     DO jj=1,idim
!!$        gg=gg+1
!!$        kk=vecti_aux(gg)
!!$        ll=order(kk)
!!$        valores(gg)=atom_hbs_num(ll)
!!$     END DO
!!$  END DO
!!$
!!$  CALL SORT_INT_ATS_1SH(num_ats,dim_vecti_aux,order,valores,interruptor,num_crits)
!!$
!!$  DEALLOCATE(valores)
!!$
!!$END SUBROUTINE SORTBYNUMHBS_ATS_1SH
!!$
!!$SUBROUTINE SORTBYNUMBS_ATS_1SH (core,num_ats,order,interruptor,num_crits,dim_vecti_aux)
!!$
!!$  IMPLICIT NONE
!!$
!!$  INTEGER,INTENT(IN)::core,num_ats
!!$  INTEGER,INTENT(INOUT)::dim_vecti_aux
!!$  INTEGER,DIMENSION(num_ats),INTENT(INOUT)::order
!!$  LOGICAL,INTENT(INOUT)::interruptor
!!$  INTEGER,INTENT(INOUT)::num_crits
!!$
!!$  INTEGER::ii,jj,gg,kk,ll,idim
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::valores
!!$
!!$  ALLOCATE(valores(dim_vecti_aux))
!!$
!!$  gg=0
!!$  DO ii=1,num_crits
!!$     gg=gg+1
!!$     idim=vecti_aux(gg)
!!$     DO jj=1,idim
!!$        gg=gg+1
!!$        kk=vecti_aux(gg)
!!$        ll=order(kk)
!!$        valores(gg)=atom_bs_num(ll)
!!$     END DO
!!$  END DO
!!$
!!$  CALL SORT_INT_ATS_1SH(num_ats,dim_vecti_aux,order,valores,interruptor,num_crits)
!!$
!!$  DEALLOCATE(valores)
!!$
!!$END SUBROUTINE SORTBYNUMBS_ATS_1SH
!!$
!!$SUBROUTINE SORTBYCATEGORY_NODES_BONDED_1SH (core,num_ats,order,interruptor,num_crits,dim_vecti_aux)
!!$
!!$  IMPLICIT NONE
!!$
!!$  INTEGER,INTENT(IN)::core,num_ats
!!$  INTEGER,INTENT(INOUT)::dim_vecti_aux
!!$  INTEGER,DIMENSION(num_ats),INTENT(INOUT)::order
!!$  LOGICAL,INTENT(INOUT)::interruptor
!!$  INTEGER,INTENT(INOUT)::num_crits
!!$
!!$  INTEGER::ii,jj,gg,kk,ll,idim
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::valores
!!$
!!$  ALLOCATE(valores(dim_vecti_aux))
!!$
!!$  gg=0
!!$  DO ii=1,num_crits
!!$     gg=gg+1
!!$     idim=vecti_aux(gg)
!!$     DO jj=1,idim
!!$        gg=gg+1
!!$        kk=vecti_aux(gg)
!!$        ll=order(kk)
!!$        valores(gg)=node_category(atom2node(ll))
!!$     END DO
!!$  END DO
!!$
!!$  CALL SORT_INT_ATS_1SH(num_ats,dim_vecti_aux,order,valores,interruptor,num_crits)
!!$
!!$  DEALLOCATE(valores)
!!$
!!$END SUBROUTINE SORTBYCATEGORY_NODES_BONDED_1SH
!!$
!!$
!!$SUBROUTINE SORTBYSUPRASETS_NODES_BONDED_1SH (core,num_ats,order,interruptor,num_crits,dim_vecti_aux)
!!$
!!$  IMPLICIT NONE
!!$
!!$  INTEGER,INTENT(IN)::core,num_ats
!!$  INTEGER,INTENT(INOUT)::dim_vecti_aux
!!$  INTEGER,DIMENSION(num_ats),INTENT(INOUT)::order
!!$  LOGICAL,INTENT(INOUT)::interruptor
!!$  INTEGER,INTENT(INOUT)::num_crits
!!$
!!$  INTEGER::ii,jj,gg,kk,ll,idim
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::valores
!!$
!!$  ALLOCATE(valores(dim_vecti_aux))
!!$
!!$  gg=0
!!$  DO ii=1,num_crits
!!$     gg=gg+1
!!$     idim=vecti_aux(gg)
!!$     DO jj=1,idim
!!$        gg=gg+1
!!$        kk=vecti_aux(gg)
!!$        ll=order(kk)
!!$        valores(gg)=node_category(atom2node(ll))
!!$     END DO
!!$  END DO
!!$
!!$  CALL SORT_INT_ATS_1SH(num_ats,dim_vecti_aux,order,valores,interruptor,num_crits)
!!$
!!$  DEALLOCATE(valores)
!!$
!!$END SUBROUTINE SORTBYSUPRASETS_NODES_BONDED_1SH
!!$


!!!#################################
!!!####  BUILD SECOND SHELL
!!!#################################




!!$SUBROUTINE build_order_ats_shell1st_2nd (core,core1st,num_ats,order,symm_ats_1sh)
!!$
!!$  IMPLICIT NONE
!!$
!!$  INTEGER,INTENT(IN)::core,core1st,num_ats
!!$  INTEGER,DIMENSION(num_ats),INTENT(OUT)::order,symm_ats_1sh
!!$
!!$  INTEGER::ii,jj,kk,gg,ll,num_crits
!!$  LOGICAL::interruptor
!!$
!!$  symm_ats_1sh(:)=0
!!$
!!$  jj=0
!!$  DO ii=node_run_ats(core)+1,node_run_ats(core+1)
!!$     jj=jj+1
!!$     order(jj)=ii
!!$  END DO
!!$ 
!!$  num_crits=symm_ats_crits(core)
!!$
!!$
!!$  IF (num_crits>0) THEN
!!$
!!$     ii=symm_ats_start(core)
!!$     jj=symm_ats_start(core+1)
!!$     gg=jj-ii
!!$     ALLOCATE(vecti_aux(gg))
!!$     vecti_aux(:)=symm_ats((ii+1):jj)
!!$
!!$     interruptor=.true.
!!$     IF (interruptor.eqv..true.) CALL SORTBYCORE_ATS_1SH (core,core1st,num_ats,order,interruptor,num_crits,gg)
!!$     IF (interruptor.eqv..true.) CALL SORTBYNUMHBS_ATS_1SH(core,num_ats,order,interruptor,num_crits,gg)
!!$     IF (interruptor.eqv..true.) CALL SORTBYNUMBS_ATS_1SH(core,num_ats,order,interruptor,num_crits,gg)
!!$     IF (interruptor.eqv..true.) CALL SORTBYCATEGORYHBS_ATS_1SH(core,num_ats,order,interruptor,num_crits,gg)
!!$     IF (interruptor.eqv..true.) CALL SORTBYCATEGORYBS_ATS_1SH(core,num_ats,order,interruptor,num_crits,gg)
!!$
!!$     IF (interruptor.eqv..true.) THEN
!!$        gg=0
!!$        DO ii=1,num_crits
!!$           gg=gg+1
!!$           kk=vecti_aux(gg)
!!$           DO jj=1,kk
!!$              gg=gg+1
!!$              ll=vecti_aux(gg)
!!$              symm_ats_1sh(ll)=ii
!!$           END DO
!!$        END DO
!!$        DEALLOCATE(vecti_aux)
!!$     END IF
!!$
!!$  END IF
!!$
!!$END SUBROUTINE build_order_ats_shell1st_2nd
!!$
!!$
!!$!!#### SORTING 1ST SHELL:
!!$
!!$
!!$SUBROUTINE SORTBYCORE_ATS_1SH(core,core1st,num_ats,order,interruptor,num_crits,dim_vecti_aux)
!!$
!!$  IMPLICIT NONE
!!$
!!$  INTEGER,INTENT(IN)::core,core1st,num_ats
!!$  INTEGER,INTENT(INOUT)::dim_vecti_aux
!!$  INTEGER,DIMENSION(num_ats),INTENT(INOUT)::order
!!$  LOGICAL,INTENT(INOUT)::interruptor
!!$  INTEGER,INTENT(INOUT)::num_crits
!!$
!!$  INTEGER::ii,jj,gg,kk,ll,mm,nn,pp,idim
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::valores
!!$     
!!$  ALLOCATE(valores(dim_vecti_aux))
!!$  valores(:)=0
!!$     
!!$  gg=0
!!$  DO ii=1,num_crits
!!$     gg=gg+1
!!$     idim=vecti_aux(gg)
!!$     DO jj=1,idim
!!$        gg=gg+1
!!$        kk=vecti_aux(gg)
!!$        ll=order(kk)
!!$        DO mm=T_hbs_start(ll)+1,T_hbs_start(ll+1)
!!$           nn=T_hbs_ind(mm)
!!$           IF (atom2node(nn)==core1st) valores(gg)=1
!!$        END DO
!!$        DO mm=T_bs_start(ll)+1,T_bs_start(ll+1)
!!$           nn=T_bs_ind(mm)
!!$           IF (atom2node(nn)==core1st) valores(gg)=1
!!$        END DO
!!$     END DO
!!$  END DO
!!$  
!!$  CALL SORT_INT_ATS_1SH(num_ats,dim_vecti_aux,order,valores,interruptor,num_crits)
!!$     
!!$  DEALLOCATE(valores)
!!$
!!$END SUBROUTINE SORTBYCORE_ATS_1SH
!!$
!!$
!!$SUBROUTINE SORTBYCATEGORYHBS_ATS_1SH(core,num_ats,order,interruptor,num_crits,dim_vecti_aux)
!!$
!!$  IMPLICIT NONE
!!$
!!$  INTEGER,INTENT(IN)::core,num_ats
!!$  INTEGER,INTENT(INOUT)::dim_vecti_aux
!!$  INTEGER,DIMENSION(num_ats),INTENT(INOUT)::order
!!$  LOGICAL,INTENT(INOUT)::interruptor
!!$  INTEGER,INTENT(INOUT)::num_crits
!!$
!!$  INTEGER::ii,jj,gg,kk,ll,mm,nn,pp,idim,dim_matrix
!!$  INTEGER,DIMENSION(:,:),ALLOCATABLE::valores
!!$
!!$  dim_matrix=0
!!$  gg=0
!!$  DO ii=1,num_crits
!!$     gg=gg+1
!!$     idim=vecti_aux(gg)
!!$     DO jj=1,idim
!!$        gg=gg+1
!!$        kk=vecti_aux(gg)
!!$        ll=order(kk)
!!$        mm=atom_hbs_num(ll)
!!$        IF (dim_matrix<mm) dim_matrix=mm
!!$     END DO
!!$  END DO
!!$
!!$  IF (dim_matrix>0) THEN
!!$     
!!$     ALLOCATE(valores(dim_vecti_aux,dim_matrix))
!!$     valores(:,:)=0
!!$     
!!$     gg=0
!!$     DO ii=1,num_crits
!!$        gg=gg+1
!!$        idim=vecti_aux(gg)
!!$        DO jj=1,idim
!!$           gg=gg+1
!!$           kk=vecti_aux(gg)
!!$           ll=order(kk)
!!$           pp=0
!!$           DO mm=T_hbs_start(ll)+1,T_hbs_start(ll+1)
!!$              pp=pp+1
!!$              nn=T_hbs_ind(mm)
!!$              nn=atom2node(nn)
!!$              nn=node_category(nn)
!!$              valores(gg,pp)=nn
!!$           END DO
!!$        END DO
!!$     END DO
!!$     
!!$     CALL SORT_INT_MATRIX_ATS_1SH(num_ats,dim_vecti_aux,dim_matrix,order,valores,interruptor,num_crits)
!!$     
!!$     DEALLOCATE(valores)
!!$
!!$  END IF
!!$
!!$END SUBROUTINE SORTBYCATEGORYHBS_ATS_1SH
!!$
!!$SUBROUTINE SORTBYCATEGORYBS_ATS_1SH(core,num_ats,order,interruptor,num_crits,dim_vecti_aux)
!!$
!!$  IMPLICIT NONE
!!$
!!$  INTEGER,INTENT(IN)::core,num_ats
!!$  INTEGER,INTENT(INOUT)::dim_vecti_aux
!!$  INTEGER,DIMENSION(num_ats),INTENT(INOUT)::order
!!$  LOGICAL,INTENT(INOUT)::interruptor
!!$  INTEGER,INTENT(INOUT)::num_crits
!!$
!!$  INTEGER::ii,jj,gg,kk,ll,mm,nn,pp,idim,dim_matrix
!!$  INTEGER,DIMENSION(:,:),ALLOCATABLE::valores
!!$
!!$  dim_matrix=0
!!$  gg=0
!!$  DO ii=1,num_crits
!!$     gg=gg+1
!!$     idim=vecti_aux(gg)
!!$     DO jj=1,idim
!!$        gg=gg+1
!!$        kk=vecti_aux(gg)
!!$        ll=order(kk)
!!$        mm=atom_hbs_num(ll)
!!$        IF (dim_matrix<mm) dim_matrix=mm
!!$     END DO
!!$  END DO
!!$
!!$  IF (dim_matrix>0) THEN
!!$     
!!$     ALLOCATE(valores(dim_vecti_aux,dim_matrix))
!!$     valores(:,:)=0
!!$     
!!$     gg=0
!!$     DO ii=1,num_crits
!!$        gg=gg+1
!!$        idim=vecti_aux(gg)
!!$        DO jj=1,idim
!!$           gg=gg+1
!!$           kk=vecti_aux(gg)
!!$           ll=order(kk)
!!$           pp=0
!!$           DO mm=T_bs_start(ll)+1,T_bs_start(ll+1)
!!$              pp=pp+1
!!$              nn=T_bs_ind(mm)
!!$              nn=atom2node(nn)
!!$              nn=node_category(nn)
!!$              valores(gg,pp)=nn
!!$           END DO
!!$        END DO
!!$     END DO
!!$     
!!$     CALL SORT_INT_MATRIX_ATS_1SH(num_ats,dim_vecti_aux,dim_matrix,order,valores,interruptor,num_crits)
!!$     
!!$     DEALLOCATE(valores)
!!$
!!$  END IF
!!$
!!$END SUBROUTINE SORTBYCATEGORYBS_ATS_1SH
!!$
!!$SUBROUTINE SORTBYNUMHBS_NODES_1SH (core,num_ats,order,interruptor,num_crits,dim_vecti_aux)
!!$
!!$  IMPLICIT NONE
!!$
!!$  INTEGER,INTENT(IN)::core,num_ats
!!$  INTEGER,INTENT(INOUT)::dim_vecti_aux
!!$  INTEGER,DIMENSION(num_ats),INTENT(INOUT)::order
!!$  LOGICAL,INTENT(INOUT)::interruptor
!!$  INTEGER,INTENT(INOUT)::num_crits
!!$
!!$  INTEGER::ii,jj,gg,kk,ll,idim
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::valores
!!$
!!$  ALLOCATE(valores(dim_vecti_aux))
!!$
!!$  gg=0
!!$  DO ii=1,num_crits
!!$     gg=gg+1
!!$     idim=vecti_aux(gg)
!!$     DO jj=1,idim
!!$        gg=gg+1
!!$        kk=vecti_aux(gg)
!!$        ll=order(kk)
!!$        valores(gg)=node_hbs_num(atom2node(ll))
!!$     END DO
!!$  END DO
!!$
!!$  CALL SORT_INT_ATS_1SH(num_ats,dim_vecti_aux,order,valores,interruptor,num_crits)
!!$
!!$  DEALLOCATE(valores)
!!$
!!$END SUBROUTINE SORTBYNUMHBS_NODES_1SH
!!$
!!$
!!$SUBROUTINE SORTBYNUMBS_NODES_1SH (core,num_ats,order,interruptor,num_crits,dim_vecti_aux)
!!$
!!$  IMPLICIT NONE
!!$
!!$  INTEGER,INTENT(IN)::core,num_ats
!!$  INTEGER,INTENT(INOUT)::dim_vecti_aux
!!$  INTEGER,DIMENSION(num_ats),INTENT(INOUT)::order
!!$  LOGICAL,INTENT(INOUT)::interruptor
!!$  INTEGER,INTENT(INOUT)::num_crits
!!$
!!$  INTEGER::ii,jj,gg,kk,ll,idim
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::valores
!!$
!!$  ALLOCATE(valores(dim_vecti_aux))
!!$
!!$  gg=0
!!$  DO ii=1,num_crits
!!$     gg=gg+1
!!$     idim=vecti_aux(gg)
!!$     DO jj=1,idim
!!$        gg=gg+1
!!$        kk=vecti_aux(gg)
!!$        ll=order(kk)
!!$        valores(gg)=node_bs_num(atom2node(ll))
!!$     END DO
!!$  END DO
!!$
!!$  CALL SORT_INT_ATS_1SH(num_ats,dim_vecti_aux,order,valores,interruptor,num_crits)
!!$
!!$  DEALLOCATE(valores)
!!$
!!$END SUBROUTINE SORTBYNUMBS_NODES_1SH
!!$
!!$
!!$
!!$SUBROUTINE SORTBYCORE_NODES_BONDED_1SH (core,core1st,num_ats,order,interruptor,num_crits,dim_vecti_aux)
!!$
!!$  IMPLICIT NONE
!!$
!!$  INTEGER,INTENT(IN)::core,core1st,num_ats
!!$  INTEGER,INTENT(INOUT)::dim_vecti_aux
!!$  INTEGER,DIMENSION(num_ats),INTENT(INOUT)::order
!!$  LOGICAL,INTENT(INOUT)::interruptor
!!$  INTEGER,INTENT(INOUT)::num_crits
!!$
!!$  INTEGER::ii,jj,gg,kk,ll,idim
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::valores
!!$
!!$  ALLOCATE(valores(dim_vecti_aux))
!!$  valores(:)=0
!!$
!!$  gg=0
!!$  DO ii=1,num_crits
!!$     gg=gg+1
!!$     idim=vecti_aux(gg)
!!$     DO jj=1,idim
!!$        gg=gg+1
!!$        kk=vecti_aux(gg)
!!$        ll=order(kk)
!!$        IF (atom2node(ll)==core1st) valores(gg)=1
!!$     END DO
!!$  END DO
!!$
!!$  CALL SORT_INT_ATS_1SH(num_ats,dim_vecti_aux,order,valores,interruptor,num_crits)
!!$
!!$  DEALLOCATE(valores)
!!$
!!$END SUBROUTINE SORTBYCORE_NODES_BONDED_1SH
!!$
!!$
!!$
!!$
!!$
!!$SUBROUTINE build_order_ats_bonded_shell1st_2nd (core,core1st,num_ats,order,symm_ats_bonded_1sh)
!!$ 
!!$  IMPLICIT NONE
!!$ 
!!$  INTEGER,INTENT(IN)::core,core1st,num_ats
!!$  INTEGER,DIMENSION(num_ats),INTENT(INOUT)::order
!!$  INTEGER,DIMENSION(num_ats),INTENT(OUT)::symm_ats_bonded_1sh
!!$ 
!!$  INTEGER::ii,jj,kk,gg,ll,num_crits
!!$  LOGICAL::interruptor
!!$ 
!!$  symm_ats_bonded_1sh(:)=0
!!$
!!$  num_crits=1
!!$  gg=num_ats+1
!!$
!!$  ALLOCATE(vecti_aux(gg))
!!$  vecti_aux=(/num_ats,(ii,ii=1,num_ats,1)/)
!!$ 
!!$
!!$  interruptor=.true.
!!$
!!$  IF (interruptor.eqv..true.) CALL SORTBYCORE_NODES_BONDED_1SH (core,core1st,num_ats,order,interruptor,num_crits,gg)
!!$  IF (interruptor.eqv..true.) CALL SORTBYCATEGORY_NODES_BONDED_1SH (core,num_ats,order,interruptor,num_crits,gg)
!!$
!!$  !IF (interruptor.eqv..true.) CALL SORTBYNUMHBS_NODES_1SH(core,num_ats,order,interruptor,num_crits,gg)
!!$  !IF (interruptor.eqv..true.) CALL SORTBYNUMBS_NODES_1SH(core,num_ats,order,interruptor,num_crits,gg)
!!$
!!$  IF (interruptor.eqv..true.) THEN
!!$     gg=0
!!$     DO ii=1,num_crits
!!$        gg=gg+1
!!$        kk=vecti_aux(gg)
!!$        DO jj=1,kk
!!$           gg=gg+1
!!$           ll=vecti_aux(gg)
!!$           symm_ats_bonded_1sh(ll)=ii
!!$        END DO
!!$     END DO
!!$     DEALLOCATE(vecti_aux)
!!$  END IF
!!$
!!$
!!$END SUBROUTINE build_order_ats_bonded_shell1st_2nd
!!$
!!$
!!$!#######################################
!!$!#######################################
!!$!#######################################
!!$!#######################################
!!$!#######################################
!!$!#######################################
!!$
!!$SUBROUTINE build_shell2nd (core)
!!$
!!$  IMPLICIT NONE
!!$ 
!!$  TYPE iarray_pointer
!!$     INTEGER,DIMENSION(:),POINTER::p1
!!$  END TYPE iarray_pointer
!!$
!!$  INTEGER,INTENT(IN)::core
!!$ 
!!$  INTEGER::nn,nn2,mm,ii,jj,gg,kk,ll
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::order_ats_1sh,symm_ats_1sh,num_ats_bonded_1sh
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::order_ats_bonded_1sh,symm_ats_bonded_1sh
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::order_nodes_1sh,order_nodes_bonded_1sh
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::boxi,boxi2
!!$
!!$  INTEGER::core2
!!$  INTEGER::nnn,nnn2,mmm,iii
!!$  !INTEGER,DIMENSION(:),ALLOCATABLE::order_ats_2nd,symm_ats_2nd,num_ats_bonded_2nd
!!$  !INTEGER,DIMENSION(:),ALLOCATABLE::order_ats_bonded_2nd,symm_ats_bonded_2nd
!!$  !INTEGER,DIMENSION(:),ALLOCATABLE::order_nodes_2nd,order_nodes_bonded_2nd
!!$  TYPE(iarray_pointer),DIMENSION(:),POINTER::order_ats_2nd,symm_ats_2nd,num_ats_bonded_2nd
!!$  TYPE(iarray_pointer),DIMENSION(:),POINTER::order_ats_bonded_2nd,symm_ats_bonded_2nd
!!$  TYPE(iarray_pointer),DIMENSION(:),POINTER::order_nodes_2nd,order_nodes_bonded_2nd
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::enes,enes2,emes
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::aux_order_ats_2nd,aux_symm_ats_2nd,aux_num_ats_bonded_2nd
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::aux_order_ats_bonded_2nd,aux_symm_ats_bonded_2nd
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::aux_order_nodes_2nd,aux_order_nodes_bonded_2nd
!!$
!!$
!!$  IF (ALLOCATED(mss_ind_atoms)) DEALLOCATE(mss_ind_atoms)
!!$  IF (ALLOCATED(mss_ind_nodes)) DEALLOCATE(mss_ind_nodes)
!!$  IF (ALLOCATED(mss_symm))      DEALLOCATE(mss_symm)
!!$  IF (ALLOCATED(mss))           DEALLOCATE(mss)
!!$
!!$  !core
!!$
!!$  nn=atomspernode(core)
!!$  nn2=2*nn
!!$ 
!!$  ALLOCATE(order_ats_1sh(nn),order_nodes_1sh(nn),symm_ats_1sh(nn),num_ats_bonded_1sh(nn2))
!!$
!!$  CALL build_order_ats_shell1st(core,nn,order_ats_1sh,symm_ats_1sh)
!!$
!!$  gg=0
!!$  DO ii=1,nn
!!$     jj=order_ats_1sh(ii)
!!$     num_ats_bonded_1sh(gg+1)=atom_hbs_num(jj)
!!$     num_ats_bonded_1sh(gg+2)=atom_bs_num(jj)
!!$     gg=gg+2
!!$  END DO
!!$
!!$  mm=SUM(num_ats_bonded_1sh(:),DIM=1)
!!$  
!!$  ALLOCATE(order_ats_bonded_1sh(mm),order_nodes_bonded_1sh(mm),symm_ats_bonded_1sh(mm))
!!$
!!$  gg=0
!!$  mm=0
!!$  symm_ats_bonded_1sh(:)=0
!!$  DO ii=1,nn
!!$     jj=order_ats_1sh(ii)
!!$     kk=num_ats_bonded_1sh(gg+1)
!!$     IF (kk>0) THEN
!!$        ll=T_hbs_start(jj)+1
!!$        IF (kk==1) THEN
!!$           mm=mm+1
!!$           order_ats_bonded_1sh(mm)=T_hbs_ind(ll)
!!$        ELSE
!!$           ALLOCATE(boxi(kk),boxi2(kk))
!!$           boxi(:)=T_hbs_ind((ll):(ll+kk-1))
!!$           CALL build_order_ats_bonded_shell1st(core,kk,boxi,boxi2)
!!$           order_ats_bonded_1sh((mm+1):(mm+kk))=boxi(:)
!!$           symm_ats_bonded_1sh((mm+1):(mm+kk))=boxi2(:)
!!$           DEALLOCATE(boxi,boxi2)
!!$           mm=mm+kk
!!$        END IF
!!$     END IF
!!$     kk=num_ats_bonded_1sh(gg+2)
!!$     IF (kk>0) THEN
!!$        ll=T_bs_start(jj)+1
!!$        IF (kk==1) THEN
!!$           mm=mm+1
!!$           order_ats_bonded_1sh(mm)=T_bs_ind(ll)
!!$        ELSE
!!$           ALLOCATE(boxi(kk),boxi2(kk))
!!$           boxi(:)=T_bs_ind((ll):(ll+kk-1))
!!$           CALL build_order_ats_bonded_shell1st(core,kk,boxi,boxi2)
!!$           order_ats_bonded_1sh((mm+1):(mm+kk))=boxi(:)
!!$           symm_ats_bonded_1sh((mm+1):(mm+kk))=boxi2(:)
!!$           DEALLOCATE(boxi,boxi2)
!!$           mm=mm+kk
!!$        END IF
!!$     END IF
!!$     gg=gg+2
!!$  END DO
!!$
!!$  !!! for the second shell
!!$
!!$  order_nodes_1sh(:)        =atom2node(order_ats_1sh(:))
!!$  order_nodes_bonded_1sh(:) =atom2node(order_ats_bonded_1sh(:))
!!$
!!$
!!$  ALLOCATE(order_ats_2nd(mm),symm_ats_2nd(mm),num_ats_bonded_2nd(mm))
!!$  ALLOCATE(order_ats_bonded_2nd(mm),symm_ats_bonded_2nd(mm))
!!$  ALLOCATE(order_nodes_2nd(mm),order_nodes_bonded_2nd(mm))
!!$  ALLOCATE(enes(mm),enes2(mm),emes(mm))
!!$
!!$  DO iii=1,mm
!!$     core2=order_nodes_bonded_1sh(iii)
!!$
!!$     nnn=atomspernode(core2)
!!$     nnn2=2*nnn
!!$
!!$     ALLOCATE(aux_order_ats_2nd(nnn),aux_order_nodes_2nd(nnn),aux_symm_ats_2nd(nnn),aux_num_ats_bonded_2nd(nnn2))
!!$
!!$     CALL build_order_ats_shell1st_2nd(core2,core,nnn,aux_order_ats_2nd,aux_symm_ats_2nd)
!!$
!!$     gg=0
!!$     DO ii=1,nnn
!!$        jj=aux_order_ats_2nd(ii)
!!$        aux_num_ats_bonded_2nd(gg+1)=atom_hbs_num(jj)
!!$        aux_num_ats_bonded_2nd(gg+2)=atom_bs_num(jj)
!!$        gg=gg+2
!!$     END DO
!!$
!!$     mmm=SUM(aux_num_ats_bonded_2nd(:),DIM=1)
!!$  
!!$     ALLOCATE(aux_order_ats_bonded_2nd(mmm),aux_order_nodes_bonded_2nd(mmm),aux_symm_ats_bonded_2nd(mmm))
!!$
!!$     gg=0
!!$     mmm=0
!!$     aux_symm_ats_bonded_2nd(:)=0
!!$     DO ii=1,nnn
!!$        jj=aux_order_ats_2nd(ii)
!!$        kk=aux_num_ats_bonded_2nd(gg+1)
!!$        IF (kk>0) THEN
!!$           ll=T_hbs_start(jj)+1
!!$           IF (kk==1) THEN
!!$              mmm=mmm+1
!!$              aux_order_ats_bonded_2nd(mmm)=T_hbs_ind(ll)
!!$           ELSE
!!$              ALLOCATE(boxi(kk),boxi2(kk))
!!$              boxi(:)=T_hbs_ind((ll):(ll+kk-1))
!!$              CALL build_order_ats_bonded_shell1st_2nd(core2,core,kk,boxi,boxi2)
!!$              aux_order_ats_bonded_2nd((mmm+1):(mmm+kk))=boxi(:)
!!$              aux_symm_ats_bonded_2nd((mmm+1):(mmm+kk))=boxi2(:)
!!$              DEALLOCATE(boxi,boxi2)
!!$              mmm=mmm+kk
!!$           END IF
!!$        END IF
!!$        kk=aux_num_ats_bonded_2nd(gg+2)
!!$        IF (kk>0) THEN
!!$           ll=T_bs_start(jj)+1
!!$           IF (kk==1) THEN
!!$              mmm=mmm+1
!!$              aux_order_ats_bonded_2nd(mmm)=T_bs_ind(ll)
!!$           ELSE
!!$              ALLOCATE(boxi(kk),boxi2(kk))
!!$              boxi(:)=T_bs_ind((ll):(ll+kk-1))
!!$              CALL build_order_ats_bonded_shell1st_2nd(core2,core,kk,boxi,boxi2)
!!$              aux_order_ats_bonded_2nd((mmm+1):(mmm+kk))=boxi(:)
!!$              aux_symm_ats_bonded_2nd((mmm+1):(mmm+kk))=boxi2(:)
!!$              DEALLOCATE(boxi,boxi2)
!!$              mmm=mmm+kk
!!$           END IF
!!$        END IF
!!$        gg=gg+2
!!$     END DO
!!$     
!!$     aux_order_nodes_2nd(:)        =atom2node(aux_order_ats_2nd(:))
!!$     aux_order_nodes_bonded_2nd(:) =atom2node(aux_order_ats_bonded_2nd(:))
!!$     enes(iii)=nnn
!!$     enes2(iii)=nnn2
!!$     emes(iii)=mmm
!!$
!!$     ALLOCATE(order_ats_2nd(iii)%p1(nnn),order_nodes_2nd(iii)%p1(nnn),symm_ats_2nd(iii)%p1(nnn),num_ats_bonded_2nd(iii)%p1(nnn2))
!!$     ALLOCATE(order_ats_bonded_2nd(iii)%p1(mmm),order_nodes_bonded_2nd(iii)%p1(mmm),symm_ats_bonded_2nd(iii)%p1(mmm))
!!$
!!$     order_ats_2nd(iii)%p1(:)=aux_order_ats_2nd(:)
!!$     order_nodes_2nd(iii)%p1(:)=aux_order_nodes_2nd(:)
!!$     symm_ats_2nd(iii)%p1(:)=aux_symm_ats_2nd(:)
!!$     num_ats_bonded_2nd(iii)%p1(:)=aux_num_ats_bonded_2nd(:)
!!$     order_ats_bonded_2nd(iii)%p1(:)=aux_order_ats_bonded_2nd(:)
!!$     order_nodes_bonded_2nd(iii)%p1(:)=aux_order_nodes_bonded_2nd(:)
!!$     symm_ats_bonded_2nd(iii)%p1(:)=aux_symm_ats_bonded_2nd(:)
!!$     
!!$     DEALLOCATE(aux_order_ats_bonded_2nd,aux_order_nodes_bonded_2nd,aux_symm_ats_bonded_2nd)
!!$     DEALLOCATE(aux_order_ats_2nd,aux_order_nodes_2nd,aux_symm_ats_2nd,aux_num_ats_bonded_2nd)
!!$
!!$  END DO
!!$  !!!!
!!$
!!$  !Translate
!!$  order_ats_1sh(:)          =trad2py_atom(order_ats_1sh(:))
!!$  order_ats_bonded_1sh(:)   =trad2py_atom(order_ats_bonded_1sh(:))
!!$  order_nodes_1sh(:)        =trad2py_node(order_nodes_1sh(:))
!!$  order_nodes_bonded_1sh(:) =trad2py_node(order_nodes_bonded_1sh(:))
!!$
!!$  DO iii=1,mm
!!$     order_ats_2nd(iii)%p1(:)          =trad2py_atom(order_ats_2nd(iii)%p1(:))
!!$     order_ats_bonded_2nd(iii)%p1(:)   =trad2py_atom(order_ats_bonded_2nd(iii)%p1(:))
!!$     order_nodes_2nd(iii)%p1(:)        =trad2py_node(order_nodes_2nd(iii)%p1(:))
!!$     order_nodes_bonded_2nd(iii)%p1(:) =trad2py_node(order_nodes_bonded_2nd(iii)%p1(:))
!!$  END DO
!!$
!!$  !Build
!!$
!!$  nnn=SUM(enes(:),DIM=1)
!!$  nnn2=SUM(enes2(:),DIM=1)
!!$  mmm=SUM(emes(:),DIM=1)
!!$
!!$  ALLOCATE(mss_ind_atoms(1+nn+nn2+mm+mm+nnn+nnn2+mmm),mss_symm(1+nn+nn2+mm+mm+nnn+nnn2+mmm))
!!$  ALLOCATE(mss_ind_nodes(1+nn+nn2+mm+mm+nnn+nnn2+mmm),mss(1+nn+nn2+mm+mm+nnn+nnn2+mmm))
!!$
!!$
!!$  mss_symm(:)=0
!!$  mss_ind_atoms(1)=nn
!!$  mss_ind_nodes(1)=nn
!!$  ii=2
!!$  jj=1+nn
!!$  mss_ind_atoms(ii:jj)=order_ats_1sh(:)
!!$  mss_ind_nodes(ii:jj)=order_nodes_1sh(:)
!!$  mss_symm(ii:jj)=symm_ats_1sh(:)
!!$  ii=jj+1
!!$  jj=jj+nn2
!!$  mss_ind_atoms(ii:jj)=num_ats_bonded_1sh(:)
!!$  mss_ind_nodes(ii:jj)=num_ats_bonded_1sh(:)
!!$  ii=jj+1
!!$  jj=jj+mm
!!$  mss_ind_atoms(ii:jj)=order_ats_bonded_1sh(:)
!!$  mss_ind_nodes(ii:jj)=order_nodes_bonded_1sh(:)
!!$  mss_symm(ii:jj)=symm_ats_bonded_1sh(:)
!!$  DO iii=1,mm
!!$     ii=jj+1
!!$     jj=jj+1
!!$     mss_ind_atoms(ii)=enes(iii)
!!$     mss_ind_nodes(ii)=enes(iii)
!!$     ii=jj+1
!!$     jj=jj+enes(iii)
!!$     mss_ind_atoms(ii:jj)=order_ats_2nd(iii)%p1(:)
!!$     mss_ind_nodes(ii:jj)=order_nodes_2nd(iii)%p1(:)
!!$     mss_symm(ii:jj)=symm_ats_2nd(iii)%p1(:)
!!$     ii=jj+1
!!$     jj=jj+enes2(iii)
!!$     mss_ind_atoms(ii:jj)=num_ats_bonded_2nd(iii)%p1(:)
!!$     mss_ind_nodes(ii:jj)=num_ats_bonded_2nd(iii)%p1(:)
!!$     ii=jj+1
!!$     jj=jj+emes(iii)
!!$     mss_ind_atoms(ii:jj)=order_ats_bonded_2nd(iii)%p1(:)
!!$     mss_ind_nodes(ii:jj)=order_nodes_bonded_2nd(iii)%p1(:)
!!$     mss_symm(ii:jj)=symm_ats_bonded_2nd(iii)%p1(:)
!!$     DEALLOCATE(order_ats_2nd(iii)%p1,order_ats_bonded_2nd(iii)%p1)
!!$     DEALLOCATE(order_nodes_2nd(iii)%p1,order_nodes_bonded_2nd(iii)%p1)
!!$     DEALLOCATE(symm_ats_2nd(iii)%p1,symm_ats_bonded_2nd(iii)%p1,num_ats_bonded_2nd(iii)%p1)
!!$  END DO
!!$
!!$
!!$  DEALLOCATE(order_ats_1sh,order_ats_bonded_1sh)
!!$  DEALLOCATE(order_nodes_1sh,order_nodes_bonded_1sh)
!!$  DEALLOCATE(symm_ats_1sh,symm_ats_bonded_1sh,num_ats_bonded_1sh)
!!$  DEALLOCATE(order_ats_2nd,order_ats_bonded_2nd)
!!$  DEALLOCATE(order_nodes_2nd,order_nodes_bonded_2nd)
!!$  DEALLOCATE(symm_ats_2nd,symm_ats_bonded_2nd,num_ats_bonded_2nd)
!!$  DEALLOCATE(enes,enes2,emes)
!!$
!!$
!!$END SUBROUTINE build_shell2nd



END MODULE GLOB
