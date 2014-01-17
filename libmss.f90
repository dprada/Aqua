MODULE GLOB

INTEGER::num_nodes,num_atoms
INTEGER,DIMENSION(:),ALLOCATABLE::node_run_ats,atom2node,trad2py_node,trad2py_atom,atomspernode
INTEGER,DIMENSION(:),ALLOCATABLE::symm_ats_crits,symm_ats_start,symm_ats

INTEGER,DIMENSION(:),ALLOCATABLE::vecti_aux

INTEGER::T_num_hbs,T_num_bs
INTEGER,DIMENSION(:),ALLOCATABLE::T_hbs_start,T_bs_start,T_hbs_ind,T_bs_ind
INTEGER,DIMENSION(:),ALLOCATABLE::T_hbs_num,T_bs_num

INTEGER,DIMENSION(:),ALLOCATABLE::mss_ind_atoms,mss_ind_nodes,mss_symm,mss

CONTAINS

SUBROUTINE load_topol(inode_run_ats,iatom2node,itrad2py_node,itrad2py_atom,&
     isymm_ats_start,isymm_ats_crits,isymm_ats,inum_nodes,inum_atoms,isymm_ats_dim)

  IMPLICIT NONE
  INTEGER,INTENT(IN)::inum_nodes,inum_atoms,isymm_ats_dim
  INTEGER,DIMENSION(inum_nodes+1),INTENT(IN)::inode_run_ats
  INTEGER,DIMENSION(inum_atoms),INTENT(IN)::iatom2node
  INTEGER,DIMENSION(inum_nodes),INTENT(IN)::itrad2py_node
  INTEGER,DIMENSION(inum_atoms),INTENT(IN)::itrad2py_atom
  INTEGER,DIMENSION(inum_nodes+1),INTENT(IN)::isymm_ats_start
  INTEGER,DIMENSION(inum_nodes),INTENT(IN)::isymm_ats_crits
  INTEGER,DIMENSION(isymm_ats_dim),INTENT(IN)::isymm_ats

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

END SUBROUTINE load_topol

SUBROUTINE load_net(iT_hbs_start,iT_bs_start,iT_hbs_ind,iT_bs_ind,iT_num_hbs,iT_num_bs,inum_atoms)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::inum_atoms,iT_num_hbs,iT_num_bs
  INTEGER,DIMENSION(inum_atoms+1),INTENT(IN)::iT_hbs_start,iT_bs_start
  INTEGER,DIMENSION(iT_num_hbs),INTENT(IN)::iT_hbs_ind
  INTEGER,DIMENSION(iT_num_bs),INTENT(IN)::iT_bs_ind

  INTEGER::ii

  IF (ALLOCATED(T_hbs_start))  DEALLOCATE(T_hbs_start)
  IF (ALLOCATED(T_bs_start))   DEALLOCATE(T_bs_start)
  IF (ALLOCATED(T_hbs_ind))    DEALLOCATE(T_hbs_ind)
  IF (ALLOCATED(T_bs_ind))     DEALLOCATE(T_bs_ind)
  IF (ALLOCATED(T_hbs_num))    DEALLOCATE(T_hbs_num)
  IF (ALLOCATED(T_bs_num))     DEALLOCATE(T_bs_num)

  T_num_hbs=iT_num_hbs
  T_num_bs=iT_num_bs

  ALLOCATE(T_hbs_start(num_atoms+1),T_bs_start(num_atoms+1))
  ALLOCATE(T_hbs_ind(T_num_hbs),T_bs_ind(T_num_bs))
  ALLOCATE(T_hbs_num(num_atoms),T_bs_num(num_atoms))

  T_hbs_start(:)=iT_hbs_start(:)
  T_bs_start(:)=iT_bs_start(:)
  T_hbs_ind(:)=iT_hbs_ind(:)
  T_bs_ind(:)=iT_bs_ind(:)

  DO ii=1,num_atoms
     T_hbs_num(ii)=T_hbs_start(ii+1)-T_hbs_start(ii)
     T_bs_num(ii)=T_bs_start(ii+1)-T_bs_start(ii)
  END DO

END SUBROUTINE load_net

SUBROUTINE build_shell1st (core)
 
  IMPLICIT NONE
 
  INTEGER,INTENT(IN)::core
 
  INTEGER::nn,nn2,mm,mm2,ii,jj,gg,kk,ll
  INTEGER,DIMENSION(:),ALLOCATABLE::order_ats_1sh,symm_ats_1sh,num_ats_bonded_1sh
  INTEGER,DIMENSION(:),ALLOCATABLE::order_ats_bonded_1sh,symm_ats_bonded_1sh
  INTEGER,DIMENSION(:),ALLOCATABLE::order_nodes_1sh,order_nodes_bonded_1sh
  INTEGER,DIMENSION(:),ALLOCATABLE::boxi,boxi2

  IF (ALLOCATED(mss_ind_atoms)) DEALLOCATE(mss_ind_atoms)
  IF (ALLOCATED(mss_ind_nodes)) DEALLOCATE(mss_ind_nodes)
  IF (ALLOCATED(mss_symm))      DEALLOCATE(mss_symm)
  IF (ALLOCATED(mss))           DEALLOCATE(mss)

  nn=atomspernode(core)
  nn2=2*nn
 
  ALLOCATE(order_ats_1sh(nn),order_nodes_1sh(nn),symm_ats_1sh(nn),num_ats_bonded_1sh(nn2))

  CALL build_order_ats_shell1st(core,nn,order_ats_1sh,symm_ats_1sh)

  gg=0
  DO ii=1,nn
     jj=order_ats_1sh(ii)
     num_ats_bonded_1sh(gg+1)=T_hbs_num(jj)
     num_ats_bonded_1sh(gg+2)=T_bs_num(jj)
     gg=gg+2
  END DO

  mm=SUM(num_ats_bonded_1sh(:),DIM=1)
  
  ALLOCATE(order_ats_bonded_1sh(mm),order_nodes_bonded_1sh(mm),symm_ats_bonded_1sh(mm))

  gg=0
  mm=0
  symm_ats_bonded_1sh(:)=0
  DO ii=1,nn
     jj=order_ats_1sh(ii)
     kk=num_ats_bonded_1sh(gg+1)
     IF (kk>0) THEN
        ll=T_hbs_start(jj)+1
        IF (kk==1) THEN
           mm=mm+1
           order_ats_bonded_1sh(mm)=T_hbs_ind(ll)
        ELSE
           ALLOCATE(boxi(kk),boxi2(kk))
           boxi(:)=T_hbs_ind((ll):(ll+kk-1))
           !CALL build_order_ats_bonded_shell1st(core,kk,boxi,boxi2)
           boxi2(:)=0
           order_ats_bonded_1sh((mm+1):(mm+kk))=boxi(:)
           symm_ats_bonded_1sh((mm+1):(mm+kk))=boxi2(:)
           DEALLOCATE(boxi,boxi2)
           mm=mm+kk
        END IF
     END IF
     kk=num_ats_bonded_1sh(gg+2)
     IF (kk>0) THEN
        ll=T_bs_start(jj)+1
        IF (kk==1) THEN
           mm=mm+1
           order_ats_bonded_1sh(mm)=T_bs_ind(ll)
        ELSE
           ALLOCATE(boxi(kk),boxi2(kk))
           boxi(:)=T_bs_ind((ll):(ll+kk-1))
           !CALL build_order_ats_bonded_shell1st(core,kk,boxi,boxi2)
           boxi2(:)=0
           order_ats_bonded_1sh((mm+1):(mm+kk))=boxi(:)
           symm_ats_bonded_1sh((mm+1):(mm+kk))=boxi2(:)
           DEALLOCATE(boxi,boxi2)
           mm=mm+kk
        END IF
     END IF
     gg=gg+2
  END DO

  ALLOCATE(mss_ind_atoms(1+nn+nn2+mm),mss_symm(1+nn+nn2+mm))
  ALLOCATE(mss_ind_nodes(1+nn+nn2+mm),mss(1+nn+nn2+mm))

  !Translate
  order_nodes_1sh(:)        =atom2node(order_ats_1sh(:))
  order_nodes_bonded_1sh(:) =atom2node(order_ats_bonded_1sh(:))
  order_ats_1sh(:)          =trad2py_atom(order_ats_1sh(:))
  order_ats_bonded_1sh(:)   =trad2py_atom(order_ats_bonded_1sh(:))
  order_nodes_1sh(:)        =trad2py_node(order_nodes_1sh(:))
  order_nodes_bonded_1sh(:) =trad2py_node(order_nodes_bonded_1sh(:))


  !Build

  mss_symm(:)=0
  mss_ind_atoms(1)=nn
  mss_ind_nodes(1)=nn
  ii=2
  jj=1+nn
  mss_ind_atoms(ii:jj)=order_ats_1sh(:)
  mss_ind_nodes(ii:jj)=order_nodes_1sh(:)
  mss_symm(ii:jj)=symm_ats_1sh(:)
  ii=jj+1
  jj=jj+nn2
  mss_ind_atoms(ii:jj)=num_ats_bonded_1sh(:)
  mss_ind_nodes(ii:jj)=num_ats_bonded_1sh(:)
  ii=jj+1
  jj=jj+mm
  mss_ind_atoms(ii:jj)=order_ats_bonded_1sh(:)
  mss_ind_nodes(ii:jj)=order_nodes_bonded_1sh(:)
  mss_symm(ii:jj)=symm_ats_bonded_1sh(:)

  DEALLOCATE(order_ats_1sh,order_ats_bonded_1sh)
  DEALLOCATE(order_nodes_1sh,order_nodes_bonded_1sh)
  DEALLOCATE(symm_ats_1sh,symm_ats_bonded_1sh,num_ats_bonded_1sh)



END SUBROUTINE build_shell1st


SUBROUTINE build_order_ats_shell1st (core,num_ats,order,symm_ats_1sh)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::core,num_ats
  INTEGER,DIMENSION(num_ats),INTENT(OUT)::order,symm_ats_1sh

  INTEGER::ii,jj,kk,gg,ll,num_crits
  LOGICAL::interruptor

  symm_ats_1sh(:)=0

  jj=0
  DO ii=node_run_ats(core)+1,node_run_ats(core+1)
     jj=jj+1
     order(jj)=ii
  END DO
 
  num_crits=symm_ats_crits(core)


  IF (num_crits>0) THEN

     ii=symm_ats_start(core)
     jj=symm_ats_start(core+1)
     gg=jj-ii
     ALLOCATE(vecti_aux(gg))
     vecti_aux(:)=symm_ats((ii+1):jj)

     interruptor=.true.
     IF (interruptor.eqv..true.) CALL SORTBYNUMHBS_ATS_1SH(core,num_ats,order,interruptor,num_crits,gg)
     IF (interruptor.eqv..true.) CALL SORTBYNUMBS_ATS_1SH(core,num_ats,order,interruptor,num_crits,gg)

     IF (interruptor.eqv..true.) THEN
        gg=0
        DO ii=1,num_crits
           gg=gg+1
           kk=vecti_aux(gg)
           DO jj=1,kk
              gg=gg+1
              ll=vecti_aux(gg)
              symm_ats_1sh(ll)=ii
           END DO
        END DO
        DEALLOCATE(vecti_aux)
     END IF

  END IF

END SUBROUTINE build_order_ats_shell1st


!!#### SORTING 1ST SHELL:

SUBROUTINE SORTBYNUMHBS_ATS_1SH (core,num_ats,order,interruptor,num_crits,dim_vecti_aux)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::core,num_ats
  INTEGER,INTENT(INOUT)::dim_vecti_aux
  INTEGER,DIMENSION(num_ats),INTENT(INOUT)::order
  LOGICAL,INTENT(INOUT)::interruptor
  INTEGER,INTENT(INOUT)::num_crits

  INTEGER::ii,jj,gg
  INTEGER,DIMENSION(:),ALLOCATABLE::valores

  ALLOCATE(valores(dim_vecti_aux))

  gg=0
  DO ii=1,num_crits
     gg=gg+1
     idim=vecti_aux(gg)
     DO jj=1,idim
        gg=gg+1
        kk=vecti_aux(gg)
        ll=order(kk)
        valores(gg)=T_hbs_num(ll)
     END DO
  END DO

  CALL SORT_INT_ATS_1SH(num_ats,dim_vecti_aux,order,valores,interruptor,num_crits)

  DEALLOCATE(valores)

END SUBROUTINE SORTBYNUMHBS_ATS_1SH

SUBROUTINE SORTBYNUMBS_ATS_1SH (core,num_ats,order,interruptor,num_crits,dim_vecti_aux)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::core,num_ats
  INTEGER,INTENT(INOUT)::dim_vecti_aux
  INTEGER,DIMENSION(num_ats),INTENT(INOUT)::order
  LOGICAL,INTENT(INOUT)::interruptor
  INTEGER,INTENT(INOUT)::num_crits

  INTEGER::ii,jj,gg
  INTEGER,DIMENSION(:),ALLOCATABLE::valores

  ALLOCATE(valores(dim_vecti_aux))

  gg=0
  DO ii=1,num_crits
     gg=gg+1
     idim=vecti_aux(gg)
     DO jj=1,idim
        gg=gg+1
        kk=vecti_aux(gg)
        ll=order(kk)
        valores(gg)=T_bs_num(ll)
     END DO
  END DO

  CALL SORT_INT_ATS_1SH(num_ats,dim_vecti_aux,order,valores,interruptor,num_crits)

  DEALLOCATE(valores)

END SUBROUTINE SORTBYNUMBS_ATS_1SH


SUBROUTINE SORT_INT_ATS_1SH (num_ats,dim_vecti_aux,order,valores,interruptor,num_crits)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::num_ats
  INTEGER,INTENT(INOUT)::dim_vecti_aux
  INTEGER,DIMENSION(num_ats),INTENT(INOUT)::order
  INTEGER,DIMENSION(dim_vecti_aux),INTENT(IN)::valores
  LOGICAL,INTENT(INOUT)::interruptor
  INTEGER,INTENT(INOUT)::num_crits

  INTEGER::ii,jj,kk,ll,gg,idim,new_num_crits,tope
  INTEGER,DIMENSION(:),ALLOCATABLE::val_aux,ind_aux,order_aux
  INTEGER,DIMENSION(:),ALLOCATABLE::vals,inds
  INTEGER,DIMENSION(:),ALLOCATABLE::new_symm_aux,cajon
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro

  new_num_crits=0

  gg=0
  DO ii=1,num_crits
     gg=gg+1
     idim=vecti_aux(gg)
     ALLOCATE(val_aux(idim),ind_aux(idim),order_aux(idim))
     DO jj=1,idim
        gg=gg+1
        kk=vecti_aux(gg)
        ll=order(kk)
        ind_aux(jj)=kk
        order_aux(jj)=ll
        val_aux(jj)=valores(gg)
     END DO
     ALLOCATE(filtro(idim),vals(idim),inds(idim))
     filtro=.TRUE.
     DO jj=1,idim
        kk=MAXLOC(val_aux,DIM=1,MASK=filtro(:))
        filtro(kk)=.FALSE.
        ll=ind_aux(jj)
        inds(jj)=ll
        vals(jj)=val_aux(kk)
        order(ll)=order_aux(kk)
     END DO
     interruptor=.FALSE.
     DO jj=2,idim
        IF (vals(jj-1)==vals(jj)) THEN
           IF (interruptor==.FALSE.) THEN
              order_aux(1)=inds(jj-1)
              order_aux(2)=inds(jj)
              ll=2
              interruptor=.TRUE.
           ELSE
              ll=ll+1
              order_aux(ll)=jj
           END IF
        ELSE
           IF (interruptor==.TRUE.) THEN
              interruptor=.FALSE.
              new_num_crits=new_num_crits+1
              IF (new_num_crits==1) THEN
                 tope=ll+1
                 ALLOCATE(new_symm_aux(tope))
                 new_symm_aux(1)=ll
                 new_symm_aux(2:tope)=order_aux(1:ll)
              ELSE
                 ALLOCATE(cajon(tope))
                 cajon(:)=new_symm_aux(:)
                 DEALLOCATE(new_symm_aux)
                 ALLOCATE(new_symm_aux(tope+1+ll))
                 new_symm_aux(1:tope)=cajon(:)
                 new_symm_aux(tope+1)=ll
                 new_symm_aux((tope+2):(tope+ll))=order_aux(1:ll)
                 tope=tope+1+ll
                 DEALLOCATE(cajon)
              END IF
           END IF
        END IF
     END DO
     IF (interruptor==.TRUE.) THEN
        interruptor=.FALSE.
        new_num_crits=new_num_crits+1
        IF (new_num_crits==1) THEN
           tope=ll+1
           ALLOCATE(new_symm_aux(tope))
           new_symm_aux(1)=ll
           new_symm_aux(2:tope)=order_aux(1:ll)
        ELSE
           ALLOCATE(cajon(tope))
           cajon(:)=new_symm_aux(:)
           DEALLOCATE(new_symm_aux)
           ALLOCATE(new_symm_aux(tope+1+ll))
           new_symm_aux(1:tope)=cajon(:)
           tope=tope+1
           new_symm_aux(tope)=ll
           new_symm_aux((tope+1):(tope+ll))=order_aux(1:ll)
           tope=tope+ll
           DEALLOCATE(cajon)
        END IF
     END IF
     DEALLOCATE(val_aux,ind_aux,order_aux)
     DEALLOCATE(filtro,vals,inds)
  END DO

  num_crits=new_num_crits
  DEALLOCATE(vecti_aux)
  IF (num_crits==0) THEN
     interruptor=.FALSE.
     dim_vecti_aux=0
  ELSE
     interruptor=.TRUE.
     ALLOCATE(vecti_aux(tope))
     dim_vecti_aux=tope
     vecti_aux(:)=new_symm_aux(:)
     DEALLOCATE(new_symm_aux)
  END IF

END SUBROUTINE SORT_INT_ATS_1SH

SUBROUTINE build_order_ats_bonded_shell1st (core,num_ats,order,symm_ats_bonded_1sh)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::core,num_ats
  INTEGER,DIMENSION(num_ats),INTENT(INOUT)::order
  INTEGER,DIMENSION(num_ats),INTENT(OUT)::symm_ats_bonded_1sh

  INTEGER::ii,jj,kk,gg,ll,num_crits
  LOGICAL::interruptor

  symm_ats_bonded_1sh(:)=0
  !### POR HACER (CONTINUAR AQUI)
  num_crits=1

     ii=symm_ats_start(core)
     jj=symm_ats_start(core+1)
     gg=jj-ii
     ALLOCATE(vecti_aux(gg))
     vecti_aux(:)=symm_ats((ii+1):jj)

     interruptor=.true.
     IF (interruptor.eqv..true.) CALL SORTBYNUMHBS_ATS_1SH(core,num_ats,order,interruptor,num_crits,gg)
     IF (interruptor.eqv..true.) CALL SORTBYNUMBS_ATS_1SH(core,num_ats,order,interruptor,num_crits,gg)

     IF (interruptor.eqv..true.) THEN
        gg=0
        DO ii=1,num_crits
           gg=gg+1
           kk=vecti_aux(gg)
           DO jj=1,kk
              gg=gg+1
              ll=vecti_aux(gg)
              symm_ats_1sh(ll)=ii
           END DO
        END DO
        DEALLOCATE(vecti_aux)
     END IF

  END IF

END SUBROUTINE build_order_ats_bonded_shell1st


END MODULE GLOB
