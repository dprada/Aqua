MODULE GLOB

INTEGER::num_nodes,num_atoms
INTEGER,DIMENSION(:),ALLOCATABLE::node_run_ats,atom2node,trad2py_node,trad2py_atom,atomspernode
INTEGER,DIMENSION(:),ALLOCATABLE::symm_ats_crits,symm_ats_start,symm_ats
INTEGER,DIMENSION(:),ALLOCATABLE::symm_aux

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
 
  INTEGER::nn,ii,jj,gg
  INTEGER,DIMENSION(:),ALLOCATABLE::order_ats_1sh,symm_ats_1sh


  IF (ALLOCATED(mss_ind_atoms)) DEALLOCATE(mss_ind_atoms)
  IF (ALLOCATED(mss_ind_nodes)) DEALLOCATE(mss_ind_nodes)
  IF (ALLOCATED(mss_symm))      DEALLOCATE(mss_symm)
  IF (ALLOCATED(mss))           DEALLOCATE(mss)

  nn=atomspernode(core)
 
  ALLOCATE(order_ats_1sh(nn),symm_ats_1sh(nn))
  
  CALL build_order_ats_shell1st(core,nn,order_ats_1sh,symm_ats_1sh)

  gg=1+nn
  ALLOCATE(mss_ind_atoms(gg+2*nn),mss_symm(gg+2*nn))
  mss_symm(:)=0
  mss_ind_atoms(1)=nn
  mss_ind_atoms(2:gg)=order_ats_1sh(:)
  
  DO ii=1,nn
     jj=order_ats_1sh(ii)
     mss_symm(1+ii)=symm_ats_1sh(ii)
     mss_ind_atoms(gg+1)=T_hbs_num(jj)
     mss_ind_atoms(gg+2)=T_bs_num(jj)
     gg=gg+2
  END DO


  DEALLOCATE(order_ats_1sh,symm_ats_1sh)
 
END SUBROUTINE build_shell1st


SUBROUTINE build_order_ats_shell1st (core,num_ats,order,symm_ats_1sh)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::core,num_ats
  INTEGER,DIMENSION(num_ats),INTENT(OUT)::order,symm_ats_1sh

  INTEGER::ii,jj,gg,num_crits
  LOGICAL::interruptor

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
     ALLOCATE(symm_aux(gg))
     symm_aux(:)=symm_ats((ii+1):jj)

     print*,symm_aux(:)

     interruptor=.true.

     !CALL averaver(symm)
     ! 
     !print*,symm(:)

!!$     ALLOCATE(filtro(num_ats))
!!$     
!!$     ggg=symm_ats_start(core)
!!$
!!$     DO ii=1,num_crits
!!$        
!!$        DO jj=1,num_ats
!!$           ggg=ggg+1
!!$           filtro(jj)=symm_ats(ggg)
!!$        END DO
!!$        interruptor=.true.
!!$
!!$        IF (interruptor.eqv..true.) CALL SORTBYNUMHBS_ATS_1SH(core,num_ats,order,filtro,interruptor)
!!$     
!!$
!!$        DO jj=1,num_ats
!!$           IF (filtro(jj).eqv..true.) THEN
!!$              symm(jj)=symm(jj)+1
!!$           END IF
!!$        END DO
!!$
!!$     END DO

     DEALLOCATE(symm_aux)

  END IF

END SUBROUTINE build_order_ats_shell1st

!!#### SORTING:

!!$SUBROUTINE SORTINTARRAY_1SH (quedan,vect_int,holes,num_ats,order,filtro)
!!$
!!$  INTEGER,INTENT(IN)::quedan,num_ats
!!$  INTEGER,DIMENSION(quedan),INTENT(IN)::vect_int,holes
!!$  INTEGER,DIMENSION(num_ats),INTENT(INOUT)::order
!!$  LOGICAL,DIMENSION(num_ats),INTENT(INOUT)::filtro
!!$
!!$  LOGICAL,DIMENSION(:),ALLOCATABLE::new_filtro
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::new_order,new_vals
!!$
!!$  ALLOCATE(new_filtro(quedan),new_vals(quedan),new_order(quedan))
!!$
!!$  DO ii=1,quedan
!!$     jj=MAXLOC(vect_int,DIM=1,MASK=new_filtro)
!!$     new_filtro(jj)=.FALSE.     
!!$     new_order(ii)=order(holes(jj))
!!$     new_vals(ii)=vect_int(jj)
!!$  END DO
!!$  DO ii=1,quedan-1
!!$     IF (new_vals(ii)==new_vals(ii+1)) THEN
!!$        new_filtro(ii)=.TRUE.
!!$        new_filtro(ii+1)=.TRUE.
!!$     END IF
!!$  END DO
!!$
!!$  DO ii=1,quedan
!!$     jj=holes(ii)
!!$     order(jj)=new_order(ii)
!!$     IF (new_vals(ii)/=new_vals(ii+1)) THEN
!!$        filtro(jj)=.FALSE.
!!$
!!$
!!$  
!!$
!!$SUBROUTINE SORTINTARRAY_1SH (vect_int,num_ats,order,filtro)
!!$
!!$  INTEGER,INTENT(IN)::core,num_ats
!!$  INTEGER,DIMENSION(),INTENT(IN)::vect_int
!!$  INTEGER,DIMENSION(num_ats),INTENT(INOUT)::order
!!$  LOGICAL,DIMENSION(num_ats),INTENT(INOUT)::filtro
!!$
!!$  INTEGER:: ii,jj,gg,val
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::sorted,holes,new_order
!!$  LOGICAL,DIMENSION(:),ALLOCATABLE::new_filtro
!!$
!!$  gg=COUNT(filtro)
!!$
!!$  ALLOCATE(sorted(gg),holes(gg),new_filtro(num_ats),new_order(num_ats))
!!$
!!$  gg=0
!!$  DO ii=1,num_ats
!!$     IF (filtro(ii).eqv..true.) THEN
!!$        gg=gg+1
!!$        holes(gg)=ii
!!$     END IF
!!$  END DO
!!$
!!$  new_order(:)=order(:)
!!$  new_filtro(:)=.FALSE.
!!$  jj=MAXLOC(vect_int,DIM=1,MASK=filtro)
!!$  val_old=vect_int(jj)
!!$  ind_old=holes(1)
!!$  new_order(ind_old)=order(jj)
!!$  filtro(jj)=.FALSE.
!!$
!!$  DO ii=2,gg
!!$     jj=MAXLOC(vect_int,DIM=1,MASK=filtro)
!!$     val_new=vect_int(jj)
!!$     ind_new=holes(jj)
!!$     new_order(ind_new)=order(jj)
!!$     filtro(jj)=.FALSE.
!!$     IF (val_new==val_old) THEN
!!$        new_filtro(ind_new)=.TRUE.
!!$        new_filtro(ind_old)=.TRUE.
!!$     END IF
!!$     val_old=val_new
!!$     ind_old=ind_new
!!$  END DO
!!$
!!$  filtro(:)=new_filtro(:)
!!$  IF (COUNT(filtro)==0) interruptor=.FALSE.
!!$
!!$  DEALLOCATE(sorted,holes,new_filtro)
!!$
!!$END SUBROUTINE SORTINTARRAY_1SH
!!$
!!$
!!$SUBROUTINE SORTBYNUMHBS_ATS_1SH (core,num_ats,order,filtro,interruptor,quedan)
!!$
!!$  INTEGER,INTENT(IN)::core,num_ats
!!$  INTEGER,DIMENSION(num_ats),INTENT(INOUT)::order
!!$  LOGICAL,DIMENSION(num_ats),INTENT(INOUT)::filtro
!!$  LOGICAL,INTENT(INOUT)::interruptor
!!$  INTEGER,INTENT(INOUT)::quedan
!!$
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::val_aux,holes
!!$
!!$  ALLOCATE(val_aux(quedan),)
!!$
!!$  gg=0
!!$  DO ii=1,num_ats
!!$     IF (filtro(ii).eqv..TRUE.) THEN
!!$        gg=gg+1
!!$        val_aux(gg)=T_hbs_num(order(ii))
!!$        holes(gg)=ii
!!$     END IF
!!$  END DO
!!$
!!$  CALL SORTINTARRAY_1SH(quedan,val_aux,holes,num_ats,order,filtro)
!!$
!!$  quedan=COUNT(filtro)
!!$  IF (quedan==0) interruptor=.FALSE.


END MODULE GLOB
