MODULE GLOB

INTEGER::num_nodes,num_atoms
INTEGER,DIMENSION(:),ALLOCATABLE::node_run_ats,atom2node,trad2py_node,trad2py_atom,atomspernode

INTEGER::T_num_hbs,T_num_bs
INTEGER,DIMENSION(:),ALLOCATABLE::T_hbs_start,T_bs_start,T_hbs_ind,T_bs_ind
INTEGER,DIMENSION(:),ALLOCATABLE::T_hbs_num,T_bs_num

INTEGER,DIMENSION(:),ALLOCATABLE::mss_ind_atoms,mss_ind_nodes,mss_symm,mss

CONTAINS

SUBROUTINE load_topol(inode_run_ats,iatom2node,itrad2py_node,itrad2py_atom,inum_nodes,inum_atoms)

  IMPLICIT NONE
  INTEGER,INTENT(IN)::inum_nodes,inum_atoms
  INTEGER,DIMENSION(inum_nodes+1),INTENT(IN)::inode_run_ats
  INTEGER,DIMENSION(inum_atoms),INTENT(IN)::iatom2node
  INTEGER,DIMENSION(inum_nodes),INTENT(IN)::itrad2py_node
  INTEGER,DIMENSION(inum_atoms),INTENT(IN)::itrad2py_atom

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
  INTEGER,DIMENSION(:),ALLOCATABLE::order


  IF (ALLOCATED(mss_ind_atoms)) DEALLOCATE(mss_ind_atoms)
  IF (ALLOCATED(mss_ind_nodes)) DEALLOCATE(mss_ind_nodes)
  IF (ALLOCATED(mss_symm))      DEALLOCATE(mss_symm)
  IF (ALLOCATED(mss))           DEALLOCATE(mss)

  nn=atomspernode(core)
 
  ALLOCATE(order(nn))
  
  jj=0
  DO ii=node_run_ats(core)+1,node_run_ats(core+1)
     jj=jj+1
     order(jj)=ii
  END DO
 
  gg=1+nn
  ALLOCATE(mss_ind_atoms(gg+2*nn))
  mss_ind_atoms(1)=nn
  mss_ind_atoms(2:gg)=order(:)
  
  DO ii=1,nn
     jj=order(ii)
     mss_ind_atoms(gg+1)=T_hbs_num(jj)
     mss_ind_atoms(gg+2)=T_bs_num(jj)
     gg=gg+2
  END DO

  DEALLOCATE(order)
 
END SUBROUTINE build_shell1st


END MODULE GLOB
