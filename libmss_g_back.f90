MODULE GLOB

INTEGER::definition_hbs
INTEGER::num_tot_hbs,num_tot_bs
INTEGER,DIMENSION(:,:),ALLOCATABLE::list_tot_hbs,list_tot_bs

CONTAINS

SUBROUTINE support_up (num_hbs,num_bs)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::num_hbs,num_bs

  IF (ALLOCATED(list_tot_hbs))     DEALLOCATE(list_tot_hbs)
  IF (ALLOCATED(list_tot_bs))      DEALLOCATE(list_tot_bs)

  ALLOCATE(list_tot_hbs(num_hbs,2),list_tot_bs(num_bs,2))

  num_tot_hbs=0 
  num_tot_bs=0

END SUBROUTINE support_up


SUBROUTINE add_hb_support(at1,at2)

  INTEGER,INTENT(IN)::at1,at2

  num_tot_hbs=num_tot_hbs+1
  list_tot_hbs(num_tot_hbs,:)=(/at1,at2/)

END SUBROUTINE add_hb_support

SUBROUTINE add_b_support(at1,at2)

  INTEGER,INTENT(IN)::at1,at2

  num_tot_bs=num_tot_bs+1
  list_tot_bs(num_tot_bs,:)=(/at1,at2/)


END SUBROUTINE add_b_support

SUBROUTINE support_down(atom2node,category,sets_per_set,&
  & nodes_per_set,num_categories,num_atoms,num_nodes,num_sets,codes_atom,codes_node,code_atom)

  INTEGER,INTENT(IN)::num_atoms,num_nodes,num_sets,num_categories
  INTEGER,DIMENSION(num_atoms),INTENT(IN)::atom2node
  INTEGER,DIMENSION(num_sets),INTENT(IN)::sets_per_set,nodes_per_set
  INTEGER,DIMENSION(num_nodes,3),INTENT(IN)::category
  INTEGER,DIMENSION(num_atoms,4),INTENT(OUT)::codes_atom
  INTEGER,DIMENSION(num_nodes,4),INTENT(OUT)::codes_node
  INTEGER,DIMENSION(num_atoms),INTENT(OUT)::code_atom

  INTEGER::ii,jj,kk,ll,at1,at2,nd1,nd2
  INTEGER::dim,dim_sets,top_sets
  INTEGER,DIMENSION(3)::cat
  INTEGER,DIMENSION(:),ALLOCATABLE::offset_sets
  LOGICAL,DIMENSION(:,:,:),ALLOCATABLE::aux_node_set_hb,aux_node_set_b
  LOGICAL,DIMENSION(:,:,:),ALLOCATABLE::aux_atom_set_hb,aux_atom_set_b
  INTEGER,DIMENSION(:,:),ALLOCATABLE  ::sup_atom_hb_1sh,sup_node_hb_1sh
  INTEGER,DIMENSION(:,:),ALLOCATABLE  ::sup_atom_b_1sh,sup_node_b_1sh
  INTEGER,DIMENSION(:,:),ALLOCATABLE  ::sup_atom_hb_2sh,sup_node_hb_2sh
  INTEGER,DIMENSION(:,:),ALLOCATABLE  ::sup_atom_b_2sh,sup_node_b_2sh

  dim=1+num_categories+num_sets
  DO ii=1,num_sets
     dim=dim+nodes_per_set(ii)
  END DO
  dim_sets=MAXVAL(sets_per_set)
  top_sets=num_categories-num_sets

  ALLOCATE(aux_node_set_hb(num_nodes,num_sets,dim_sets))
  ALLOCATE(aux_node_set_b (num_nodes,num_sets,dim_sets))
  ALLOCATE(aux_atom_set_hb(num_atoms,num_sets,dim_sets))
  ALLOCATE(aux_atom_set_b (num_atoms,num_sets,dim_sets))

  ALLOCATE(offset_sets(num_sets))

  offset_sets(1)=1+num_categories+num_sets
  DO ii=2,num_sets
     offset_sets(ii)=offset_sets(ii-1)+nodes_per_set(ii-1)
  END DO

  ALLOCATE(sup_atom_hb_1sh (num_atoms,dim))
  ALLOCATE(sup_atom_b_1sh  (num_atoms,dim)) 
  ALLOCATE(sup_node_hb_1sh (num_nodes,dim))
  ALLOCATE(sup_node_b_1sh  (num_nodes,dim)) 
  ALLOCATE(sup_atom_hb_2sh (num_atoms,dim))
  ALLOCATE(sup_atom_b_2sh  (num_atoms,dim)) 
  ALLOCATE(sup_node_hb_2sh (num_nodes,dim))
  ALLOCATE(sup_node_b_2sh  (num_nodes,dim)) 

  aux_node_set_hb=.FALSE.
  aux_node_set_b=.FALSE.
  aux_atom_set_hb=.FALSE.
  aux_atom_set_b=.FALSE.
  sup_atom_hb_1sh=0
  sup_atom_b_1sh=0
  sup_node_hb_1sh=0
  sup_node_b_1sh=0
  sup_atom_hb_2sh=0
  sup_atom_b_2sh=0
  sup_node_hb_2sh=0
  sup_node_b_2sh=0

 
  DO ii=1,num_tot_hbs

     at1=list_tot_hbs(ii,1)
     at2=list_tot_hbs(ii,2)
     nd1=atom2node(at1)
     nd2=atom2node(at2)

     !Num_hbs
     sup_atom_hb_1sh(at1,1)=sup_atom_hb_1sh(at1,1)+1
 
     !Num_nodes per type
     cat=category(nd2,:)+1
     jj=1+cat(1)
     sup_atom_hb_1sh(at1,jj)=sup_atom_hb_1sh(at1,jj)+1
     jj=jj-1
     IF (jj>top_sets) THEN
        kk=jj-top_sets
        aux_atom_set_hb(at1,kk,cat(2))=.TRUE.
        aux_node_set_hb(nd1,kk,cat(2))=.TRUE.
        ll=offset_sets(kk)+cat(3)
        sup_atom_hb_1sh(at1,ll)=sup_atom_hb_1sh(at1,ll)+1
     END IF

  END DO

  DO ii=1,num_tot_bs

     at1=list_tot_bs(ii,1)
     at2=list_tot_bs(ii,2)
     nd1=atom2node(at1)
     nd2=atom2node(at2)

     !Num_hbs
     sup_atom_b_1sh(at1,1)=sup_atom_b_1sh(at1,1)+1
 
     !Num_nodes per type
     cat=category(nd2,:)+1
     jj=1+cat(1)
     sup_atom_b_1sh(at1,jj)=sup_atom_b_1sh(at1,jj)+1
     jj=jj-1
     IF (jj>top_sets) THEN
        kk=jj-top_sets
        aux_atom_set_b(at1,kk,cat(2))=.TRUE.
        aux_node_set_b(nd1,kk,cat(2))=.TRUE.
        ll=offset_sets(kk)+cat(3)
        sup_atom_b_1sh(at1,ll)=sup_atom_b_1sh(at1,ll)+1
     END IF

  END DO

  DO ii=1,num_atoms
     ll=atom2node(ii)
     DO jj=1,num_sets
        kk=1+num_categories+jj
        sup_atom_hb_1sh(ii,kk)=COUNT(aux_atom_set_hb(ii,jj,:),DIM=1)
        sup_atom_b_1sh(ii,kk) =COUNT(aux_atom_set_b(ii,jj,:),DIM=1)
     END DO
     sup_node_hb_1sh(ll,:)= sup_node_hb_1sh(ll,:)+sup_atom_hb_1sh(ii,:)
     sup_node_b_1sh(ll,:) = sup_node_b_1sh(ll,:)+sup_atom_b_1sh(ii,:)
  END DO

  DO ii=1,num_nodes
     DO jj=1,num_sets
        kk=1+num_categories+jj
        sup_node_hb_1sh(ii,kk)=COUNT(aux_node_set_hb(ii,jj,:),DIM=1)
        sup_node_b_1sh(ii,kk) =COUNT(aux_node_set_b(ii,jj,:),DIM=1)
     END DO
  END DO

  aux_atom_set_hb=.FALSE.
  DO ii=1,num_tot_hbs

     at1=list_tot_hbs(ii,1)
     at2=list_tot_hbs(ii,2)
     nd1=atom2node(at1)
     nd2=atom2node(at2)

     sup_atom_hb_2sh(at1,:)=sup_atom_hb_2sh(at1,:)+sup_node_hb_1sh(nd2,:)
     DO jj=1,num_sets
        aux_atom_set_hb(at1,jj,:)=aux_atom_set_hb(at1,jj,:).OR.aux_node_set_hb(nd2,jj,:)
     END DO

  END DO

  aux_atom_set_b=.FALSE.
  DO ii=1,num_tot_bs

     at1=list_tot_bs(ii,1)
     at2=list_tot_bs(ii,2)
     nd1=atom2node(at1)
     nd2=atom2node(at2)

     sup_atom_b_2sh(at1,:)=sup_atom_b_2sh(at1,:)+sup_node_b_1sh(nd2,:)
     DO jj=1,num_sets
        aux_atom_set_b(at1,jj,:)=aux_atom_set_b(at1,jj,:).OR.aux_node_set_b(nd2,jj,:)
     END DO

  END DO

  aux_node_set_hb=.FALSE.
  aux_node_set_b=.FALSE.
  DO ii=1,num_atoms
     ll=atom2node(ii)
     DO jj=1,num_sets
        kk=1+num_categories+jj
        sup_atom_hb_2sh(ii,kk)=COUNT(aux_atom_set_hb(ii,jj,:),DIM=1)
        sup_atom_b_2sh(ii,kk) =COUNT(aux_atom_set_b(ii,jj,:),DIM=1)
        aux_node_set_hb(ll,jj,:)=aux_node_set_hb(ll,jj,:).OR.aux_atom_set_hb(ii,jj,:)
        aux_node_set_b(ll,jj,:) =aux_node_set_b(ll,jj,:).OR.aux_atom_set_b(ii,jj,:)
     END DO
     sup_node_hb_2sh(ll,:)= sup_node_hb_2sh(ll,:)+sup_atom_hb_2sh(ii,:)
     sup_node_b_2sh(ll,:) = sup_node_b_2sh(ll,:)+sup_atom_b_2sh(ii,:)
  END DO

  DO ii=1,num_nodes
     DO jj=1,num_sets
        kk=1+num_categories+jj
        sup_node_hb_2sh(ii,kk)=COUNT(aux_node_set_hb(ii,jj,:),DIM=1)
        sup_node_b_2sh(ii,kk) =COUNT(aux_node_set_b(ii,jj,:),DIM=1)
     END DO
  END DO

  !codigos

  CALL order_support (sup_atom_hb_1sh,num_atoms,dim,codes_atom(:,1))
  CALL order_support (sup_atom_b_1sh,num_atoms,dim, codes_atom(:,2))
  CALL order_support (sup_atom_hb_2sh,num_atoms,dim,codes_atom(:,3))
  CALL order_support (sup_atom_b_2sh,num_atoms,dim, codes_atom(:,4))
  CALL order_support (codes_atom,num_atoms,4,code_atom)

  CALL order_support (sup_node_hb_1sh,num_nodes,dim,codes_node(:,1))
  CALL order_support (sup_node_b_1sh,num_nodes,dim, codes_node(:,2))
  CALL order_support (sup_node_hb_2sh,num_nodes,dim,codes_node(:,3))
  CALL order_support (sup_node_b_2sh,num_nodes,dim, codes_node(:,4))


  DEALLOCATE(aux_node_set_hb,aux_node_set_b)
  DEALLOCATE(aux_atom_set_hb,aux_atom_set_b)
  DEALLOCATE(sup_atom_hb_1sh,sup_atom_b_1sh)
  DEALLOCATE(sup_atom_hb_2sh,sup_atom_b_2sh)
  DEALLOCATE(sup_node_hb_1sh,sup_node_b_1sh)
  DEALLOCATE(sup_node_hb_2sh,sup_node_b_2sh)
  DEALLOCATE(offset_sets)

END SUBROUTINE support_down


SUBROUTINE order_support (box,num_obj,dim,order)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::num_obj,dim
  INTEGER,DIMENSION(num_obj,dim),INTENT(IN)::box
  INTEGER,DIMENSION(num_obj),INTENT(OUT)::order

  INTEGER::ii,jj,kk,factor
  INTEGER,DIMENSION(:),ALLOCATABLE::trad
  LOGICAL,DIMENSION(:),ALLOCATABLE::veo

  order(:)=0

  DO ii=1,dim
     factor=MAXVAL(box(:,ii),DIM=1)
     order(:)=order(:)*(factor+1)+box(:,ii)
     factor=MAXVAL(order(:),DIM=1)
     ALLOCATE(trad(0:factor),veo(0:factor))
     veo=.FALSE.
     DO jj=1,num_obj
        veo(order(jj))=.TRUE.
     END DO
     kk=0
     DO jj=0,factor
        IF (veo(jj).eqv..TRUE.) THEN
           trad(jj)=kk
           kk=kk+1
        END IF
     END DO
     DO jj=1,num_obj
        order(jj)=trad(order(jj))
     END DO
     DEALLOCATE(trad,veo)
  END DO

END SUBROUTINE order_support



!!SUBROUTINE breaking_symmmetry_centrality(ind_atoms,ind_nodes,symm,center,len_mss,new_ind_atoms,new_ind_nodes,new_symm)
!! 
!!  INTEGER::len_mss,center
!!  INTEGER,DIMENSION(len_mss),INTENT(IN)::ind_atoms,ind_nodes,symm
!!  INTEGER,DIMENSION(len_mss),INTENT(OUT)::new_ind_atoms,new_ind_nodes,new_symm
!!  
!!  natoms=ind_nodes(1)
!!  ALLOCATE(nbs(natoms,2))
!! 
!!  jj=1
!!  DO ii=1,natoms
!!     nbs(ii,1)=ind_nodes(jj) 
!!     nbs(ii,2)=ind_nodes(jj+1) 
!!     jj=jj+2
!!  END DO
!! 
!!  
!! 
!! 
!!END SUBROUTINE breaking_symmmetry_centrality

SUBROUTINE breaking_symmetry_1st_2(criterium,orden,support,num_atoms,neworden,new_symm)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::num_atoms
  INTEGER,DIMENSION(num_atoms),INTENT(IN)::criterium
  INTEGER,DIMENSION(num_atoms),INTENT(IN)::orden
  INTEGER,DIMENSION(num_atoms),INTENT(IN)::support
  INTEGER,DIMENSION(num_atoms),INTENT(OUT)::neworden
  INTEGER,DIMENSION(num_atoms),INTENT(OUT)::new_symm

  INTEGER::ii,jj,gg,kk,remains
  INTEGER,DIMENSION(:),ALLOCATABLE::holes
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro


  ALLOCATE(holes(num_atoms))
  ALLOCATE(filtro(num_atoms))

  remains=0
  neworden(:)=orden(:)
  filtro(:)=.FALSE.
  DO ii=1,num_atoms
     IF (criterium(ii)==1) THEN
        filtro(ii)=.TRUE.
        remains=remains+1
        holes(remains)=ii
     END IF
  END DO

  new_symm(:)=0
  gg=-1
  DO ii=1,remains
     jj=MAXLOC(support,DIM=1,MASK=filtro)
     kk=support(jj)
     IF ((gg==kk).AND.(gg>0)) THEN
        new_symm(holes(ii-1))=1
        new_symm(holes(ii))=1
     END IF
     gg=kk
     filtro(jj)=.FALSE.
     neworden(holes(ii))=orden(jj)
  END DO

  DEALLOCATE(filtro,holes)


END SUBROUTINE breaking_symmetry_1st_2


SUBROUTINE breaking_symmetry_1st(criterium,orden,support,num_atoms,num_crit,neworden,new_symm)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::num_atoms,num_crit
  INTEGER,DIMENSION(num_atoms),INTENT(IN)::criterium
  INTEGER,DIMENSION(num_atoms),INTENT(IN)::orden
  INTEGER,DIMENSION(num_atoms,num_crit),INTENT(IN)::support
  INTEGER,DIMENSION(num_atoms),INTENT(OUT)::neworden
  INTEGER,DIMENSION(num_atoms),INTENT(OUT)::new_symm

  INTEGER::ii,jj,gg,kk,eff_num_crit,remains
  INTEGER,DIMENSION(:),ALLOCATABLE::potencias,valores,base,holes
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro
  INTEGER,DIMENSION(:,:),ALLOCATABLE::eff_support

  ALLOCATE(filtro(num_crit))
  filtro=.False.
  DO ii=1,num_crit
     IF (SUM(support(:,ii))>0) THEN
        filtro(ii)=.True.
     END IF
  END DO
  eff_num_crit=COUNT(filtro)

  ALLOCATE(eff_support(num_atoms,eff_num_crit),base(eff_num_crit),valores(num_atoms))

  eff_num_crit=0
  DO ii=1,num_crit
     IF (filtro(ii).eqv..True.) THEN
        eff_num_crit=eff_num_crit+1
        valores(:)=support(:,ii)
        eff_support(:,eff_num_crit)=valores(:)
        base(eff_num_crit)=MAXVAL(valores(:),DIM=1)+1
     END IF
  END DO

  DEALLOCATE(filtro)

  ALLOCATE(potencias(eff_num_crit))
  ALLOCATE(filtro(num_atoms))


  potencias(eff_num_crit)=1
  DO ii=(eff_num_crit-1),1,-1
     potencias(ii)=potencias(ii+1)*base(ii+1)
  END DO

  DO ii=1,num_atoms
     valores(ii)=dot_product(potencias,eff_support(ii,:))
  END DO

  ALLOCATE(holes(num_atoms))

  remains=0
  neworden(:)=orden(:)
  filtro(:)=.FALSE.
  DO ii=1,num_atoms
     IF (criterium(ii)==1) THEN
        filtro(ii)=.TRUE.
        remains=remains+1
        holes(remains)=ii
     END IF
  END DO

  new_symm(:)=0
  gg=-1
  DO ii=1,remains
     jj=MAXLOC(valores,DIM=1,MASK=filtro)
     kk=valores(jj)
     IF ((gg==kk).AND.(gg>0)) THEN
        new_symm(holes(ii-1))=1
        new_symm(holes(ii))=1
     END IF
     gg=kk
     filtro(jj)=.FALSE.
     neworden(holes(ii))=orden(jj)
  END DO

  DEALLOCATE(filtro,potencias,eff_support,base,valores,holes)

END SUBROUTINE breaking_symmetry_1st

SUBROUTINE breaking_symmetry_2nd(orden1,orden2,support,num_atoms,num_crit,neworden1,neworden2,new_symm)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::num_atoms,num_crit
  INTEGER,DIMENSION(num_atoms),INTENT(IN)::orden1,orden2
  INTEGER,DIMENSION(num_atoms,num_crit),INTENT(IN)::support
  INTEGER,DIMENSION(num_atoms),INTENT(OUT)::neworden1,neworden2
  INTEGER,DIMENSION(num_atoms),INTENT(OUT)::new_symm

  INTEGER::ii,jj,gg,kk,eff_num_crit
  INTEGER,DIMENSION(:),ALLOCATABLE::potencias,valores,base
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro
  INTEGER,DIMENSION(:,:),ALLOCATABLE::eff_support


  ALLOCATE(filtro(num_crit))
  filtro=.False.
  DO ii=1,num_crit
     IF (SUM(support(:,ii))>0) THEN
        filtro(ii)=.True.
     END IF
  END DO
  eff_num_crit=COUNT(filtro)

  ALLOCATE(eff_support(num_atoms,eff_num_crit),base(eff_num_crit),valores(num_atoms))

  eff_num_crit=0
  DO ii=1,num_crit
     IF (filtro(ii).eqv..True.) THEN
        eff_num_crit=eff_num_crit+1
        valores(:)=support(:,ii)
        eff_support(:,eff_num_crit)=valores(:)
        base(eff_num_crit)=MAXVAL(valores(:),DIM=1)+1
     END IF
  END DO

  DEALLOCATE(filtro)

  ALLOCATE(potencias(eff_num_crit))
  ALLOCATE(filtro(num_atoms))

  filtro(:)=.True.

  potencias(eff_num_crit)=1
  DO ii=(eff_num_crit-1),1,-1
     potencias(ii)=potencias(ii+1)*base(ii+1)
  END DO

  DO ii=1,num_atoms
     valores(ii)=dot_product(potencias,eff_support(ii,:))
  END DO

  new_symm(:)=0
  gg=-1
  DO ii=1,num_atoms
     jj=MAXLOC(valores,DIM=1,MASK=filtro)
     kk=valores(jj)
     IF (gg==kk) THEN
        new_symm(ii)=1
        new_symm(ii-1)=1
     END IF
     gg=kk
     filtro(jj)=.FALSE.
     neworden1(ii)=orden1(jj)
     neworden2(ii)=orden2(jj)
  END DO

  DEALLOCATE(filtro,potencias,eff_support,base,valores)

END SUBROUTINE breaking_symmetry_2nd



SUBROUTINE ind_wat_limit_4_nosim (mss,aux,hbs,hbdists,num_wats,num_atoms,num_hbs)


  IMPLICIT NONE

  INTEGER,INTENT(IN)::num_atoms,num_wats,num_hbs
  INTEGER,DIMENSION(num_atoms,3),INTENT(IN)::aux
  INTEGER,DIMENSION(num_hbs,3),INTENT(IN)::hbs
  DOUBLE PRECISION,DIMENSION(num_hbs),INTENT(IN)::hbdists
  INTEGER,DIMENSION(num_wats,17),INTENT(OUT)::mss
  
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::dist_first_shell
  INTEGER,DIMENSION(:,:),ALLOCATABLE::first_shell
  INTEGER,DIMENSION(:),ALLOCATABLE::bonds_o
  INTEGER,DIMENSION(:),ALLOCATABLE::microstate
  
  INTEGER::ii,jj,kk,bb,gg
  INTEGER::ind_oh,hi,ind_o
  
  mss=0
  ALLOCATE(dist_first_shell(num_wats,4),first_shell(num_wats,4),bonds_o(num_wats))
  ALLOCATE(microstate(17))

  first_shell=0
  dist_first_shell=0.0d0
  bonds_o=0


  IF (ANY((/2,3/)==definition_hbs)) THEN
     DO ii=1,num_hbs
        ind_oh=aux(hbs(ii,1)+1,1)+1
        hi=aux(hbs(ii,2)+1,2)
        ind_o=aux(hbs(ii,3)+1,1)+1
        IF (first_shell(ind_oh,hi)==0) THEN
           first_shell(ind_oh,hi)=ind_o
           dist_first_shell(ind_oh,hi)=hbdists(ii)
        ELSE
           IF (dist_first_shell(ind_oh,hi)>hbdists(ii)) THEN
              first_shell(ind_oh,hi)=ind_o
              dist_first_shell(ind_oh,hi)=hbdists(ii)
           END IF
        END IF
        IF (bonds_o(ind_o)==0) THEN
           first_shell(ind_o,3)=ind_oh
           dist_first_shell(ind_o,3)=hbdists(ii)
           bonds_o(ind_o)=1
        ELSE
           IF (bonds_o(ind_o)==1) THEN
              IF (dist_first_shell(ind_o,3)<hbdists(ii)) THEN
                 first_shell(ind_o,4)=ind_oh
                 dist_first_shell(ind_o,4)=hbdists(ii)
              ELSE
                 first_shell(ind_o,4)=first_shell(ind_o,3)
                 dist_first_shell(ind_o,4)=dist_first_shell(ind_o,3)
                 first_shell(ind_o,3)=ind_oh
                 dist_first_shell(ind_o,3)=hbdists(ii)
              END IF
              bonds_o(ind_o)=2
           ELSE
              IF (bonds_o(ind_o)==2) THEN
                 IF (dist_first_shell(ind_o,3)>hbdists(ii)) THEN
                    first_shell(ind_o,4)=first_shell(ind_o,3)
                    dist_first_shell(ind_o,4)=dist_first_shell(ind_o,3)
                    first_shell(ind_o,3)=ind_oh
                    dist_first_shell(ind_o,3)=hbdists(ii)
                 ELSE
                    IF (dist_first_shell(ind_o,4)>hbdists(ii)) THEN
                       first_shell(ind_o,4)=ind_oh
                       dist_first_shell(ind_o,4)=hbdists(ii)
                    END IF
                 END IF
              END IF
           END IF
        END IF
     END DO
  ELSE IF (ANY((/1/)==definition_hbs)) THEN
     DO ii=1,num_hbs
        ind_oh=aux(hbs(ii,1)+1,1)+1
        hi=aux(hbs(ii,2)+1,2)
        ind_o=aux(hbs(ii,3)+1,1)+1
        IF (first_shell(ind_oh,hi)==0) THEN
           first_shell(ind_oh,hi)=ind_o
           dist_first_shell(ind_oh,hi)=hbdists(ii)
        ELSE
           IF (dist_first_shell(ind_oh,hi)<hbdists(ii)) THEN
              first_shell(ind_oh,hi)=ind_o
              dist_first_shell(ind_oh,hi)=hbdists(ii)
           END IF
        END IF
        IF (bonds_o(ind_o)==0) THEN
           first_shell(ind_o,3)=ind_oh
           dist_first_shell(ind_o,3)=hbdists(ii)
           bonds_o(ind_o)=1
        ELSE
           IF (bonds_o(ind_o)==1) THEN
              IF (dist_first_shell(ind_o,3)>hbdists(ii)) THEN
                 first_shell(ind_o,4)=ind_oh
                 dist_first_shell(ind_o,4)=hbdists(ii)
              ELSE
                 first_shell(ind_o,4)=first_shell(ind_o,3)
                 dist_first_shell(ind_o,4)=dist_first_shell(ind_o,3)
                 first_shell(ind_o,3)=ind_oh
                 dist_first_shell(ind_o,3)=hbdists(ii)
              END IF
              bonds_o(ind_o)=2
           ELSE
              IF (bonds_o(ind_o)==2) THEN
                 IF (dist_first_shell(ind_o,3)<hbdists(ii)) THEN
                    first_shell(ind_o,4)=first_shell(ind_o,3)
                    dist_first_shell(ind_o,4)=dist_first_shell(ind_o,3)
                    first_shell(ind_o,3)=ind_oh
                    dist_first_shell(ind_o,3)=hbdists(ii)
                 ELSE
                    IF (dist_first_shell(ind_o,4)<hbdists(ii)) THEN
                       first_shell(ind_o,4)=ind_oh
                       dist_first_shell(ind_o,4)=hbdists(ii)
                    END IF
                 END IF
              END IF
           END IF
        END IF
     END DO
  END IF

  DO kk=1,num_wats

     microstate=0
     microstate(1)=kk
     
     microstate(2:5)=(/ first_shell(kk,:) /)

     gg=6
     DO ii=2,3
        bb=microstate(ii)
        IF (bb>0) THEN
           microstate(gg:gg+1)=(/ first_shell(bb,1:2) /)
           IF ((first_shell(bb,3)==kk).or.(first_shell(bb,4)==kk)) THEN
              IF (first_shell(bb,3)==kk) THEN
                 microstate(gg+2)=first_shell(bb,4)
              ELSE
                 microstate(gg+2)=first_shell(bb,3)
              END IF
              IF ((first_shell(bb,3)==kk).and.(first_shell(bb,4)==kk)) THEN
                 !print*,'problema ssn1' !!! Comprobar esto
                 microstate(gg+2)=first_shell(bb,3)
              END IF
           ELSE
              microstate(gg+2)=first_shell(bb,3)
           END IF
           gg=gg+3
        ELSE
           microstate(gg:gg+2)=0
           gg=gg+3
        END IF
     END DO

     DO ii=4,5
        bb=microstate(ii)
        IF (bb>0) THEN
           microstate(gg+1:gg+2)=(/ first_shell(bb,3:4) /)
           IF ((first_shell(bb,1)==kk).or.(first_shell(bb,2)==kk)) THEN
              IF (first_shell(bb,1)==kk) THEN
                 microstate(gg)=first_shell(bb,2)
              ELSE
                 microstate(gg)=first_shell(bb,1)
              END IF
              IF ((first_shell(bb,1)==kk).and.(first_shell(bb,2)==kk)) THEN
                 ! print*,'pproblema ssn2' !!! Comprobar esto
                 ! print*,a,'|',first_shell(b,:)
                 microstate(gg)=first_shell(bb,1)
                 ! print'(I5,A,4I5,A,3I5,A,3I5,A,3I5,A,3I5)',microstate(1),'||',microstate(2:5),'||',microstate(6:8),&
                 ! & '||',microstate(9:11),'||',microstate(12:14),'||',microstate(15:17)
              END IF
           ELSE
              microstate(gg)=first_shell(bb,1) !!!! esto esta mal, hay que elegir entre el 1 y el 2
              !print*,'siiiiiiiiiiiiiiiiii'
           END IF
           gg=gg+3
        ELSE
           microstate(gg:gg+2)=0
           gg=gg+3
        END IF
     END DO
     
     !print'(I5,A,4I5,A,3I5,A,3I5,A,3I5,A,3I5)',microstate(1),'||',microstate(2:5),'||',microstate(6:8),&
     ! & '||',microstate(9:11),'||',microstate(12:14),'||',microstate(15:17)

     !Quito indices de mols

     !!filtro=.true.
     !!aux=microstate
     !!mss_ind_wat=microstate
     !!DO i=1,17
     !!   IF (aux(i)<=0) filtro(i)=.false.
     !!END DO
     !! 
     !!shell_w(a,:)=0
     !!g=0
     !! 
     !!DO j=1,17
     !!   IF (filtro(j).eqv..true.) THEN
     !!      b=aux(j)
     !!      g=g+1
     !!      shell_w(a,g)=b
     !!      DO i=1,17
     !!         IF (filtro(i).eqv..true.) THEN
     !!            IF (aux(i)==b) THEN
     !!               microstate(i)=j
     !!               filtro(i)=.false.
     !!            END IF
     !!         END IF
     !!      END DO
     !!   END IF
     !!END DO

     mss(kk,:)=microstate

  END DO


  DEALLOCATE(dist_first_shell,first_shell,bonds_o)
  DEALLOCATE(microstate)


END SUBROUTINE ind_wat_limit_4_nosim

SUBROUTINE ind_wat_limit_4_nosim_prot (mss,aux,filt_water,hbs,hbdists,num_wats,num_atoms,num_hbs)


  IMPLICIT NONE

  INTEGER,INTENT(IN)::num_atoms,num_wats,num_hbs
  INTEGER,DIMENSION(num_atoms,3),INTENT(IN)::aux
  LOGICAL,DIMENSION(num_atoms),INTENT(IN)::filt_water
  INTEGER,DIMENSION(num_hbs,3),INTENT(IN)::hbs
  DOUBLE PRECISION,DIMENSION(num_hbs),INTENT(IN)::hbdists
  INTEGER,DIMENSION(num_wats,17),INTENT(OUT)::mss
  
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::dist_first_shell
  INTEGER,DIMENSION(:,:),ALLOCATABLE::first_shell
  INTEGER,DIMENSION(:),ALLOCATABLE::bonds_o
  INTEGER,DIMENSION(:),ALLOCATABLE::microstate
  
  INTEGER::ii,jj,kk,bb,gg
  INTEGER::ind_oh,hi,ind_o
  
  LOGICAL::interr_oh,interr_o

  mss=0
  ALLOCATE(dist_first_shell(num_wats,4),first_shell(num_wats,4),bonds_o(num_wats))
  ALLOCATE(microstate(17))

  first_shell=0
  dist_first_shell=0.0d0
  bonds_o=0


  IF (ANY((/2,3/)==definition_hbs)) THEN
     DO ii=1,num_hbs
        interr_oh=.FALSE.
        interr_o=.FALSE.
        IF (filt_water(hbs(ii,1)+1)) THEN
           interr_oh=.TRUE.
           ind_oh=aux(hbs(ii,1)+1,1)+1
           hi=aux(hbs(ii,2)+1,2)
        ELSE
           ind_oh=-(hbs(ii,2)+1)
        END IF
        IF (filt_water(hbs(ii,3)+1)) THEN
           interr_o=.TRUE.
           ind_o=aux(hbs(ii,3)+1,1)+1
        ELSE
           ind_o=-(hbs(ii,3)+1)
        END IF
        IF (interr_oh.eqv..TRUE.) THEN
           IF (first_shell(ind_oh,hi)==0) THEN
              first_shell(ind_oh,hi)=ind_o
              dist_first_shell(ind_oh,hi)=hbdists(ii)
           ELSE
              IF (dist_first_shell(ind_oh,hi)>hbdists(ii)) THEN
                 first_shell(ind_oh,hi)=ind_o
                 dist_first_shell(ind_oh,hi)=hbdists(ii)
              END IF
           END IF
        END IF
        IF (interr_o.eqv..TRUE.) THEN
           IF (bonds_o(ind_o)==0) THEN
              first_shell(ind_o,3)=ind_oh
              dist_first_shell(ind_o,3)=hbdists(ii)
              bonds_o(ind_o)=1
           ELSE
              IF (bonds_o(ind_o)==1) THEN
                 IF (dist_first_shell(ind_o,3)<hbdists(ii)) THEN
                    first_shell(ind_o,4)=ind_oh
                    dist_first_shell(ind_o,4)=hbdists(ii)
                 ELSE
                    first_shell(ind_o,4)=first_shell(ind_o,3)
                    dist_first_shell(ind_o,4)=dist_first_shell(ind_o,3)
                    first_shell(ind_o,3)=ind_oh
                    dist_first_shell(ind_o,3)=hbdists(ii)
                 END IF
                 bonds_o(ind_o)=2
              ELSE
                 IF (bonds_o(ind_o)==2) THEN
                    IF (dist_first_shell(ind_o,3)>hbdists(ii)) THEN
                       first_shell(ind_o,4)=first_shell(ind_o,3)
                       dist_first_shell(ind_o,4)=dist_first_shell(ind_o,3)
                       first_shell(ind_o,3)=ind_oh
                       dist_first_shell(ind_o,3)=hbdists(ii)
                    ELSE
                       IF (dist_first_shell(ind_o,4)>hbdists(ii)) THEN
                          first_shell(ind_o,4)=ind_oh
                          dist_first_shell(ind_o,4)=hbdists(ii)
                       END IF
                    END IF
                 END IF
              END IF
           END IF
        END IF
     END DO
  ELSE IF (ANY((/1/)==definition_hbs)) THEN
     DO ii=1,num_hbs
        interr_oh=.FALSE.
        interr_o=.FALSE.
        IF (filt_water(hbs(ii,1)+1)) THEN
           interr_oh=.TRUE.
           ind_oh=aux(hbs(ii,1)+1,1)+1
           hi=aux(hbs(ii,2)+1,2)
        ELSE
           ind_oh=-(hbs(ii,2)+1)
        END IF
        IF (filt_water(hbs(ii,3)+1)) THEN
           interr_o=.TRUE.
           ind_o=aux(hbs(ii,3)+1,1)+1
        ELSE
           ind_o=-(hbs(ii,3)+1)
        END IF
        IF (interr_oh.eqv..TRUE.) THEN
           IF (first_shell(ind_oh,hi)==0) THEN
              first_shell(ind_oh,hi)=ind_o
              dist_first_shell(ind_oh,hi)=hbdists(ii)
           ELSE
              IF (dist_first_shell(ind_oh,hi)<hbdists(ii)) THEN
                 first_shell(ind_oh,hi)=ind_o
                 dist_first_shell(ind_oh,hi)=hbdists(ii)
              END IF
           END IF
        END IF
        IF (interr_o.eqv..TRUE.) THEN
           IF (bonds_o(ind_o)==0) THEN
              first_shell(ind_o,3)=ind_oh
              dist_first_shell(ind_o,3)=hbdists(ii)
              bonds_o(ind_o)=1
           ELSE
              IF (bonds_o(ind_o)==1) THEN
                 IF (dist_first_shell(ind_o,3)>hbdists(ii)) THEN
                    first_shell(ind_o,4)=ind_oh
                    dist_first_shell(ind_o,4)=hbdists(ii)
                 ELSE
                    first_shell(ind_o,4)=first_shell(ind_o,3)
                    dist_first_shell(ind_o,4)=dist_first_shell(ind_o,3)
                    first_shell(ind_o,3)=ind_oh
                    dist_first_shell(ind_o,3)=hbdists(ii)
                 END IF
                 bonds_o(ind_o)=2
              ELSE
                 IF (bonds_o(ind_o)==2) THEN
                    IF (dist_first_shell(ind_o,3)<hbdists(ii)) THEN
                       first_shell(ind_o,4)=first_shell(ind_o,3)
                       dist_first_shell(ind_o,4)=dist_first_shell(ind_o,3)
                       first_shell(ind_o,3)=ind_oh
                       dist_first_shell(ind_o,3)=hbdists(ii)
                    ELSE
                       IF (dist_first_shell(ind_o,4)<hbdists(ii)) THEN
                          first_shell(ind_o,4)=ind_oh
                          dist_first_shell(ind_o,4)=hbdists(ii)
                       END IF
                    END IF
                 END IF
              END IF
           END IF
        END IF
     END DO
  END IF

  DO kk=1,num_wats

     microstate=0
     microstate(1)=kk
     
     microstate(2:5)=(/ first_shell(kk,:) /)

     gg=6
     DO ii=2,3
        bb=microstate(ii)
        IF (bb>0) THEN
           microstate(gg:gg+1)=(/ first_shell(bb,1:2) /)
           IF ((first_shell(bb,3)==kk).or.(first_shell(bb,4)==kk)) THEN
              IF (first_shell(bb,3)==kk) THEN
                 microstate(gg+2)=first_shell(bb,4)
              ELSE
                 microstate(gg+2)=first_shell(bb,3)
              END IF
              IF ((first_shell(bb,3)==kk).and.(first_shell(bb,4)==kk)) THEN
                 !print*,'problema ssn1' !!! Comprobar esto
                 microstate(gg+2)=first_shell(bb,3)
              END IF
           ELSE
              microstate(gg+2)=first_shell(bb,3)
           END IF
           gg=gg+3
        ELSE
           microstate(gg:gg+2)=0
           gg=gg+3
        END IF
     END DO

     DO ii=4,5
        bb=microstate(ii)
        IF (bb>0) THEN
           microstate(gg+1:gg+2)=(/ first_shell(bb,3:4) /)
           IF ((first_shell(bb,1)==kk).or.(first_shell(bb,2)==kk)) THEN
              IF (first_shell(bb,1)==kk) THEN
                 microstate(gg)=first_shell(bb,2)
              ELSE
                 microstate(gg)=first_shell(bb,1)
              END IF
              IF ((first_shell(bb,1)==kk).and.(first_shell(bb,2)==kk)) THEN
                 ! print*,'pproblema ssn2' !!! Comprobar esto
                 ! print*,a,'|',first_shell(b,:)
                 microstate(gg)=first_shell(bb,1)
                 ! print'(I5,A,4I5,A,3I5,A,3I5,A,3I5,A,3I5)',microstate(1),'||',microstate(2:5),'||',microstate(6:8),&
                 ! & '||',microstate(9:11),'||',microstate(12:14),'||',microstate(15:17)
              END IF
           ELSE
              microstate(gg)=first_shell(bb,1) !!!! esto esta mal, hay que elegir entre el 1 y el 2
              !print*,'siiiiiiiiiiiiiiiiii'
           END IF
           gg=gg+3
        ELSE
           microstate(gg:gg+2)=0
           gg=gg+3
        END IF
     END DO
     
     !print'(I5,A,4I5,A,3I5,A,3I5,A,3I5,A,3I5)',microstate(1),'||',microstate(2:5),'||',microstate(6:8),&
     ! & '||',microstate(9:11),'||',microstate(12:14),'||',microstate(15:17)

     !Quito indices de mols

     !!filtro=.true.
     !!aux=microstate
     !!mss_ind_wat=microstate
     !!DO i=1,17
     !!   IF (aux(i)<=0) filtro(i)=.false.
     !!END DO
     !! 
     !!shell_w(a,:)=0
     !!g=0
     !! 
     !!DO j=1,17
     !!   IF (filtro(j).eqv..true.) THEN
     !!      b=aux(j)
     !!      g=g+1
     !!      shell_w(a,g)=b
     !!      DO i=1,17
     !!         IF (filtro(i).eqv..true.) THEN
     !!            IF (aux(i)==b) THEN
     !!               microstate(i)=j
     !!               filtro(i)=.false.
     !!            END IF
     !!         END IF
     !!      END DO
     !!   END IF
     !!END DO

     mss(kk,:)=microstate

  END DO


  DEALLOCATE(dist_first_shell,first_shell,bonds_o)
  DEALLOCATE(microstate)


END SUBROUTINE ind_wat_limit_4_nosim_prot


SUBROUTINE remove_index_mol(mss_ind,num_mss,mss)

  IMPLICIT NONE
  INTEGER,INTENT(IN)::num_mss
  INTEGER,DIMENSION(num_mss,17),INTENT(IN)::mss_ind
  INTEGER,DIMENSION(num_mss,17),INTENT(OUT)::mss

  INTEGER::ii,jj,gg,ll,bb
  INTEGER,DIMENSION(17)::aux_mss,aux_mss2
  LOGICAL,DIMENSION(17)::filtro


  DO jj=1,num_mss
     filtro=.true.
     aux_mss2=mss_ind(jj,:)
     aux_mss=aux_mss2

     DO ii=1,17
        IF (aux_mss(ii)<=0) filtro(ii)=.false.
     END DO
  
     gg=0
  
     DO ii=1,17
        IF (filtro(ii).eqv..true.) THEN
           bb=aux_mss(ii)
           gg=gg+1
           DO ll=1,17
              IF (filtro(ll).eqv..true.) THEN
                 IF (aux_mss(ll)==bb) THEN
                    aux_mss2(ll)=ii
                    filtro(ll)=.false.
                 END IF
              END IF
           END DO
        END IF
     END DO
  
     mss(jj,:)=aux_mss2(:)

  END DO


END SUBROUTINE remove_index_mol

SUBROUTINE remove_permutations_limit_4_nosim(mss,mss_ind,num_mss)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::num_mss
  INTEGER,DIMENSION(num_mss,17),INTENT(INOUT)::mss,mss_ind

  INTEGER::i,j,g,h,ii,kk,vigilo
  INTEGER,DIMENSION(17)::microstate,ms_short2,bb,key,key_aux
  LOGICAL,DIMENSION(17)::filtro
  INTEGER,DIMENSION(4,2)::aux_permut
  LOGICAL::interruptor
  INTEGER::ceros4,ceros5,ceros,x_ellos
  INTEGER::x_primera4,x_primera5,x_primera
  INTEGER::x_core,x_core4,x_core5
  INTEGER::x_segunda4,x_segunda5,x_segunda

  DO kk=1,num_mss

     key=mss_ind(kk,:)
     key_aux=key
     ms_short2=mss(kk,:)
     microstate=ms_short2
     aux_permut(1,1)=6
     aux_permut(1,2)=7
     aux_permut(2,1)=9
     aux_permut(2,2)=10
     aux_permut(3,1)=13
     aux_permut(3,2)=14
     aux_permut(4,1)=16
     aux_permut(4,2)=17
     
     DO h=1,2
        i=aux_permut(h,1)
        j=aux_permut(h,2)
        IF (ms_short2(i)>0) THEN
           IF (ms_short2(i)>ms_short2(j)) THEN
              microstate(i)=ms_short2(j)
              microstate(j)=ms_short2(i)
              key_aux(i)=key(j)
              key_aux(j)=key(i)
              filtro=.true.
              DO g=1,17
                 IF ((microstate(g)==i).and.(filtro(g).eqv..true.)) THEN
                    microstate(g)=j
                    filtro(g)=.false.
                 END IF
              END DO
              DO g=1,17
                 IF ((microstate(g)==j).and.(filtro(g).eqv..true.)) THEN
                    microstate(g)=i
                    filtro(g)=.false.
                 END IF
              END DO
              ms_short2=microstate
              key=key_aux
           END IF
        END IF
     END DO
     key=key_aux
     ms_short2=microstate
     
     
     
     ceros4=0
     ceros5=0
     ceros=0
     x_core4=0
     x_core5=0
     x_core=0
     x_primera4=0
     x_primera5=0
     x_primera=0
     x_segunda4=0
     x_segunda5=0
     x_segunda=0
     x_ellos=0
     
     DO i=12,14
        IF (ms_short2(i)==0) THEN
           ceros4=ceros4+1
        END IF
        IF (ms_short2(i)==1) THEN
           x_core4=x_core4+1
        END IF
        IF ((ms_short2(i)<=5).and.(ms_short2(i)>=2)) THEN
           x_primera4=x_primera4+1
        END IF
        IF ((ms_short2(i)<=11).and.(ms_short2(i)>=6)) THEN
           x_segunda4=x_segunda4+1
        END IF
        IF ((ms_short2(i)/=i).and.(ms_short2(i)>11)) THEN
           x_ellos=x_ellos+1
        END IF
     END DO
     DO i=15,17
        IF (ms_short2(i)==0) THEN
           ceros5=ceros5+1
        ELSE
           IF (ms_short2(i)==1) THEN
              x_core5=x_core5+1
           ELSE
              IF ((ms_short2(i)<=5).and.(ms_short2(i)>=2)) THEN
                 x_primera5=x_primera5+1
              ELSE
                 IF ((ms_short2(i)<=11).and.(ms_short2(i)>=6)) THEN
                    x_segunda5=x_segunda5+1
                 ELSE
                    IF ((ms_short2(i)/=i).and.(ms_short2(i)>11)) THEN
                       x_ellos=x_ellos+1
                    END IF
                 END IF
              END IF
           END IF
        END IF
     END DO
     ceros=ceros4+ceros5
     x_core=x_core4+x_core5
     x_primera=x_primera4+x_primera5
     x_segunda=x_segunda4+x_segunda5
     
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     interruptor=.false.
     
     IF (ms_short2(5)==0) THEN
        interruptor=.true.
     END IF
     IF (interruptor.eqv..false.) THEN
        IF (x_core/=0) THEN
           IF (x_core==1) THEN
              IF (x_core5==1) THEN
                 CALL DOY_VUELTA(ms_short2)
                 CALL DOY_VUELTA_KEY (key,key_aux)
                 microstate=ms_short2
                 interruptor=.true.
              END IF
           ELSE
              !print*,'problema x_core',frame,mol
              !print 117, ms_short2(:)
              !STOP
           END IF
        END IF
     END IF
     
     IF (interruptor.eqv..false.) THEN
        IF (ceros4/=ceros5) interruptor=.true.
        IF (ceros5<ceros4) THEN
           CALL DOY_VUELTA(ms_short2)
           CALL DOY_VUELTA_KEY (key,key_aux)
        END IF
     END IF
     IF (interruptor.eqv..false.) THEN
        IF (x_primera4/=x_primera5) interruptor=.true.
        IF (x_primera5<x_primera4) THEN
           CALL DOY_VUELTA(ms_short2)
           CALL DOY_VUELTA_KEY (key,key_aux)
        END IF
     END IF
     IF (interruptor.eqv..false.) THEN
        IF (x_segunda4/=x_segunda5) interruptor=.true.
        IF (x_segunda5<x_segunda4) THEN
           CALL DOY_VUELTA(ms_short2)
           CALL DOY_VUELTA_KEY (key,key_aux)
        END IF
     END IF
     
     !!
     microstate=ms_short2
     key_aux=key
     DO h=3,4
        i=aux_permut(h,1)
        j=aux_permut(h,2)
        IF (ms_short2(i)/=0) THEN
           IF (ms_short2(i)>ms_short2(j)) THEN
              key_aux(i)=key(j)
              key_aux(j)=key(i)
              microstate(i)=ms_short2(j)
              microstate(j)=ms_short2(i)
              filtro=.true.
              DO g=1,17
                 IF ((microstate(g)==i).and.(filtro(g).eqv..true.)) THEN
                    microstate(g)=j
                    filtro(g)=.false.
                 END IF
              END DO
              DO g=1,17
                 IF ((microstate(g)==j).and.(filtro(g).eqv..true.)) THEN
                    microstate(g)=i
                    filtro(g)=.false.
                 END IF
              END DO
              ms_short2=microstate
              key=key_aux
           END IF
        END IF
     END DO
     ms_short2=microstate
     key=key_aux
     !!
     
     
     !sigo
     
     IF (interruptor.eqv..false.) THEN
        IF (ceros>0) THEN
           IF ((ms_short2(12)==0).and.(ms_short2(15)/=0)) THEN
              interruptor=.true.
           ELSE
              IF ((ms_short2(12)/=0).and.(ms_short2(15)==0)) THEN
                 CALL DOY_VUELTA(ms_short2)
                 CALL DOY_VUELTA_KEY (key,key_aux)
                 interruptor=.true.
              END IF
           END IF
           IF (interruptor.eqv..false.) THEN
              IF ((ms_short2(13)==0).and.(ms_short2(16)/=0)) THEN
                 interruptor=.true.
              ELSE
                 IF ((ms_short2(13)/=0).and.(ms_short2(16)==0)) THEN
                    CALL DOY_VUELTA(ms_short2)
                    CALL DOY_VUELTA_KEY (key,key_aux)
                    interruptor=.true.
                 END IF
              END IF
           END IF
           IF (interruptor.eqv..false.) THEN
              IF ((ms_short2(14)==0).and.(ms_short2(17)/=0)) THEN
                 interruptor=.true.
              ELSE
                 IF ((ms_short2(14)/=0).and.(ms_short2(17)==0)) THEN
                    CALL DOY_VUELTA(ms_short2)
                    CALL DOY_VUELTA_KEY (key,key_aux)
                    interruptor=.true.
                 END IF
              END IF
           END IF
        END IF
     END IF
     
     
     IF (interruptor.eqv..false.) THEN
        IF ((ms_short2(12)==12).and.(ms_short2(15)/=15)) THEN
           interruptor=.true.
        ELSE
           IF ((ms_short2(12)/=12).and.(ms_short2(15)==15)) THEN
              CALL DOY_VUELTA(ms_short2)
              CALL DOY_VUELTA_KEY (key,key_aux)
              interruptor=.true.
           END IF
        END IF
        IF (interruptor.eqv..false.) THEN
           IF ((ms_short2(13)==13).and.(ms_short2(16)/=16)) THEN
              interruptor=.true.
           ELSE
              IF ((ms_short2(13)/=13).and.(ms_short2(16)==16)) THEN
                 CALL DOY_VUELTA(ms_short2)
                 CALL DOY_VUELTA_KEY (key,key_aux)
                 interruptor=.true.
              END IF
           END IF
        END IF
        IF (interruptor.eqv..false.) THEN
           IF ((ms_short2(14)==14).and.(ms_short2(17)/=17)) THEN
              interruptor=.true.
           ELSE
              IF ((ms_short2(14)/=14).and.(ms_short2(17)==17)) THEN
                 CALL DOY_VUELTA(ms_short2)
                 CALL DOY_VUELTA_KEY (key,key_aux)
                 interruptor=.true.
              END IF
           END IF
        END IF
     END IF
     microstate=ms_short2
     key_aux=key
     
     !bb=0
     !! Cruces entre ellos:
     ! - Un cruce
     ! bb(12:17)=(/ 12,13,14,14,16,17 /) !SII Elijo este como representante
     ! bb(12:17)=(/ 12,13,14,13,16,17 /) !SII ---> 12,13,14,14,16,17
     ! bb(12:17)=(/ 12,13,14,15,16,12 /) !NO
     !bb(12:17)=(/ 12,13,14,15,12,17 /) !SII
     ! - Dos cruces
     ! bb(12:17)=(/ 12,13,14,14,16,12 /) !NO
     ! bb(12:17)=(/ 12,13,14,14,12,17 /) !SII Elijo este como representante
     ! bb(12:17)=(/ 12,13,14,13,16,12 /) !NO
     ! bb(12:17)=(/ 12,13,14,13,12,17 /) !SII ---> 12,13,14,14,12,17
     ! - Tres cruces da igual
     !! Cruces entre ellos con la primera capa tambien:
     ! - Un cruce entre ellos y uno con primera capa
     ! bb(12:17)=(/ 5,13,14,14,4,17 /)
     ! bb(12:17)=(/ 5,13,14,14,16,4 /)
     !
     ! Corrijo esto:
     
     IF (ms_short2(16)==12) THEN
        CALL DOY_VUELTA(ms_short2)
        CALL DOY_VUELTA_KEY (key,key_aux)
     END IF
     microstate=ms_short2
     key_aux=key
     IF (ms_short2(15)==13) THEN
        key(13)=key_aux(14)
        key(14)=key_aux(13)
        key_aux(13:14)=key(13:14)
        ms_short2(13)=microstate(14)
        ms_short2(14)=microstate(13)
        microstate(13:14)=ms_short2(13:14)   !!! Esto lo acabo de cambiar del de Roman
        filtro=.true.
        DO i=12,17
           j=ms_short2(i)
           IF ((j>11).and.(filtro(i).eqv..true.)) THEN
              ms_short2(i)=i
              filtro(i)=.false.
              DO ii=i+1,17
                 IF ((microstate(ii)==j).and.(filtro(ii).eqv..true.)) THEN
                    ms_short2(ii)=i
                    filtro(ii)=.false.
                 END IF
              END DO
           END IF
        END DO
     END IF
     microstate=ms_short2
     key_aux=key
     
     
     IF (interruptor.eqv..false.) THEN
        IF ((ms_short2(13)==3).and.(ms_short2(16)==2)) THEN
           CALL DOY_VUELTA(ms_short2)
           CALL DOY_VUELTA_KEY (key,key_aux)
           microstate=ms_short2
           key_aux=key
        END IF
     END IF
     
     
     mss_ind(kk,:)=key(:)
     mss(kk,:)=microstate(:)

  END DO

!!!!!!!
     
     !IF (0) THEN
     ! interruptor=.true.
     ! DO i=12,17
     ! IF (ms_short2(i)/=bb(i)) THEN
     ! interruptor=.false.
     ! EXIT
     ! END IF
     ! END DO
     ! IF (interruptor.eqv..true.) THEN
     ! print*,'...>',frame,mol
     ! print 117,ms_short2(:)
     ! END IF
     !END IF
     
     
     !IF (interruptor.eqv..false.) THEN
     ! print*,'TATE',frame,mol
     ! print 117,ms_short2
     !END IF
     
     
     !117 format (1I3," |",4I3," |",3I3," |",3I3," |",3I3," |",3I3)


END SUBROUTINE remove_permutations_limit_4_nosim


SUBROUTINE DOY_VUELTA (ms_short2)

  IMPLICIT NONE

  INTEGER,DIMENSION(17),INTENT(INOUT)::ms_short2
  INTEGER::i,j,ii
  INTEGER,DIMENSION(17)::microstate
  LOGICAL,DIMENSION(17)::filtro
  
  filtro=.true.
  microstate=ms_short2
  
  ms_short2(4)=microstate(5)
  ms_short2(5)=microstate(4)
  ms_short2(12:14)=microstate(15:17)
  ms_short2(15:17)=microstate(12:14)
  microstate=ms_short2
  
  DO i=1,17
     j=ms_short2(i)
     IF ((j>1).and.(filtro(i).eqv..true.)) THEN
        ms_short2(i)=i
        filtro(i)=.false.
        DO ii=i+1,17
           IF ((microstate(ii)==j).and.(filtro(ii).eqv..true.)) THEN
              ms_short2(ii)=i
              filtro(ii)=.false.
           END IF
        END DO
     END IF
  END DO

END SUBROUTINE DOY_VUELTA

SUBROUTINE DOY_VUELTA_KEY (key,key_aux)
  
  IMPLICIT NONE
  
  INTEGER,DIMENSION(17),INTENT(INOUT)::key,key_aux

  key_aux=key
  key(4)=key_aux(5)
  key(5)=key_aux(4)
  key(12:14)=key_aux(15:17)
  key(15:17)=key_aux(12:14)
  key_aux=key

END SUBROUTINE DOY_VUELTA_KEY


SUBROUTINE addbonds (tipo,mss,mss_ind,bonds,num_wats,num_bonds)


  IMPLICIT NONE

  INTEGER,INTENT(IN)::tipo,num_wats,num_bonds
  INTEGER,DIMENSION(num_wats,17),INTENT(INOUT)::mss
  INTEGER,DIMENSION(num_wats,17),INTENT(INOUT)::mss_ind
  INTEGER,DIMENSION(num_bonds),INTENT(IN)::bonds
  
  INTEGER::ii,jj,kk

  LOGICAL,DIMENSION(num_wats)::dentro
  LOGICAL,DIMENSION(17)::corrijo

  dentro=.FALSE.
  DO ii=1,num_bonds
     dentro(bonds(ii)+1)=.TRUE.
  END DO

  SELECT CASE (tipo)
     CASE (1)
        DO jj=1,num_wats
           corrijo=.FALSE.
           DO kk=1,5
              IF (dentro(mss_ind(jj,kk)).eqv..TRUE.) THEN
                 IF (kk==1) THEN
                    IF (mss_ind(jj,4)==0) THEN
                       corrijo(4)=.TRUE.
                    ELSE
                       IF (mss_ind(jj,5)==0) THEN
                          corrijo(5)=.TRUE.
                       END IF
                    END IF
                 ELSEIF (kk==2) THEN
                    IF (mss_ind(jj,8)==0) THEN
                       corrijo(8)=.TRUE.
                    END IF
                 ELSEIF (kk==3) THEN
                    IF (mss_ind(jj,11)==0) THEN
                       corrijo(11)=.TRUE.
                    END IF
                 ELSEIF (kk==4) THEN
                    IF (mss_ind(jj,13)==0) THEN
                       corrijo(13)=.TRUE.
                    ELSE
                       IF (mss_ind(jj,14)==0) THEN
                          corrijo(14)=.TRUE.
                       END IF
                    END IF
                 ELSEIF (kk==5) THEN
                    IF (mss_ind(jj,16)==0) THEN
                       corrijo(16)=.TRUE.
                    ELSE
                       IF (mss_ind(jj,17)==0) THEN
                          corrijo(17)=.TRUE.
                       END IF
                    END IF
                 END IF
              END IF
           END DO
           DO kk=1,17
              IF (corrijo(kk).eqv..TRUE.) THEN
                 mss_ind(jj,kk)=-1
                 mss(jj,kk)=-1
              END IF
           END DO
        END DO

     CASE (2)

        DO jj=1,num_wats
           DO kk=1,17
              ii=mss_ind(jj,kk)
              IF (ii>0) THEN
                 IF (dentro(ii).eqv..TRUE.) THEN
                    mss(jj,kk)=-mss(jj,kk)
                 END IF
              END IF
           END DO
        END DO


     CASE (3)

        DO jj=1,num_wats
           DO kk=1,5
              ii=mss_ind(jj,kk)
              IF (ii>0) THEN
                 IF (dentro(ii).eqv..TRUE.) THEN
                    mss(jj,kk)=-mss(jj,kk)
                 END IF
              END IF
           END DO
        END DO



     CASE DEFAULT

        PRINT*,'NO OPTION CHOOSEN'

     END SELECT

END SUBROUTINE addbonds


END MODULE GLOB
