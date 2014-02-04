!!!#################################
!!!####  COMMON VARIABLES AND
!!!####  SUBROUTINES TO UPLOAD THEM
!!!#################################

MODULE GLOB

!! TOPOLOGY

INTEGER::num_ats,num_nods,num_sets,num_supsets
INTEGER,DIMENSION(:),ALLOCATABLE::at2nod
INTEGER,DIMENSION(:),ALLOCATABLE::in_at_nod,num_ats_nod
INTEGER,DIMENSION(:),ALLOCATABLE::lev_ats,lev_nods,lev_sets,lev_supsets
INTEGER,DIMENSION(:),ALLOCATABLE::topes_supsets
LOGICAL,DIMENSION(:),ALLOCATABLE::filtro_supsets
INTEGER,DIMENSION(:),ALLOCATABLE::trad2py_nod,trad2py_at

INTEGER,DIMENSION(:),ALLOCATABLE::ats_symm,ats_symm_num_crits,ats_symm_in,ats_symm_length
INTEGER,DIMENSION(:),ALLOCATABLE::vect_aux

LOGICAL::superfiltro_sets

!! NETWORK

INTEGER::Total_num_hbs,Total_num_bs
INTEGER,DIMENSION(:),ALLOCATABLE::in_Hbs,in_Bs,Hbs,Bs
INTEGER,DIMENSION(:),ALLOCATABLE::num_Hbs_at,num_Bs_at
INTEGER,DIMENSION(:),ALLOCATABLE::num_Hbs_Bs_nod

!! MICROSTATES

INTEGER,DIMENSION(:),ALLOCATABLE::mss_ind_ats,mss_ind_nods,mss

CONTAINS

SUBROUTINE load_topol(xx_at2nod,&
     xx_trad2py_at,xx_trad2py_nod,&
     xx_symm_ats_start,xx_symm_ats_crits,xx_symm_ats,&
     xx_symm_nods,xx_symm_nods_num,&
     xx_symm_sets,xx_symm_sets_num,&
     xx_num_ats,xx_num_nods,&
     xx_symm_ats_dim,xx_symm_nods_dim,xx_symm_sets_dim)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::xx_num_ats,xx_num_nods,xx_symm_nods_num,xx_symm_sets_num
  INTEGER,INTENT(IN)::xx_symm_ats_dim,xx_symm_nods_dim,xx_symm_sets_dim
  INTEGER,DIMENSION(xx_num_ats),INTENT(IN)::xx_at2nod
  INTEGER,DIMENSION(xx_num_ats),INTENT(IN)::xx_trad2py_at
  INTEGER,DIMENSION(xx_num_nods),INTENT(IN)::xx_trad2py_nod
  INTEGER,DIMENSION(xx_num_nods),INTENT(IN)::xx_symm_ats_start
  INTEGER,DIMENSION(xx_num_nods),INTENT(IN)::xx_symm_ats_crits
  INTEGER,DIMENSION(xx_symm_ats_dim),INTENT(IN)::xx_symm_ats
  INTEGER,DIMENSION(xx_symm_nods_dim),INTENT(IN)::xx_symm_nods
  INTEGER,DIMENSION(xx_symm_sets_dim),INTENT(IN)::xx_symm_sets

  INTEGER::gg,hh,ii,jj,kk,ll,mm,nn,oo

  IF (ALLOCATED(at2nod))           DEALLOCATE(at2nod)
  IF (ALLOCATED(trad2py_at))       DEALLOCATE(trad2py_at)
  IF (ALLOCATED(trad2py_nod))      DEALLOCATE(trad2py_nod)
  IF (ALLOCATED(in_at_nod))        DEALLOCATE(in_at_nod)
  IF (ALLOCATED(num_ats_nod))      DEALLOCATE(num_ats_nod)

  IF (ALLOCATED(ats_symm))           DEALLOCATE(ats_symm)
  IF (ALLOCATED(ats_symm_num_crits))  DEALLOCATE(ats_symm_num_crits)
  IF (ALLOCATED(ats_symm_in))        DEALLOCATE(ats_symm_in)
                                    
  IF (ALLOCATED(lev_ats))          DEALLOCATE(lev_ats)
  IF (ALLOCATED(lev_nods))         DEALLOCATE(lev_nods)
  IF (ALLOCATED(lev_sets))         DEALLOCATE(lev_sets)
  IF (ALLOCATED(lev_supsets))      DEALLOCATE(lev_supsets)
  IF (ALLOCATED(topes_supsets))    DEALLOCATE(topes_supsets)
  IF (ALLOCATED(filtro_supsets))   DEALLOCATE(filtro_supsets)

  num_ats     = xx_num_ats
  num_nods    = xx_num_nods
  num_sets    = xx_num_nods
  num_supsets = xx_num_nods

  ALLOCATE(at2nod(num_ats))
  ALLOCATE(trad2py_at(num_ats),trad2py_nod(num_nods))

  at2nod(:)      = xx_at2nod(:)
  trad2py_at(:)  = xx_trad2py_at(:)
  trad2py_nod(:) = xx_trad2py_nod(:)

  ALLOCATE(num_ats_nod(num_nods),in_at_nod(num_nods))
  ALLOCATE(lev_nods(num_nods),lev_sets(num_nods),lev_supsets(num_nods))
  ALLOCATE(filtro_supsets(num_nods))


  num_ats_nod(:)=0
  DO ii=1,num_ats
     jj=at2nod(ii)
     num_ats_nod(jj)=num_ats_nod(jj)+1
  END DO

  gg=0
  DO ii=1,num_nods
     in_at_nod(ii)=gg
     gg=gg+num_ats_nod(ii)
     lev_nods(ii)=1
     lev_sets(ii)=1
     lev_supsets(ii)=ii
  END DO

  ! SYMM ATS

  ALLOCATE(ats_symm(xx_symm_ats_dim))
  ALLOCATE(ats_symm_length(num_nods),ats_symm_num_crits(num_nods),ats_symm_in(num_nods))

  ats_symm(:)=xx_symm_ats
  ats_symm_num_crits(:)=xx_symm_ats_crits(:)
  ats_symm_in(:)=xx_symm_ats_start(:)

  DO ii=1,num_nods
     jj=ats_symm_num_crits(ii)
     kk=ats_symm_in(ii)-1
     gg=0
     DO ll=1,jj
        gg=gg+1
        gg=gg+ats_symm(kk+gg)
     END DO
     ats_symm_length(ii)=gg
  END DO

  ! SYMM NODS

  gg=1
  DO ii=1,xx_symm_nods_num
     hh=xx_symm_nods(gg)
     kk=xx_symm_nods(gg+1)
     kk=lev_supsets(kk)
     DO jj=gg+1,gg+hh
        ll=xx_symm_nods(jj)
        lev_supsets(ll)=kk
     END DO
     gg=gg+hh+1
  END DO
        
  ! SYMM SETS

  superfiltro_sets=.FALSE.
  filtro_supsets(:)=.FALSE.
  gg=1
  DO ii=1,xx_symm_sets_num
     superfiltro_sets=.TRUE.
     jj=xx_symm_sets(gg)
     kk=xx_symm_sets(gg+1)
     ll=xx_symm_sets(gg+2)
     ll=lev_supsets(ll)
     gg=gg+2
     DO mm=1,jj
        DO hh=1,kk
           nn=xx_symm_sets(gg)
           lev_supsets(nn)=ll
           lev_sets(nn)=mm
           lev_nods(nn)=hh
           filtro_supsets(nn)=.TRUE.
           gg=gg+1
        END DO
     END DO
  END DO

  ! TOPES SUPSETS

  jj=1
  DO ii=1,num_nods-1
     IF (lev_supsets(ii)/=(lev_supsets(ii+1))) jj=jj+1
  END DO

  ALLOCATE(topes_supsets(jj))
  jj=1
  topes_supsets(1)=jj
  DO ii=1,num_nods-1
     IF (lev_supsets(ii)/=(lev_supsets(ii+1))) THEN
        lev_supsets(ii)=jj
        jj=jj+1
        topes_supsets(jj)=lev_supsets(ii+1)
     ELSE
        lev_supsets(ii)=jj
     END IF
  END DO
  lev_supsets(num_nods)=jj



END SUBROUTINE load_topol


SUBROUTINE load_net(xx_hbs,xx_bs,xx_num_Hbs_at,xx_num_Bs_at,&
     xx_Total_num_hbs,xx_Total_num_bs,xx_num_ats)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::xx_num_ats,xx_Total_num_hbs,xx_Total_num_bs
  INTEGER,DIMENSION(xx_num_ats),INTENT(IN)::xx_num_hbs_at,xx_num_bs_at
  INTEGER,DIMENSION(xx_Total_num_hbs),INTENT(IN)::xx_hbs
  INTEGER,DIMENSION(xx_Total_num_bs),INTENT(IN)::xx_bs

  INTEGER::ii,jj,kk

  IF (ALLOCATED(in_hbs))       DEALLOCATE(in_hbs)
  IF (ALLOCATED(num_hbs_at))   DEALLOCATE(num_hbs_at)
  IF (ALLOCATED(hbs))          DEALLOCATE(hbs)

  IF (ALLOCATED(in_bs))        DEALLOCATE(in_bs)
  IF (ALLOCATED(num_bs_at))    DEALLOCATE(num_bs_at)
  IF (ALLOCATED(bs))           DEALLOCATE(bs)

  IF (ALLOCATED(num_Hbs_Bs_nod))   DEALLOCATE(num_Hbs_Bs_nod)

  Total_num_hbs = xx_Total_num_hbs
  Total_num_bs  = xx_Total_num_bs

  ALLOCATE(in_hbs(xx_num_ats),in_bs(xx_num_ats))
  ALLOCATE(hbs(Total_num_hbs),bs(Total_num_bs))
  ALLOCATE(num_hbs_at(xx_num_ats),num_bs_at(xx_num_ats))

  num_hbs_at(:) = xx_num_hbs_at(:)
  num_bs_at(:)  = xx_num_bs_at(:)
  hbs(:)        = xx_hbs(:)
  bs(:)         = xx_bs(:)

  jj=0
  kk=0
  DO ii=1,num_ats
     in_hbs(ii) = jj
     in_bs(ii)  = kk
     jj         = jj+num_hbs_at(ii)
     kk         = kk+num_bs_at(ii)
  END DO

  ALLOCATE(num_Hbs_Bs_nod(num_nods))

  num_Hbs_Bs_nod(:)=0
  DO ii=1,num_nods
     kk=0
     DO jj=in_at_nod(ii)+1,in_at_nod(ii)+num_ats_nod(ii)
        kk=kk+num_hbs_at(jj)+num_bs_at(jj)
     END DO
     num_Hbs_Bs_nod(ii)=kk
  END DO


END SUBROUTINE load_net


!!!#################################
!!!####  GENERAL FUNCTIONS TO SORT
!!!#################################

SUBROUTINE SORT_INT_1SH (nats,dim_aux,order,valores,num_crits)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::nats
  INTEGER,INTENT(INOUT)::dim_aux
  INTEGER,DIMENSION(nats),INTENT(INOUT)::order
  INTEGER,DIMENSION(dim_aux),INTENT(IN)::valores
  INTEGER,INTENT(INOUT)::num_crits

  INTEGER::ii,jj,kk,ll,gg,idim,new_num_crits,tope
  INTEGER,DIMENSION(:),ALLOCATABLE::val_aux,ind_aux,order_aux
  INTEGER,DIMENSION(:),ALLOCATABLE::vals,inds
  INTEGER,DIMENSION(:),ALLOCATABLE::new_symm_aux,cajon
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro
  LOGICAL::interruptor

  new_num_crits=0

  gg=0
  DO ii=1,num_crits
     gg=gg+1
     idim=vect_aux(gg)
     ALLOCATE(val_aux(idim),ind_aux(idim),order_aux(idim))
     DO jj=1,idim
        gg=gg+1
        kk=vect_aux(gg)
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
           IF (interruptor.eqv..FALSE.) THEN
              order_aux(1)=inds(jj-1)
              order_aux(2)=inds(jj)
              ll=2
              interruptor=.TRUE.
           ELSE
              ll=ll+1
              order_aux(ll)=inds(jj)
           END IF
        ELSE
           IF (interruptor.eqv..True.) THEN
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
        END IF
     END DO
     IF (interruptor.eqv..True.) THEN
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
  DEALLOCATE(vect_aux)
  IF (num_crits>0) THEN
     ALLOCATE(vect_aux(tope))
     dim_aux=tope
     vect_aux(:)=new_symm_aux(:)
     DEALLOCATE(new_symm_aux)
  ELSE
     dim_aux=0
  END IF

END SUBROUTINE SORT_INT_1SH

!!!#################################
!!!####  BUILD FIRST SHELL
!!!#################################

SUBROUTINE build_shell1st (core)
 
  IMPLICIT NONE
 
  TYPE iarray_pointer
     INTEGER,DIMENSION(:),POINTER::p1
  END TYPE iarray_pointer


  INTEGER,INTENT(IN)::core
 
  INTEGER::nats,nats2,nnods,ntot
  INTEGER::aa,bb,cc,dd
  INTEGER::ii,jj,kk

  INTEGER,DIMENSION(:),ALLOCATABLE::order_core,eff_order_core,vaux
  INTEGER,DIMENSION(:,:),ALLOCATABLE::dim_tree_core
  TYPE(iarray_pointer),DIMENSION(:,:),POINTER::tree_core
  INTEGER::num_crits

  IF (ALLOCATED(mss_ind_ats))   DEALLOCATE(mss_ind_ats)
  IF (ALLOCATED(mss_ind_nods))  DEALLOCATE(mss_ind_nods)
  IF (ALLOCATED(mss))           DEALLOCATE(mss)
 
  !! Construyo Arbol

  nats=num_ats_nod(core)
  nats2=2*nats
  nnods=num_Hbs_Bs_nod(core)

  ALLOCATE(order_core(nats),eff_order_core(nats))

  ii=in_at_nod(core)
  eff_order_core=(/(jj,jj=1,nats)/)
  order_core(:)=eff_order_core(:)+ii

  ALLOCATE(tree_core(nats,2),dim_tree_core(nats,2))

  DO ii=1,nats
     jj=order_core(ii)
     aa=num_Hbs_at(jj)
     bb=num_Bs_at(jj)
     cc=in_Hbs(jj)
     dd=in_Bs(jj)
     ALLOCATE(tree_core(ii,1)%p1(aa),tree_core(ii,2)%p1(bb))
     dim_tree_core(ii,:)=(/aa,bb/)
     tree_core(ii,1)%p1(:)=Hbs(cc+1:cc+aa)
     tree_core(ii,2)%p1(:)= Bs(dd+1:dd+bb)
  END DO

  !! Quito Simetrias

  DO ii=1,nats
     DO jj=1,2
        IF (dim_tree_core(ii,jj)>1) THEN
           print*,trad2py_nod(at2nod(tree_core(ii,jj)%p1(:)))
           CALL build_order_bonded_1st(tree_core(ii,jj)%p1(:),dim_tree_core(ii,jj),num_crits)
           print*,trad2py_nod(at2nod(tree_core(ii,jj)%p1(:)))
           IF (num_crits>0) THEN
              DEALLOCATE(vect_aux)
           END IF
        END IF
     END DO
  END DO


  ntot=1+nats+nats2+nnods
  ALLOCATE(mss_ind_ats(ntot),mss_ind_nods(ntot),mss(ntot))
  mss_ind_ats(:)=0
  mss_ind_nods(:)=0
  mss(:)=0

  jj=1
  mss_ind_ats(jj)  = nats
  mss_ind_nods(jj) = nats
  ii=jj+1
  jj=jj+nats
  mss_ind_ats(ii:jj)  = order_core(:)
  mss_ind_nods(ii:jj) = at2nod(mss_ind_ats(ii:jj))
  mss_ind_ats(ii:jj)  = trad2py_at(mss_ind_ats(ii:jj))
  mss_ind_nods(ii:jj) = trad2py_nod(mss_ind_nods(ii:jj))
  ii=jj+1
  jj=jj+nats2
  mss_ind_ats(ii:jj)  = (/(dim_tree_core(eff_order_core(kk),:),kk=1,nats)/)
  mss_ind_nods(ii:jj) = mss_ind_ats(ii:jj)
  ii=jj+1
  jj=jj+nnods
  mss_ind_ats(ii:jj) = (/((/tree_core(eff_order_core(kk),1)%p1(:),tree_core(eff_order_core(kk),2)%p1(:)/),kk=1,nats)/)
  mss_ind_nods(ii:jj) = at2nod(mss_ind_ats(ii:jj))
  mss_ind_ats(ii:jj)  = trad2py_at(mss_ind_ats(ii:jj))
  mss_ind_nods(ii:jj) = trad2py_nod(mss_ind_nods(ii:jj))

END SUBROUTINE build_shell1st

SUBROUTINE order_ats_1st (core,nats,order)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::core,nats
  INTEGER,DIMENSION(nats),INTENT(OUT)::order

  INTEGER::ii,jj,num_crits,dim_aux

  jj=in_at_nod(core)
  order(:)=(/(ii,ii=jj+1,jj+nats)/)

  num_crits=ats_symm_num_crits(core)

  IF (num_crits>0) THEN

     dim_aux=ats_symm_length(core)
     ii=ats_symm_in(core)
     ALLOCATE(vect_aux(dim_aux))
     vect_aux(:)=ats_symm(ii:(ii+dim_aux-1))
     CALL SORTBYNUMHBS_ATS_1SH(core,nats,order,num_crits,dim_aux)
     IF (num_crits>0) CALL SORTBYNUMBS_ATS_1SH(core,nats,order,num_crits,dim_aux)
     IF (num_crits>0) DEALLOCATE(vect_aux)

  END IF


END SUBROUTINE order_ats_1st


SUBROUTINE build_order_bonded_1st (order,dim_ord,num_crits)

  INTEGER::dim_ord
  INTEGER,DIMENSION(:),INTENT(INOUT)::order
  INTEGER,INTENT(OUT)::num_crits

  INTEGER::ii,dim_aux
  INTEGER,DIMENSION(dim_ord)::lista_nodos
  INTEGER,DIMENSION(:),ALLOCATABLE::valores
  
  lista_nodos(:)=at2nod(order(:))
  
  ALLOCATE(vect_aux(dim_ord+1),valores(dim_ord+1))
  vect_aux(:)=(/dim_ord,(/(ii,ii=1,dim_ord)/)/)
  valores(:)=(/0,(/lev_supsets(lista_nodos(:))/)/)

  num_crits=1
  dim_aux=dim_ord+1

  CALL SORT_INT_1SH(dim_ord,dim_aux,order,valores,num_crits)

  IF (superfiltro_sets.eqv..TRUE.) THEN
     IF (ANY(filtro_supsets(lista_nodos(:)),DIM=1).eqv..TRUE.) THEN

     END IF
  END IF

  DEALLOCATE(valores)

END SUBROUTINE build_order_bonded_1st

!SUBROUTINE build_shell1st_2 (core)
! 
!  IMPLICIT NONE
! 
!  INTEGER,INTENT(IN)::core
! 
!  INTEGER::nats,nats2,nnods,ntot
!  INTEGER::ff,gg,hh,ii,jj,kk,lhbs,lbs,mhbs,mbs
!  INTEGER,DIMENSION(:),ALLOCATABLE::order_ats_1sh,num_ats_bonded_1sh,order_ats_bonded_1sh
!  INTEGER,DIMENSION(:),ALLOCATABLE::order_nods_1sh,order_nods_bonded_1sh
! 
!  IF (ALLOCATED(mss_ind_ats)) DEALLOCATE(mss_ind_ats)
!  IF (ALLOCATED(mss_ind_nods)) DEALLOCATE(mss_ind_nods)
!  IF (ALLOCATED(mss))           DEALLOCATE(mss)
! 
!  nats=num_ats_nod(core)
!  nats2=2*nats
!  nnods=num_Hbs_Bs_nod(core)
! 
!  ALLOCATE(order_ats_1sh(nats),num_ats_bonded_1sh(nats2),order_ats_bonded_1sh(nnods))
! 
!  CALL order_ats_1st(core,nats,order_ats_1sh)
! 
!  gg=0
!  ff=0
!  DO ii=1,nats
!     jj=order_ats_1sh(ii)
!     lhbs=num_hbs_at(jj)
!     lbs=num_bs_at(jj)
!     mhbs=in_hbs(jj)
!     mbs=in_bs(jj)
!     gg=gg+1
!     num_ats_bonded_1sh(gg)=lhbs
!     gg=gg+1
!     num_ats_bonded_1sh(gg)=lbs
!     DO kk=1,lhbs
!        hh=mhbs+kk
!        ff=ff+1
!        order_ats_bonded_1sh(ff)=hbs(hh)
!     END DO
!     DO kk=1,lbs
!        hh=mbs+kk
!        ff=ff+1
!        order_ats_bonded_1sh(ff)=bs(hh)
!     END DO
!  END DO
! 
! 
!  !Translate
! 
!  ALLOCATE(order_nods_1sh(nats),order_nods_bonded_1sh(nnods))
! 
!  order_nods_1sh(:)        =at2nod(order_ats_1sh(:))
!  order_nods_bonded_1sh(:) =at2nod(order_ats_bonded_1sh(:))
! 
!  order_ats_1sh(:)          =trad2py_at(order_ats_1sh(:))
!  order_ats_bonded_1sh(:)   =trad2py_at(order_ats_bonded_1sh(:))
!  order_nods_1sh(:)         =trad2py_nod(order_nods_1sh(:))
!  order_nods_bonded_1sh(:)  =trad2py_nod(order_nods_bonded_1sh(:))
! 
!  !Build
! 
!  ntot=1+nats+nats2+nnods
!  ALLOCATE(mss_ind_ats(ntot),mss_ind_nods(ntot),mss(ntot))
! 
!  mss_ind_ats(1)  = nats
!  mss_ind_nods(1) = nats
!  ii=2
!  jj=1+nats
!  mss_ind_ats(ii:jj)  = order_ats_1sh(:)
!  mss_ind_nods(ii:jj) = order_nods_1sh(:)
!  ii=jj+1
!  jj=jj+nats2
!  mss_ind_ats(ii:jj)  = num_ats_bonded_1sh(:)
!  mss_ind_nods(ii:jj) = num_ats_bonded_1sh(:)
!  ii=jj+1
!  jj=jj+nnods
!  mss_ind_ats(ii:jj)  = order_ats_bonded_1sh(:)
!  mss_ind_nods(ii:jj) = order_nods_bonded_1sh(:)
! 
!  DEALLOCATE(order_ats_1sh,num_ats_bonded_1sh,order_ats_bonded_1sh)
!  DEALLOCATE(order_nods_1sh,order_nods_bonded_1sh)
! 
!END SUBROUTINE build_shell1st_2


SUBROUTINE SORTBYNUMHBS_ATS_1SH (core,nats,order,num_crits,dim_aux)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::core,nats
  INTEGER,INTENT(INOUT)::dim_aux
  INTEGER,DIMENSION(nats),INTENT(INOUT)::order
  INTEGER,INTENT(INOUT)::num_crits

  INTEGER::ii,jj,gg,kk,ll,idim
  INTEGER,DIMENSION(:),ALLOCATABLE::valores

  ALLOCATE(valores(dim_aux))

  gg=0
  DO ii=1,num_crits
     gg=gg+1
     idim=vect_aux(gg)
     DO jj=1,idim
        gg=gg+1
        kk=vect_aux(gg)
        ll=order(kk)
        valores(gg)=num_hbs_at(ll)
     END DO
  END DO

  CALL SORT_INT_1SH(nats,dim_aux,order,valores,num_crits)

  DEALLOCATE(valores)

END SUBROUTINE SORTBYNUMHBS_ATS_1SH

SUBROUTINE SORTBYNUMBS_ATS_1SH (core,nats,order,num_crits,dim_aux)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::core,nats
  INTEGER,INTENT(INOUT)::dim_aux
  INTEGER,DIMENSION(nats),INTENT(INOUT)::order
  INTEGER,INTENT(INOUT)::num_crits

  INTEGER::ii,jj,gg,kk,ll,idim
  INTEGER,DIMENSION(:),ALLOCATABLE::valores

  ALLOCATE(valores(dim_aux))

  gg=0
  DO ii=1,num_crits
     gg=gg+1
     idim=vect_aux(gg)
     DO jj=1,idim
        gg=gg+1
        kk=vect_aux(gg)
        ll=order(kk)
        valores(gg)=num_bs_at(ll)
     END DO
  END DO

  CALL SORT_INT_1SH(nats,dim_aux,order,valores,num_crits)

  DEALLOCATE(valores)

END SUBROUTINE SORTBYNUMBS_ATS_1SH





END MODULE GLOB


