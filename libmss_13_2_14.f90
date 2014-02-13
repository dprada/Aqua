!!!#################################
!!!####  COMMON VARIABLES AND
!!!####  SUBROUTINES TO UPLOAD THEM
!!!#################################

MODULE GLOB

!! TOPOLOGY

INTEGER::num_ats,num_nods,num_sets,num_supsets
INTEGER,DIMENSION(:),ALLOCATABLE::at2nod
INTEGER,DIMENSION(:),ALLOCATABLE::in_at_nod,num_ats_nod,num_ats_nod2
INTEGER,DIMENSION(:),ALLOCATABLE::lev_ats,lev_nods,lev_sets,lev_supsets
INTEGER,DIMENSION(:),ALLOCATABLE::topes_supsets
LOGICAL,DIMENSION(:),ALLOCATABLE::filtro_supsets
INTEGER,DIMENSION(:),ALLOCATABLE::trad2py_nod,trad2py_at

INTEGER,DIMENSION(:),ALLOCATABLE::ats_symm,ats_symm_num_crits,ats_symm_in,ats_symm_length

INTEGER,DIMENSION(:),ALLOCATABLE::symm_crit
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
  IF (ALLOCATED(num_ats_nod2))     DEALLOCATE(num_ats_nod2)

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

  ALLOCATE(num_ats_nod(num_nods),num_ats_nod2(num_nods),in_at_nod(num_nods))
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
     num_ats_nod2(ii)=num_ats_nod(ii)*2
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

!!$SUBROUTINE SORT_INT_1SH (nats,dim_aux,order,valores,num_crits)
!!$
!!$  IMPLICIT NONE
!!$
!!$  INTEGER,INTENT(IN)::nats
!!$  INTEGER,INTENT(INOUT)::dim_aux
!!$  INTEGER,DIMENSION(nats),INTENT(INOUT)::order
!!$  INTEGER,DIMENSION(dim_aux),INTENT(IN)::valores
!!$  INTEGER,INTENT(INOUT)::num_crits
!!$
!!$  INTEGER::ii,jj,kk,ll,gg,idim,new_num_crits,tope
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::val_aux,ind_aux,order_aux
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::vals,inds
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::new_symm_aux,cajon
!!$  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro
!!$  LOGICAL::interruptor
!!$
!!$  new_num_crits=0
!!$
!!$  gg=0
!!$  DO ii=1,num_crits
!!$     gg=gg+1
!!$     idim=vect_aux(gg)
!!$     ALLOCATE(val_aux(idim),ind_aux(idim),order_aux(idim))
!!$     DO jj=1,idim
!!$        gg=gg+1
!!$        kk=vect_aux(gg)
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
!!$  DEALLOCATE(vect_aux)
!!$  IF (num_crits>0) THEN
!!$     ALLOCATE(vect_aux(tope))
!!$     dim_aux=tope
!!$     vect_aux(:)=new_symm_aux(:)
!!$     DEALLOCATE(new_symm_aux)
!!$  ELSE
!!$     dim_aux=0
!!$  END IF
!!$
!!$END SUBROUTINE SORT_INT_1SH


SUBROUTINE SORT_INT_MATRIX_ATS (num_ats,dim_vecti_aux,dim_matrix,order,valores,num_crits)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::num_ats,dim_matrix
  INTEGER,INTENT(INOUT)::dim_vecti_aux
  INTEGER,DIMENSION(num_ats),INTENT(INOUT)::order
  INTEGER,DIMENSION(dim_vecti_aux,dim_matrix),INTENT(IN)::valores
  INTEGER,INTENT(INOUT)::num_crits

  INTEGER::ii,jj,kk,ll,gg,nn,pp,idim,new_num_crits,tope
  INTEGER,DIMENSION(:),ALLOCATABLE::val_aux,ind_aux,order_aux
  INTEGER,DIMENSION(:),ALLOCATABLE::vals,inds
  INTEGER,DIMENSION(:),ALLOCATABLE::new_symm_aux,cajon
  INTEGER,DIMENSION(num_ats)::trad_vals
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro
  LOGICAL::interruptor

  gg=0
  DO ii=1,num_crits
     gg=gg+1
     DO jj=1,symm_crit(gg)
        gg=gg+1
        trad_vals(symm_crit(gg))=gg
     END DO
  END DO

  pp=1
  DO WHILE ((num_crits>0).AND.(pp<=dim_matrix))

     new_num_crits=0
     
     gg=0
     DO ii=1,num_crits
        gg=gg+1
        idim=symm_crit(gg)
        ALLOCATE(val_aux(idim))
        DO jj=1,idim
           val_aux(jj)=valores(trad_vals(symm_crit(gg+jj)),pp)
        END DO
        IF (ALL(val_aux(2:).eq.val_aux(1)).eqv..TRUE.) THEN
           new_num_crits=new_num_crits+1
           IF (new_num_crits==1) THEN
              tope=idim+1
              ALLOCATE(new_symm_aux(tope))
              new_symm_aux(1)=idim
              new_symm_aux(2:tope)=symm_crit((gg+1):(gg+idim))
           ELSE
              ALLOCATE(cajon(tope))
              cajon(:)=new_symm_aux(:)
              DEALLOCATE(new_symm_aux)
              ALLOCATE(new_symm_aux(tope+1+idim))
              new_symm_aux(1:tope)=cajon(:)
              tope=tope+1
              new_symm_aux(tope)=idim
              new_symm_aux((tope+1):(tope+idim))=symm_crit((gg+1):(gg+idim))
              tope=tope+idim
              DEALLOCATE(cajon)
           END IF
           gg=gg+idim
        ELSE
           ALLOCATE(ind_aux(idim),order_aux(idim))
           DO jj=1,idim
              gg=gg+1
              kk=symm_crit(gg)
              ll=order(kk)
              ind_aux(jj)=kk
              order_aux(jj)=ll
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
           DEALLOCATE(ind_aux,order_aux)
           DEALLOCATE(filtro,vals,inds)
        END IF
        DEALLOCATE(val_aux)
     END DO
     
     num_crits=new_num_crits
     DEALLOCATE(symm_crit)
     IF (num_crits==0) THEN
        dim_vecti_aux=0
     ELSE
        ALLOCATE(symm_crit(tope))
        dim_vecti_aux=tope
        symm_crit(:)=new_symm_aux(:)
        DEALLOCATE(new_symm_aux)
     END IF

     pp=pp+1
  END DO



END SUBROUTINE SORT_INT_MATRIX_ATS



!!!#################################
!!!####  BUILD FIRST SHELL
!!!#################################


SUBROUTINE build_shell1st (core)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::core

  INTEGER::aa,bb,cc,dd,ee,ii,jj,kk,ll

  INTEGER::s1_nats,s1_nats2,s1_nnods,s1_ntot,s1_ats_symm_num_crits
  INTEGER,DIMENSION(:),ALLOCATABLE::s1_ats_symm,s1_order_ats,s1_num_bonds,s1_order_bonded_ats


  s1_nats=num_ats_nod(core)
  s1_nats2=num_ats_nod2(core)
  s1_nnods=num_Hbs_Bs_nod(core)
  s1_ntot=s1_nats+s1_nats2+s1_nnods+1
  ll=ats_symm_length(core)
  ii=ats_symm_in(core)
  jj=ii+ll
  ALLOCATE(s1_ats_symm(ll))
  s1_ats_symm_num_crits=ats_symm_num_crits(core)
  s1_ats_symm(:)=ats_symm(ii:jj)
  ALLOCATE(s1_order_ats(s1_nats))
  ii=in_at_nod(core)
  s1_order_ats(:)=(/(jj,jj=1,s1_nats)/)+ii
  ALLOCATE(s1_num_bonds(s1_nats2))
  ALLOCATE(s1_order_bonded_ats(s1_nnods))
  kk=1
  ll=0
  DO ii=1,s1_nats
     jj=s1_order_ats(ii)
     aa=num_Hbs_at(jj)
     bb=num_Bs_at(jj)
     cc=aa+bb
     dd=in_Hbs(jj)
     ee=in_Bs(jj)
     s1_num_bonds(kk)=aa
     s1_num_bonds(kk+1)=bb
     kk=kk+2
     s1_order_bonded_ats((ll+1):(ll+cc))=(/Hbs(dd+1:dd+aa),Bs(ee+1:ee+bb)/)
     ll=ll+cc
  END DO

  IF (ALLOCATED(mss_ind_ats))   DEALLOCATE(mss_ind_ats)
  IF (ALLOCATED(mss_ind_nods))  DEALLOCATE(mss_ind_nods)
  IF (ALLOCATED(mss))           DEALLOCATE(mss)

  ALLOCATE(mss_ind_ats(s1_ntot),mss_ind_nods(s1_ntot),mss(s1_ntot))
  mss_ind_ats(:)=0
  mss_ind_nods(:)=0
  mss(:)=0

  jj=1
  mss_ind_ats(jj)  = s1_nats
  mss_ind_nods(jj) = s1_nats
  ii=jj+1
  jj=jj+s1_nats
  mss_ind_ats(ii:jj)  = s1_order_ats(:)
  mss_ind_nods(ii:jj) = at2nod(s1_order_ats(:))
  mss_ind_ats(ii:jj)  = trad2py_at(s1_order_ats(:))
  mss_ind_nods(ii:jj) = trad2py_nod(mss_ind_nods(ii:jj))
  ii=jj+1
  jj=jj+s1_nats2
  mss_ind_ats(ii:jj)  = s1_num_bonds(:)
  mss_ind_nods(ii:jj) = s1_num_bonds(:)
  ii=jj+1
  jj=jj+s1_nnods
  mss_ind_ats(ii:jj)  = s1_order_bonded_ats(:)
  mss_ind_nods(ii:jj) = at2nod(s1_order_bonded_ats(:))
  mss_ind_ats(ii:jj)  = trad2py_at(s1_order_bonded_ats(:))
  mss_ind_nods(ii:jj) = trad2py_nod(mss_ind_nods(ii:jj))

  DEALLOCATE(s1_ats_symm,s1_order_ats,s1_num_bonds,s1_order_bonded_ats)

END SUBROUTINE build_shell1st


SUBROUTINE build_shell2nd (core)

  IMPLICIT NONE

  TYPE iarray_pointer
     INTEGER,DIMENSION(:),POINTER::p1
  END TYPE iarray_pointer

  INTEGER,INTENT(IN)::core

  INTEGER::aa,bb,cc,dd,ee,ii,jj,kk,ll,mm,nn,gg,iii,jjj,aaa,bbb,ccc,ggg
  INTEGER::core2

  INTEGER::s1_nats,s1_nats2,s1_nnods,s1_ntot,s1_ats_symm_num_crits
  INTEGER,DIMENSION(:),ALLOCATABLE::s1_ats_symm,s1_order_ats,s1_num_bonds,s1_order_bonded_ats

  INTEGER::si_nats,si_nats2,si_nnods
  INTEGER,DIMENSION(:),ALLOCATABLE::s2_nats,s2_nats2,s2_nnods,s2_ntot,s2_ats_symm_num_crits
  TYPE(iarray_pointer),DIMENSION(:),POINTER::s2_ats_symm,s2_order_ats,s2_num_bonds,s2_order_bonded_ats

  INTEGER::num_crits,dim_symm_crits,dim_mat
  INTEGER,DIMENSION(:),ALLOCATABLE::aux_order,aux_nod,aux_at
  INTEGER,DIMENSION(:,:),ALLOCATABLE::valores,box
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro,filtro2

  INTEGER::ntot

  s1_nats=num_ats_nod(core)
  s1_nats2=num_ats_nod2(core)
  s1_nnods=num_Hbs_Bs_nod(core)
  s1_ntot=s1_nats+s1_nats2+s1_nnods+1
  ll=ats_symm_length(core)
  ii=ats_symm_in(core)
  jj=ii+ll
  ALLOCATE(s1_ats_symm(ll))
  s1_ats_symm_num_crits=ats_symm_num_crits(core)
  s1_ats_symm(:)=ats_symm(ii:jj)
  ALLOCATE(s1_order_ats(s1_nats))
  ii=in_at_nod(core)
  s1_order_ats(:)=(/(jj,jj=1,s1_nats)/)+ii
  ALLOCATE(s1_num_bonds(s1_nats2))
  ALLOCATE(s1_order_bonded_ats(s1_nnods))
  kk=1
  ll=0
  DO ii=1,s1_nats
     jj=s1_order_ats(ii)
     aa=num_Hbs_at(jj)
     bb=num_Bs_at(jj)
     cc=aa+bb
     dd=in_Hbs(jj)
     ee=in_Bs(jj)
     s1_num_bonds(kk)=aa
     s1_num_bonds(kk+1)=bb
     kk=kk+2
     s1_order_bonded_ats((ll+1):(ll+cc))=(/Hbs(dd+1:dd+aa),Bs(ee+1:ee+bb)/)
     ll=ll+cc
  END DO

  ALLOCATE(s2_nats(s1_nnods),s2_nats2(s1_nnods),s2_nnods(s1_nnods),s2_ntot(s1_nnods),s2_ats_symm_num_crits(s1_nnods))
  ALLOCATE(s2_ats_symm(s1_nnods),s2_order_ats(s1_nnods),s2_num_bonds(s1_nnods),s2_order_bonded_ats(s1_nnods))

  DO iii=1,s1_nnods
     core2=s1_order_bonded_ats(iii)
     core2=at2nod(core2)
     si_nats=num_ats_nod(core2)
     si_nats2=num_ats_nod2(core2)
     si_nnods=num_Hbs_Bs_nod(core2)
     s2_nats(iii)=si_nats
     s2_nats2(iii)=si_nats2
     s2_nnods(iii)=si_nnods
     s2_ntot(iii)=si_nats+si_nats2+si_nnods+1
     ll=ats_symm_length(core)
     ii=ats_symm_in(core)
     jj=ii+ll
     ALLOCATE(s2_ats_symm(iii)%p1(ll))
     s2_ats_symm_num_crits(iii)=ats_symm_num_crits(core2)
     s2_ats_symm(iii)%p1(:)=ats_symm(ii:jj)
     ALLOCATE(s2_order_ats(iii)%p1(si_nats))
     ii=in_at_nod(core2)
     s2_order_ats(iii)%p1(:)=(/(jj,jj=1,si_nats)/)+ii
     ALLOCATE(s2_num_bonds(iii)%p1(si_nats2))
     ALLOCATE(s2_order_bonded_ats(iii)%p1(si_nnods))
     kk=1
     ll=0
     DO ii=1,si_nats
        jj=s2_order_ats(iii)%p1(ii)
        aa=num_Hbs_at(jj)
        bb=num_Bs_at(jj)
        cc=aa+bb
        dd=in_Hbs(jj)
        ee=in_Bs(jj)
        s2_num_bonds(iii)%p1(kk)=aa
        s2_num_bonds(iii)%p1(kk+1)=bb
        kk=kk+2
        s2_order_bonded_ats(iii)%p1((ll+1):(ll+cc))=(/Hbs(dd+1:dd+aa),Bs(ee+1:ee+bb)/)
        ll=ll+cc
     END DO
  END DO


  !!! Quito simetrias en 2nd shell order bonded

  DO iii=1,s1_nnods

     jj=0
     DO jjj=1,s2_nats2(iii)
        si_nats=s2_num_bonds(iii)%p1(jjj)
        ii=jj+1
        jj=jj+si_nats
        IF (si_nats>1) THEN
           ALLOCATE(aux_at(si_nats))
           aux_at(:)=s2_order_bonded_ats(iii)%p1(ii:jj)
           CALL bueno(core,aux_at,si_nats)
           s2_order_bonded_ats(iii)%p1(ii:jj)=aux_at(:)

           mm=ii
           DO nn=1,si_nats              
              s2_order_bonded_ats(iii)%p1(mm)=aux_at(aux_order(nn))
              mm=mm+1
           END DO
           DEALLOCATE(aux_order,aux_nod,aux_at,filtro)
           IF (ALLOCATED(symm_crit)) DEALLOCATE(symm_crit)
        END IF
     END DO

  END DO

  !!! Quito simetrias en 2nd shell order ats (Y si lo engancho con el hecho de no tener antes armado el microestado, puedo ahorrar calculo)
  !!! Deberia armar el microestado justo despues de este paso.

  !DO iii=1,s1_nnods


  !END DO

  ! Quito simetrias order bonded 1st

  ! Quito simetrias order ats 1st


  ntot=s1_ntot+SUM(s2_ntot(:),DIM=1)

  IF (ALLOCATED(mss_ind_ats))   DEALLOCATE(mss_ind_ats)
  IF (ALLOCATED(mss_ind_nods))  DEALLOCATE(mss_ind_nods)
  IF (ALLOCATED(mss))           DEALLOCATE(mss)

  ALLOCATE(mss_ind_ats(ntot),mss_ind_nods(ntot),mss(ntot))
  mss_ind_ats(:)=0
  mss_ind_nods(:)=0
  mss(:)=0

  jj=1
  mss_ind_ats(jj)  = s1_nats
  mss_ind_nods(jj) = s1_nats
  ii=jj+1
  jj=jj+s1_nats
  mss_ind_ats(ii:jj)  = s1_order_ats(:)
  mss_ind_nods(ii:jj) = at2nod(s1_order_ats(:))
  mss_ind_ats(ii:jj)  = trad2py_at(s1_order_ats(:))
  mss_ind_nods(ii:jj) = trad2py_nod(mss_ind_nods(ii:jj))
  ii=jj+1
  jj=jj+s1_nats2
  mss_ind_ats(ii:jj)  = s1_num_bonds(:)
  mss_ind_nods(ii:jj) = s1_num_bonds(:)
  ii=jj+1
  jj=jj+s1_nnods
  mss_ind_ats(ii:jj)  = s1_order_bonded_ats(:)
  mss_ind_nods(ii:jj) = at2nod(s1_order_bonded_ats(:))
  mss_ind_ats(ii:jj)  = trad2py_at(s1_order_bonded_ats(:))
  mss_ind_nods(ii:jj) = trad2py_nod(mss_ind_nods(ii:jj))

  DO iii=1,s1_nnods
     jj=jj+1
     si_nats=s2_nats(iii) 
     mss_ind_ats(jj)  = si_nats
     mss_ind_nods(jj) = si_nats
     ii=jj+1
     jj=jj+si_nats
     mss_ind_ats(ii:jj)  = s2_order_ats(iii)%p1(:)
     mss_ind_nods(ii:jj) = at2nod(mss_ind_ats(ii:jj))
     mss_ind_ats(ii:jj)  = trad2py_at(mss_ind_ats(ii:jj))
     mss_ind_nods(ii:jj) = trad2py_nod(mss_ind_nods(ii:jj))
     ii=jj+1
     jj=jj+s2_nats2(iii)
     mss_ind_ats(ii:jj)  = s2_num_bonds(iii)%p1(:)
     mss_ind_nods(ii:jj) = mss_ind_ats(ii:jj)
     ii=jj+1
     jj=jj+s2_nnods(iii)
     mss_ind_ats(ii:jj)  = s2_order_bonded_ats(iii)%p1(:)
     mss_ind_nods(ii:jj) = at2nod(mss_ind_ats(ii:jj))
     mss_ind_ats(ii:jj)  = trad2py_at(mss_ind_ats(ii:jj))
     mss_ind_nods(ii:jj) = trad2py_nod(mss_ind_nods(ii:jj))
     DEALLOCATE(s2_ats_symm(iii)%p1,s2_order_ats(iii)%p1,s2_num_bonds(iii)%p1,s2_order_bonded_ats(iii)%p1)
  END DO


  DEALLOCATE(s1_ats_symm,s1_order_ats,s1_num_bonds,s1_order_bonded_ats)
  DEALLOCATE(s2_nats,s2_nats2,s2_nnods,s2_ntot,s2_ats_symm_num_crits)
  DEALLOCATE(s2_ats_symm,s2_order_ats,s2_num_bonds,s2_order_bonded_ats)

END SUBROUTINE build_shell2nd

SUBROUTINE bueno()


  num_crits=1
  dim_symm_crits=si_nats+1
  ALLOCATE(aux_order(si_nats),aux_nod(si_nats),aux_at(si_nats),filtro(si_nats),symm_crit(dim_symm_crits))
  aux_order(:)=(/(aa,aa=1,si_nats)/)
  symm_crit(:)=(/si_nats,aux_order(:)/)
  aux_nod(:)=at2nod(aux_at(:))
  IF (ANY(aux_nod(:).eq.core).eqv..TRUE.) THEN ! SI EL CORE ESTA EN ORDER
     IF (superfiltro_sets.eqv..TRUE.) THEN !SI EXISTEN SUPSETS
        IF (filtro_supsets(core).eqv..TRUE.) THEN !SI EL CORE ESTA EN SUPSETS
           dim_mat=3
           ALLOCATE(valores(dim_symm_crits,dim_mat))
           valores(:,:)=0
           aa=lev_supsets(core)
           bb=lev_sets(core)
           cc=core
           gg=0
           DO nn=1,num_crits
              gg=gg+1
              DO mm=1,symm_crit(gg)
                 gg=gg+1
                 kk=aux_nod(symm_crit(gg))
                 aaa=lev_supsets(kk)
                 bbb=lev_sets(kk)
                 ccc=kk
                 IF (aa==aaa) THEN
                    valores(gg,1)=1
                    IF (bb==bbb) THEN
                       valores(gg,2)=1
                       IF (cc==ccc) THEN
                          valores(gg,3)=1
                       END IF
                    END IF
                 END IF
              END DO
           END DO
           CALL SORT_INT_MATRIX_ATS (si_nats,dim_symm_crits,dim_mat,aux_order,valores,num_crits)
           DEALLOCATE(valores)
        ELSE !SI EL CORE NO ESTA EN SUPSETS
           dim_mat=2
           ALLOCATE(valores(dim_symm_crits,dim_mat))
           valores(:,:)=0
           aa=lev_supsets(core)
           bb=core
           gg=0
           DO nn=1,num_crits
              gg=gg+1
              DO mm=1,symm_crit(gg)
                 gg=gg+1
                 kk=aux_nod(symm_crit(gg))
                 aaa=lev_supsets(kk)
                 bbb=kk
                 IF (aa==aaa) THEN
                    valores(gg,1)=1
                    IF (bb==bbb) THEN
                       valores(gg,2)=1
                    END IF
                 END IF
              END DO
           END DO
           CALL SORT_INT_MATRIX_ATS (si_nats,dim_symm_crits,dim_mat,aux_order,valores,num_crits)
           DEALLOCATE(valores)
        END IF
     ELSE ! SI NO EXISTEN SUPSETS
        dim_mat=2
        ALLOCATE(valores(dim_symm_crits,dim_mat))
        valores(:,:)=0
        aa=lev_supsets(core)
        bb=core
        gg=0
        DO nn=1,num_crits
           gg=gg+1
           DO mm=1,symm_crit(gg)
              gg=gg+1
              kk=aux_nod(symm_crit(gg))
              aaa=lev_supsets(kk)
              bbb=kk
              IF (aa==aaa) THEN
                 valores(gg,1)=1
                 IF (bb==bbb) THEN
                    valores(gg,2)=1
                 END IF
              END IF
           END DO
        END DO
        CALL SORT_INT_MATRIX_ATS (si_nats,dim_symm_crits,dim_mat,aux_order,valores,num_crits)
        DEALLOCATE(valores)
     END IF
  END IF
  !! aqui podria poner mas condiciones como si hay loops en 2nd shell
  IF (num_crits>0) THEN
     IF (superfiltro_sets.eqv..TRUE.) THEN
        filtro(:)=.FALSE.
        gg=0
        DO nn=1,num_crits
           gg=gg+1
           DO mm=1,symm_crit(gg)
              gg=gg+1
              filtro(symm_crit(gg))=.TRUE.
           END DO
        END DO
        IF (ANY(filtro_supsets(PACK(aux_nod(:),MASK=filtro))).eqv..TRUE.) THEN
           dim_mat=3
           ALLOCATE(valores(dim_symm_crits,dim_mat),filtro2(dim_symm_crits))
           filtro2=.TRUE.
           gg=0
           DO nn=1,num_crits
              gg=gg+1
              filtro2(gg)=.FALSE.
              DO mm=1,symm_crit(gg)
                 gg=gg+1
                 kk=aux_nod(symm_crit(gg))
                 valores(gg,1)=lev_supsets(kk)
                 valores(gg,2)=lev_sets(kk)
                 valores(gg,3)=lev_nods(kk)
              END DO
           END DO
           aa=MAXVAL(valores(:,1),DIM=1,MASK=filtro2)
           bb=MAXVAL(valores(:,2),DIM=1,MASK=filtro2)
           DEALLOCATE(filtro2)
           ALLOCATE(box(aa,bb))
           gg=0
           DO nn=1,num_crits
              gg=gg+1
              box(:,:)=0
              ggg=gg
              DO mm=1,symm_crit(gg)
                 gg=gg+1
                 aa=valores(gg,1)
                 bb=valores(gg,2)
                 box(aa,bb)=box(aa,bb)+1
              END DO
              gg=ggg
              DO mm=1,symm_crit(gg)
                 gg=gg+1
                 aa=valores(gg,1)
                 bb=valores(gg,2)
                 valores(gg,2)=box(aa,bb)
              END DO
           END DO
           CALL SORT_INT_MATRIX_ATS (si_nats,dim_symm_crits,dim_mat,aux_order,valores,num_crits)
           DEALLOCATE(valores,box)
        ELSE
           dim_mat=1
           ALLOCATE(valores(dim_symm_crits,dim_mat))
           gg=0
           DO nn=1,num_crits
              gg=gg+1
              DO mm=1,symm_crit(gg)
                 gg=gg+1
                 kk=aux_nod(symm_crit(gg))
                 valores(gg,1)=lev_supsets(kk)
              END DO
           END DO
           CALL SORT_INT_MATRIX_ATS (si_nats,dim_symm_crits,dim_mat,aux_order,valores,num_crits)
           DEALLOCATE(valores)
        END IF
     ELSE
        dim_mat=1
        ALLOCATE(valores(dim_symm_crits,dim_mat))
        gg=0
        DO nn=1,num_crits
           gg=gg+1
           DO mm=1,symm_crit(gg)
              gg=gg+1
              kk=aux_nod(symm_crit(gg))
              valores(gg,1)=lev_supsets(kk)
           END DO
        END DO
        CALL SORT_INT_MATRIX_ATS (si_nats,dim_symm_crits,dim_mat,aux_order,valores,num_crits)
        DEALLOCATE(valores)
     END IF
  END IF
  !aqui reordeno y recoloco

  aux_at(:)=aux_at(aux_order(nn))

  DEALLOCATE(aux_order,aux_nod,filtro)
  IF (ALLOCATED(symm_crit)) DEALLOCATE(symm_crit)

END SUBROUTINE bueno

!!$SUBROUTINE build_shell1st_old (core)
!!$ 
!!$  IMPLICIT NONE
!!$ 
!!$  TYPE iarray_pointer
!!$     INTEGER,DIMENSION(:),POINTER::p1
!!$  END TYPE iarray_pointer
!!$
!!$
!!$  INTEGER,INTENT(IN)::core
!!$ 
!!$  INTEGER::nats,nats2,nnods,ntot
!!$  INTEGER::aa,bb,cc,dd
!!$  INTEGER::ii,jj,kk
!!$
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::order_core,eff_order_core,vaux
!!$  INTEGER,DIMENSION(:,:),ALLOCATABLE::dim_tree_core
!!$  TYPE(iarray_pointer),DIMENSION(:,:),POINTER::tree_core
!!$  INTEGER::num_crits
!!$
!!$  IF (ALLOCATED(mss_ind_ats))   DEALLOCATE(mss_ind_ats)
!!$  IF (ALLOCATED(mss_ind_nods))  DEALLOCATE(mss_ind_nods)
!!$  IF (ALLOCATED(mss))           DEALLOCATE(mss)
!!$ 
!!$  !! Construyo Arbol
!!$
!!$  nats=num_ats_nod(core)
!!$  nats2=2*nats
!!$  nnods=num_Hbs_Bs_nod(core)
!!$
!!$  ALLOCATE(order_core(nats),eff_order_core(nats))
!!$
!!$  ii=in_at_nod(core)
!!$  eff_order_core=(/(jj,jj=1,nats)/)
!!$  order_core(:)=eff_order_core(:)+ii
!!$
!!$  ALLOCATE(tree_core(nats,2),dim_tree_core(nats,2))
!!$
!!$  DO ii=1,nats
!!$     jj=order_core(ii)
!!$     aa=num_Hbs_at(jj)
!!$     bb=num_Bs_at(jj)
!!$     cc=in_Hbs(jj)
!!$     dd=in_Bs(jj)
!!$     ALLOCATE(tree_core(ii,1)%p1(aa),tree_core(ii,2)%p1(bb))
!!$     dim_tree_core(ii,:)=(/aa,bb/)
!!$     tree_core(ii,1)%p1(:)=Hbs(cc+1:cc+aa)
!!$     tree_core(ii,2)%p1(:)= Bs(dd+1:dd+bb)
!!$  END DO
!!$
!!$  !! Quito Simetrias
!!$
!!$  DO ii=1,nats
!!$     DO jj=1,2
!!$        IF (dim_tree_core(ii,jj)>1) THEN
!!$           CALL build_order_bonded_1st(tree_core(ii,jj)%p1(:),dim_tree_core(ii,jj),num_crits)
!!$           IF (num_crits>0) THEN
!!$              DEALLOCATE(vect_aux)
!!$           END IF
!!$        END IF
!!$     END DO
!!$  END DO
!!$
!!$
!!$  ntot=1+nats+nats2+nnods
!!$  ALLOCATE(mss_ind_ats(ntot),mss_ind_nods(ntot),mss(ntot))
!!$  mss_ind_ats(:)=0
!!$  mss_ind_nods(:)=0
!!$  mss(:)=0
!!$
!!$  jj=1
!!$  mss_ind_ats(jj)  = nats
!!$  mss_ind_nods(jj) = nats
!!$  ii=jj+1
!!$  jj=jj+nats
!!$  mss_ind_ats(ii:jj)  = order_core(:)
!!$  mss_ind_nods(ii:jj) = at2nod(mss_ind_ats(ii:jj))
!!$  mss_ind_ats(ii:jj)  = trad2py_at(mss_ind_ats(ii:jj))
!!$  mss_ind_nods(ii:jj) = trad2py_nod(mss_ind_nods(ii:jj))
!!$  ii=jj+1
!!$  jj=jj+nats2
!!$  mss_ind_ats(ii:jj)  = (/(dim_tree_core(eff_order_core(kk),:),kk=1,nats)/)
!!$  mss_ind_nods(ii:jj) = mss_ind_ats(ii:jj)
!!$  ii=jj+1
!!$  jj=jj+nnods
!!$  mss_ind_ats(ii:jj) = (/((/tree_core(eff_order_core(kk),1)%p1(:),tree_core(eff_order_core(kk),2)%p1(:)/),kk=1,nats)/)
!!$  mss_ind_nods(ii:jj) = at2nod(mss_ind_ats(ii:jj))
!!$  mss_ind_ats(ii:jj)  = trad2py_at(mss_ind_ats(ii:jj))
!!$  mss_ind_nods(ii:jj) = trad2py_nod(mss_ind_nods(ii:jj))
!!$
!!$END SUBROUTINE build_shell1st_old
!!$
!!$
!!$SUBROUTINE build_shell2nd_old (core)
!!$ 
!!$  IMPLICIT NONE
!!$ 
!!$  TYPE iarray_pointer
!!$     INTEGER,DIMENSION(:),POINTER::p1
!!$  END TYPE iarray_pointer
!!$
!!$
!!$  INTEGER,INTENT(IN)::core
!!$ 
!!$  INTEGER::nats,nats2,nnods,ntot
!!$  INTEGER::aa,bb,cc,dd
!!$  INTEGER::ii,jj,kk
!!$
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::order_core,eff_order_core,vaux
!!$  INTEGER,DIMENSION(:,:),ALLOCATABLE::dim_tree_core
!!$  TYPE(iarray_pointer),DIMENSION(:,:),POINTER::tree_core
!!$  INTEGER::num_crits
!!$
!!$  IF (ALLOCATED(mss_ind_ats))   DEALLOCATE(mss_ind_ats)
!!$  IF (ALLOCATED(mss_ind_nods))  DEALLOCATE(mss_ind_nods)
!!$  IF (ALLOCATED(mss))           DEALLOCATE(mss)
!!$ 
!!$  !! Comienzo 1st shell
!!$
!!$  nats=num_ats_nod(core)
!!$  nats2=2*nats
!!$  nnods=num_Hbs_Bs_nod(core)
!!$
!!$  ALLOCATE(order_core(nats),eff_order_core(nats))
!!$
!!$  ii=in_at_nod(core)
!!$  eff_order_core=(/(jj,jj=1,nats)/)
!!$  order_core(:)=eff_order_core(:)+ii
!!$
!!$  ALLOCATE(tree_core(nats,2),dim_tree_core(nats,2))
!!$
!!$  DO ii=1,nats
!!$     jj=order_core(ii)
!!$     aa=num_Hbs_at(jj)
!!$     bb=num_Bs_at(jj)
!!$     cc=in_Hbs(jj)
!!$     dd=in_Bs(jj)
!!$     ALLOCATE(tree_core(ii,1)%p1(aa),tree_core(ii,2)%p1(bb))
!!$     dim_tree_core(ii,:)=(/aa,bb/)
!!$     tree_core(ii,1)%p1(:)=Hbs(cc+1:cc+aa)
!!$     tree_core(ii,2)%p1(:)= Bs(dd+1:dd+bb)
!!$  END DO
!!$
!!$  !! Quito Simetrias
!!$
!!$  DO ii=1,nats
!!$     DO jj=1,2
!!$        IF (dim_tree_core(ii,jj)>1) THEN
!!$           CALL build_order_bonded_1st(tree_core(ii,jj)%p1(:),dim_tree_core(ii,jj),num_crits)
!!$           IF (num_crits>0) THEN
!!$              DEALLOCATE(vect_aux)
!!$           END IF
!!$        END IF
!!$     END DO
!!$  END DO
!!$
!!$  !! Comienzo 2nd shell
!!$
!!$
!!$
!!$
!!$  ntot=1+nats+nats2+nnods
!!$  ALLOCATE(mss_ind_ats(ntot),mss_ind_nods(ntot),mss(ntot))
!!$  mss_ind_ats(:)=0
!!$  mss_ind_nods(:)=0
!!$  mss(:)=0
!!$
!!$  jj=1
!!$  mss_ind_ats(jj)  = nats
!!$  mss_ind_nods(jj) = nats
!!$  ii=jj+1
!!$  jj=jj+nats
!!$  mss_ind_ats(ii:jj)  = order_core(:)
!!$  mss_ind_nods(ii:jj) = at2nod(mss_ind_ats(ii:jj))
!!$  mss_ind_ats(ii:jj)  = trad2py_at(mss_ind_ats(ii:jj))
!!$  mss_ind_nods(ii:jj) = trad2py_nod(mss_ind_nods(ii:jj))
!!$  ii=jj+1
!!$  jj=jj+nats2
!!$  mss_ind_ats(ii:jj)  = (/(dim_tree_core(eff_order_core(kk),:),kk=1,nats)/)
!!$  mss_ind_nods(ii:jj) = mss_ind_ats(ii:jj)
!!$  ii=jj+1
!!$  jj=jj+nnods
!!$  mss_ind_ats(ii:jj) = (/((/tree_core(eff_order_core(kk),1)%p1(:),tree_core(eff_order_core(kk),2)%p1(:)/),kk=1,nats)/)
!!$  mss_ind_nods(ii:jj) = at2nod(mss_ind_ats(ii:jj))
!!$  mss_ind_ats(ii:jj)  = trad2py_at(mss_ind_ats(ii:jj))
!!$  mss_ind_nods(ii:jj) = trad2py_nod(mss_ind_nods(ii:jj))
!!$
!!$END SUBROUTINE build_shell2nd_old
!!$  


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



!!$SUBROUTINE order_ats_1st (core,nats,order)
!!$
!!$  IMPLICIT NONE
!!$
!!$  INTEGER,INTENT(IN)::core,nats
!!$  INTEGER,DIMENSION(nats),INTENT(OUT)::order
!!$
!!$  INTEGER::ii,jj,num_crits,dim_aux
!!$
!!$  jj=in_at_nod(core)
!!$  order(:)=(/(ii,ii=jj+1,jj+nats)/)
!!$
!!$  num_crits=ats_symm_num_crits(core)
!!$
!!$  IF (num_crits>0) THEN
!!$
!!$     dim_aux=ats_symm_length(core)
!!$     ii=ats_symm_in(core)
!!$     ALLOCATE(vect_aux(dim_aux))
!!$     vect_aux(:)=ats_symm(ii:(ii+dim_aux-1))
!!$     CALL SORTBYNUMHBS_ATS_1SH(core,nats,order,num_crits,dim_aux)
!!$     IF (num_crits>0) CALL SORTBYNUMBS_ATS_1SH(core,nats,order,num_crits,dim_aux)
!!$     IF (num_crits>0) DEALLOCATE(vect_aux)
!!$
!!$  END IF
!!$
!!$END SUBROUTINE order_ats_1st
!!$
!!$
!!$SUBROUTINE build_order_bonded_1st (order,dim_ord,num_crits)
!!$
!!$  INTEGER::dim_ord
!!$  INTEGER,DIMENSION(:),INTENT(INOUT)::order
!!$  INTEGER,INTENT(OUT)::num_crits
!!$
!!$  INTEGER::ii,dim_aux
!!$  INTEGER,DIMENSION(dim_ord)::lista_nodos
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::valores
!!$
!!$  INTEGER::dim_sop
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::sop
!!$  LOGICAL,DIMENSION(:),ALLOCATABLE::filt_sop
!!$  
!!$  lista_nodos(:)=at2nod(order(:))
!!$  
!!$  ALLOCATE(vect_aux(dim_ord+1),valores(dim_ord+1))
!!$  vect_aux(:)=(/dim_ord,(/(ii,ii=1,dim_ord)/)/)
!!$  valores(:)=(/0,(/lev_supsets(lista_nodos(:))/)/)
!!$
!!$  num_crits=1
!!$  dim_aux=dim_ord+1
!!$
!!$  CALL SORT_INT_1SH(dim_ord,dim_aux,order,valores,num_crits)
!!$
!!$  IF (superfiltro_sets.eqv..TRUE.) THEN
!!$     IF (num_crits>0) THEN
!!$        
!!$        !ALLOCATE(filt_sop(dim_aux))
!!$        !filt_sop
!!$
!!$     !IF (ANY(filtro_supsets(lista_nodos(:)),DIM=1).eqv..TRUE.) THEN
!!$     !IF (COUNT(filtro_supsets(lista_nodos(:)),DIM=1)>1) THEN
!!$        !print*,'si'
!!$        !dim_sop=dim_aux
!!$        !ALLOCATE(sop(dim_aux))
!!$        !sop(:)=vect_aux(:)
!!$        ! 
!!$        !print*, num_crits
!!$        !print*, vect_aux(:)
!!$
!!$     END IF
!!$  END IF
!!$
!!$  DEALLOCATE(valores)
!!$
!!$END SUBROUTINE build_order_bonded_1st
!!$
!!$
!!$
!!$SUBROUTINE SORTBYNUMHBS_ATS_1SH (core,nats,order,num_crits,dim_aux)
!!$
!!$  IMPLICIT NONE
!!$
!!$  INTEGER,INTENT(IN)::core,nats
!!$  INTEGER,INTENT(INOUT)::dim_aux
!!$  INTEGER,DIMENSION(nats),INTENT(INOUT)::order
!!$  INTEGER,INTENT(INOUT)::num_crits
!!$
!!$  INTEGER::ii,jj,gg,kk,ll,idim
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::valores
!!$
!!$  ALLOCATE(valores(dim_aux))
!!$
!!$  gg=0
!!$  DO ii=1,num_crits
!!$     gg=gg+1
!!$     idim=vect_aux(gg)
!!$     DO jj=1,idim
!!$        gg=gg+1
!!$        kk=vect_aux(gg)
!!$        ll=order(kk)
!!$        valores(gg)=num_hbs_at(ll)
!!$     END DO
!!$  END DO
!!$
!!$  CALL SORT_INT_1SH(nats,dim_aux,order,valores,num_crits)
!!$
!!$  DEALLOCATE(valores)
!!$
!!$END SUBROUTINE SORTBYNUMHBS_ATS_1SH
!!$
!!$SUBROUTINE SORTBYNUMBS_ATS_1SH (core,nats,order,num_crits,dim_aux)
!!$
!!$  IMPLICIT NONE
!!$
!!$  INTEGER,INTENT(IN)::core,nats
!!$  INTEGER,INTENT(INOUT)::dim_aux
!!$  INTEGER,DIMENSION(nats),INTENT(INOUT)::order
!!$  INTEGER,INTENT(INOUT)::num_crits
!!$
!!$  INTEGER::ii,jj,gg,kk,ll,idim
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::valores
!!$
!!$  ALLOCATE(valores(dim_aux))
!!$
!!$  gg=0
!!$  DO ii=1,num_crits
!!$     gg=gg+1
!!$     idim=vect_aux(gg)
!!$     DO jj=1,idim
!!$        gg=gg+1
!!$        kk=vect_aux(gg)
!!$        ll=order(kk)
!!$        valores(gg)=num_bs_at(ll)
!!$     END DO
!!$  END DO
!!$
!!$  CALL SORT_INT_1SH(nats,dim_aux,order,valores,num_crits)
!!$
!!$  DEALLOCATE(valores)
!!$
!!$END SUBROUTINE SORTBYNUMBS_ATS_1SH





END MODULE GLOB


