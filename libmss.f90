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
INTEGER,DIMENSION(:),ALLOCATABLE::trad2py_nod,trad2py_at


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
                                    
  IF (ALLOCATED(lev_ats))          DEALLOCATE(lev_ats)
  IF (ALLOCATED(lev_nods))         DEALLOCATE(lev_nods)
  IF (ALLOCATED(lev_sets))         DEALLOCATE(lev_sets)
  IF (ALLOCATED(lev_supsets))      DEALLOCATE(lev_supsets)
                                    
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
  ALLOCATE(lev_ats(num_ats),lev_nods(num_ats),lev_sets(num_ats),lev_supsets(num_ats))

  num_ats_nod(:)=0
  DO ii=1,num_ats
     jj=at2nod(ii)
     num_ats_nod(jj)=num_ats_nod(jj)+1
  END DO

  jj=0
  DO ii=1,num_nods
     in_at_nod(ii)=jj
     DO kk=1,num_ats_nod(ii)
        jj=jj+1
        lev_ats(jj)=kk
        lev_nods(jj)=1
        lev_sets(jj)=1
        lev_supsets(jj)=ii
     END DO
  END DO

  ! SYMM ATS

  DO ii=1,num_nods
     jj=xx_symm_ats_start(ii)
     DO kk=1,xx_symm_ats_crits(ii)
        hh=xx_symm_ats(jj)
        gg=xx_symm_ats(jj+1)
        DO ll=jj+1,jj+hh
           mm=in_at_nod(ii)+xx_symm_ats(ll)
           lev_ats(mm)=gg
        END DO
        jj=jj+hh+1
     END DO
  END DO

  ! SYMM NODS

  gg=1
  DO ii=1,xx_symm_nods_num
     hh=xx_symm_nods(gg)
     kk=xx_symm_nods(gg+1)
     kk=in_at_nod(kk)+1
     kk=lev_supsets(kk)
     DO jj=gg+1,gg+hh
        ll=xx_symm_nods(jj)
        DO mm=in_at_nod(ll)+1,in_at_nod(ll)+num_ats_nod(ll)
           lev_supsets(mm)=kk
        END DO
     END DO
     gg=gg+hh+1
  END DO
        
  ! SYMM SETS

  gg=1
  DO ii=1,xx_symm_sets_num
     jj=xx_symm_sets(gg)
     kk=xx_symm_sets(gg+1)
     ll=xx_symm_sets(gg+2)
     ll=in_at_nod(ll)+1
     ll=lev_supsets(ll)
     gg=gg+2
     DO mm=1,jj
        DO hh=1,kk
           nn=xx_symm_sets(gg)
           DO oo=in_at_nod(nn)+1,in_at_nod(nn)+num_ats_nod(nn)
              lev_supsets(oo)=ll
              lev_sets(oo)=mm
              lev_nods(oo)=hh
           END DO
           gg=gg+1
        END DO
     END DO
  END DO

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



!!!#################################
!!!####  BUILD FIRST SHELL
!!!#################################

SUBROUTINE build_shell1st (core)
 
  IMPLICIT NONE
 
  INTEGER,INTENT(IN)::core
 
  INTEGER::nats,nats2,nnods,ntot
  INTEGER::ff,gg,hh,ii,jj,kk,lhbs,lbs,mhbs,mbs
  INTEGER,DIMENSION(:),ALLOCATABLE::order_ats_1sh,num_ats_bonded_1sh,order_ats_bonded_1sh
  INTEGER,DIMENSION(:),ALLOCATABLE::order_nods_1sh,order_nods_bonded_1sh

  IF (ALLOCATED(mss_ind_ats)) DEALLOCATE(mss_ind_ats)
  IF (ALLOCATED(mss_ind_nods)) DEALLOCATE(mss_ind_nods)
  IF (ALLOCATED(mss))           DEALLOCATE(mss)
 
  nats=num_ats_nod(core)
  nats2=2*nats
  nnods=num_Hbs_Bs_nod(core)
 
  ALLOCATE(order_ats_1sh(nats),num_ats_bonded_1sh(nats2),order_ats_bonded_1sh(nnods))
 
  !CALL order_ats_1st(core,nn,order_ats_1sh,symm_ats_1sh)
  DO ii=1,nats
     order_ats_1sh(ii)=in_at_nod(core)+ii
  END DO
 
  gg=0
  ff=0
  DO ii=1,nats
     jj=order_ats_1sh(ii)
     lhbs=num_hbs_at(jj)
     lbs=num_bs_at(jj)
     mhbs=in_hbs(jj)
     mbs=in_bs(jj)
     gg=gg+1
     num_ats_bonded_1sh(gg)=lhbs
     gg=gg+1
     num_ats_bonded_1sh(gg)=lbs
     DO kk=1,lhbs
        hh=mhbs+kk
        ff=ff+1
        order_ats_bonded_1sh(ff)=hbs(hh)
     END DO
     DO kk=1,lbs
        hh=mbs+kk
        ff=ff+1
        order_ats_bonded_1sh(ff)=bs(hh)
     END DO
  END DO

 
  !Translate

  ALLOCATE(order_nods_1sh(nats),order_nods_bonded_1sh(nnods))

  order_nods_1sh(:)        =at2nod(order_ats_1sh(:))
  order_nods_bonded_1sh(:) =at2nod(order_ats_bonded_1sh(:))

  order_ats_1sh(:)          =trad2py_at(order_ats_1sh(:))
  order_ats_bonded_1sh(:)   =trad2py_at(order_ats_bonded_1sh(:))
  order_nods_1sh(:)         =trad2py_nod(order_nods_1sh(:))
  order_nods_bonded_1sh(:)  =trad2py_nod(order_nods_bonded_1sh(:))

  !Build

  ntot=1+nats+nats2+nnods
  ALLOCATE(mss_ind_ats(ntot),mss_ind_nods(ntot),mss(ntot))

  mss_ind_ats(1)  = nats
  mss_ind_nods(1) = nats
  ii=2
  jj=1+nats
  mss_ind_ats(ii:jj)  = order_ats_1sh(:)
  mss_ind_nods(ii:jj) = order_nods_1sh(:)
  ii=jj+1
  jj=jj+nats2
  mss_ind_ats(ii:jj)  = num_ats_bonded_1sh(:)
  mss_ind_nods(ii:jj) = num_ats_bonded_1sh(:)
  ii=jj+1
  jj=jj+nnods
  mss_ind_ats(ii:jj)  = order_ats_bonded_1sh(:)
  mss_ind_nods(ii:jj) = order_nods_bonded_1sh(:)

  DEALLOCATE(order_ats_1sh,num_ats_bonded_1sh,order_ats_bonded_1sh)
  DEALLOCATE(order_nods_1sh,order_nods_bonded_1sh)
 
END SUBROUTINE build_shell1st



END MODULE GLOB


