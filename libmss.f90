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





CONTAINS

SUBROUTINE load_topol(xx_at2nod,&
     xx_trad2py_at,xx_trad2py_nod,&
     xx_symm_ats_start,xx_symm_ats_crits,xx_symm_ats,&
     xx_num_ats,xx_num_nods,xx_symm_ats_dim)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::xx_num_ats,xx_num_nods,xx_symm_ats_dim
  INTEGER,DIMENSION(xx_num_ats),INTENT(IN)::xx_at2nod
  INTEGER,DIMENSION(xx_num_ats),INTENT(IN)::xx_trad2py_at
  INTEGER,DIMENSION(xx_num_nods),INTENT(IN)::xx_trad2py_nod
  INTEGER,DIMENSION(xx_num_nods),INTENT(IN)::xx_symm_ats_start
  INTEGER,DIMENSION(xx_num_nods),INTENT(IN)::xx_symm_ats_crits
  INTEGER,DIMENSION(xx_symm_ats_dim),INTENT(IN)::xx_symm_ats

  INTEGER::gg,hh,ii,jj,kk,ll,mm

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


END SUBROUTINE load_net


END MODULE GLOB


