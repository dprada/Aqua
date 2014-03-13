!!!###################################################################################################
!!!####  VARIABLES
!!!###################################################################################################

MODULE GLOB

  IMPLICIT NONE

TYPE p_symm
   INTEGER::num_crits,length
   INTEGER,DIMENSION(:),ALLOCATABLE::crit
END TYPE p_symm

TYPE p_bonded
   INTEGER,DIMENSION(:),ALLOCATABLE::ats,nods,aux2sh
   INTEGER,DIMENSION(:),ALLOCATABLE::order,wsymm
   INTEGER,DIMENSION(:),ALLOCATABLE::lev_supsets
   INTEGER::num

   LOGICAL::unsolved_loops
   INTEGER,DIMENSION(:),ALLOCATABLE::loops_in,loops_same,loops_ant,loops_post
   LOGICAL::with_loops_in,with_loops_same,with_loops_ant,with_loops_post

   INTEGER::num_permutations
   INTEGER,DIMENSION(:,:),ALLOCATABLE::permutations

   TYPE(p_symm)::symm
END TYPE p_bonded

TYPE p_at
   TYPE(p_bonded),DIMENSION(2)::bonded
   LOGICAL::empty
   INTEGER::ind
END TYPE p_at

TYPE p_shell

   TYPE(p_at),DIMENSION(:),ALLOCATABLE::ats
   INTEGER,DIMENSION(:),ALLOCATABLE::order
   INTEGER::nats,nats2,nnods,ntot
   INTEGER::ind
   INTEGER,DIMENSION(:),ALLOCATABLE::wsymm
   TYPE(p_symm)::symm
   INTEGER::lev_supsets

   INTEGER::num_permutations
   INTEGER,DIMENSION(:,:),ALLOCATABLE::permutations

   LOGICAL::unsolved_loops
   INTEGER::nnods_norepe,nnods_repe
   INTEGER,DIMENSION(:),ALLOCATABLE::nods_norepe,nods_repe,nods_cant_repe

END TYPE p_shell

TYPE p_trans_set
   INTEGER::tope,num_nods
   INTEGER::count_sets
END type p_trans_set

TYPE p_translation
   TYPE(p_trans_set),DIMENSION(:),ALLOCATABLE::supset
   INTEGER::num_supsets
   LOGICAL,DIMENSION(:),ALLOCATABLE::total_filtro
   INTEGER,DIMENSION(:),ALLOCATABLE::dict
END type p_translation

!f2py   intent(hide)::list_shells
TYPE(p_shell),DIMENSION(:),ALLOCATABLE,TARGET::list_shells

!! TOPOLOGY

!f2py   intent(hide)::num_ats,num_nods,num_supsets
!f2py   intent(hide)::at2nod,lev_supsets
!f2py   intent(hide)::trad2py_nod,trad2py_at
INTEGER::num_ats,num_nods,num_supsets
INTEGER,DIMENSION(:),ALLOCATABLE::at2nod
INTEGER,DIMENSION(:),ALLOCATABLE::lev_supsets
INTEGER,DIMENSION(:),ALLOCATABLE::trad2py_nod,trad2py_at



!!! AUXILIAR

!f2py   intent(hide)::symm,translator
TYPE(p_symm)::symm
TYPE(p_translation),TARGET::translator

!f2py intent(hide)::pff2,precalc2,precalc2_num
!f2py intent(hide)::pff3,precalc3,precalc3_num
!f2py intent(hide)::pff4,precalc4,precalc4_num
INTEGER,PARAMETER::precalc2_num=2
INTEGER,PARAMETER::precalc3_num=6
INTEGER,PARAMETER::precalc4_num=24
INTEGER,DIMENSION(2*2)::&
     pff2=(/1,2, 2,1/)
INTEGER,DIMENSION(3*6)::&
     pff3=(/1,2,3, 1,3,2, 2,1,3,&
            2,3,1, 3,1,2, 3,2,1/)
INTEGER,DIMENSION(4*24)::&
     pff4=(/1,2,3,4, 1,2,4,3, 1,3,2,4, 1,3,4,2, 1,4,2,3, 1,4,3,2,&
            2,1,3,4, 2,1,4,3, 2,3,1,4, 2,3,4,1, 2,4,1,3, 2,4,3,1,&
            3,1,2,4, 3,1,4,2, 3,2,1,4, 3,2,4,1, 3,4,1,2, 3,4,2,1,&
            4,1,2,3, 4,1,3,2, 4,2,1,3, 4,2,3,1, 4,3,1,2, 4,3,2,1/)
INTEGER,DIMENSION(2,2)::precalc2
INTEGER,DIMENSION(3,6)::precalc3
INTEGER,DIMENSION(4,24)::precalc4
equivalence(pff2,precalc2)
equivalence(pff3,precalc3)
equivalence(pff4,precalc4)

!!! OUTPUT

INTEGER,DIMENSION(:),ALLOCATABLE::mss_ind_ats,mss_ind_nods,mss,mss_symm
!f2py intent(hide)::mss_template
LOGICAL,DIMENSION(:),ALLOCATABLE::mss_template

CONTAINS

!!!###################################################################################################
!!!####  AUXILIAR TOPOL AND NET SUBROUTINES
!!!###################################################################################################

SUBROUTINE SHELL_DOWN(shell)

  IMPLICIT NONE

  TYPE(p_shell),INTENT(INOUT)::shell

  shell%nnods=0
  shell%ntot=0
  shell%unsolved_loops=.FALSE.
  IF (ALLOCATED(shell%nods_norepe)) THEN
     DEALLOCATE(shell%nods_norepe)
     DEALLOCATE(shell%nods_repe)
     DEALLOCATE(shell%nods_cant_repe)
  END IF


END SUBROUTINE SHELL_DOWN

SUBROUTINE AT_DOWN(at)
 
  IMPLICIT NONE
 
  TYPE(p_at),TARGET,INTENT(INOUT)::at
  TYPE(p_bonded),POINTER::bonded
  INTEGER::ii

  DO ii=1,2
     bonded=>at%bonded(ii)
     IF (ALLOCATED(bonded%ats)) THEN
        DEALLOCATE(bonded%ats)
        DEALLOCATE(bonded%nods)
        DEALLOCATE(bonded%order)
        DEALLOCATE(bonded%lev_supsets)
        DEALLOCATE(bonded%loops_in)
        DEALLOCATE(bonded%loops_same)
        DEALLOCATE(bonded%loops_post)
        DEALLOCATE(bonded%loops_ant)
        DEALLOCATE(bonded%aux2sh)
        DEALLOCATE(bonded%wsymm)
        DEALLOCATE(bonded%symm%crit)
     END IF
  END DO

  NULLIFY(bonded)

END SUBROUTINE AT_DOWN

SUBROUTINE AT_UP(at)
 
  IMPLICIT NONE
 
  TYPE(p_at),TARGET,INTENT(INOUT)::at
  TYPE(p_bonded),POINTER::bonded
  INTEGER::ii,kk,num

  at%empty=.TRUE.

  DO ii=1,2
     bonded=>at%bonded(ii)
     num=bonded%num
     ALLOCATE(bonded%ats(num))
     ALLOCATE(bonded%nods(num))
     ALLOCATE(bonded%order(num))
     ALLOCATE(bonded%lev_supsets(num))
     ALLOCATE(bonded%loops_in(num))
     ALLOCATE(bonded%loops_same(num))
     ALLOCATE(bonded%loops_post(num))
     ALLOCATE(bonded%loops_ant(num))
     ALLOCATE(bonded%aux2sh(num)) 
     ALLOCATE(bonded%wsymm(num)) 
     ALLOCATE(bonded%symm%crit(num+1)) !!!
     bonded%order                =(/(kk,kk=1,num)/)
     bonded%loops_in(:)          = 0
     bonded%loops_same(:)        = 0
     bonded%loops_post(:)        = 0
     bonded%loops_ant(:)         = 0
     bonded%aux2sh(:)            = 0
     bonded%wsymm(:)             = 0
     bonded%symm%length          = num+1
     IF (num>1) THEN
        bonded%symm%num_crits    = 1
        bonded%symm%crit(:)      = (/num,bonded%order(:)/)
     ELSE
        bonded%symm%num_crits    = 0
        bonded%symm%crit(:)      = 0
     END IF
     bonded%unsolved_loops            =.FALSE.
     bonded%with_loops_in             =.FALSE.
     bonded%with_loops_same           =.FALSE.
     bonded%with_loops_post           =.FALSE.
     bonded%with_loops_ant            =.FALSE.
  END DO

  NULLIFY(bonded)

END SUBROUTINE AT_UP

SUBROUTINE DOWNLOAD_SYMM(loc_symm)
 
  IMPLICIT NONE
 
  TYPE(p_symm),INTENT(INOUT)::loc_symm
 
  loc_symm%num_crits=symm%num_crits
  loc_symm%length=symm%length
 
  DEALLOCATE(loc_symm%crit)
  IF (symm%num_crits>0) THEN
     ALLOCATE(loc_symm%crit(loc_symm%length))
     loc_symm%crit(:)=symm%crit(:)
  ELSE
     ALLOCATE(loc_symm%crit(0))
  END IF
 
  DEALLOCATE(symm%crit)
 
END SUBROUTINE DOWNLOAD_SYMM

SUBROUTINE UPLOAD_SYMM(loc_symm)
 
  IMPLICIT NONE
 
  TYPE(p_symm),INTENT(IN)::loc_symm
 
  symm%num_crits=loc_symm%num_crits
  symm%length=loc_symm%length
 
  ALLOCATE(symm%crit(symm%length))
  symm%crit(:)=loc_symm%crit(:)
 
END SUBROUTINE UPLOAD_SYMM



SUBROUTINE RESET_TRANSLATOR()
 
  IMPLICIT NONE
 
  INTEGER::ii
 
  DO ii=1,translator%num_supsets
     translator%supset(ii)%count_sets=0
  END DO
 
  translator%total_filtro=.FALSE.

END SUBROUTINE RESET_TRANSLATOR


SUBROUTINE BUILD_LOOPS_IN(shell)

  IMPLICIT NONE

  TYPE(p_shell),TARGET,INTENT(INOUT)::shell

  TYPE(p_at),POINTER::at
  TYPE(p_bonded),POINTER::bonded

  INTEGER::ii,jj,kk,ll,mm
  INTEGER,DIMENSION(:),ALLOCATABLE::carro

  ALLOCATE(carro(shell%nnods))
  carro(:)=(/((/shell%ats(kk)%bonded(1)%nods(:),shell%ats(kk)%bonded(2)%nods(:)/),kk=1,shell%nats)/)

  CALL BUILDING_REPE_LISTS (carro,shell%nnods,shell%nods_norepe,shell%nnods_norepe,shell%nods_repe,&
       shell%nods_cant_repe,shell%nnods_repe)

  IF (shell%nnods_repe>0) THEN
     DO mm=1,shell%nats
        at=>shell%ats(mm)
        DO ll=1,2
           IF (at%bonded(ll)%num>0) THEN
              bonded=>at%bonded(ll)
              DO ii=1,shell%nnods_repe
                 jj=shell%nods_repe(ii)
                 DO kk=1,bonded%num
                    IF (jj==bonded%nods(kk)) THEN
                       bonded%loops_in(kk)=shell%nods_cant_repe(ii)
                       bonded%with_loops_in=.TRUE.
                    END IF
                 END DO
              END DO
           END IF
        END DO
     END DO
  END IF

  DEALLOCATE(carro)
  NULLIFY(at,bonded)

END SUBROUTINE BUILD_LOOPS_IN


!!!###################################################################################################
!!!####  LOADING TOPOL AND NET
!!!###################################################################################################


SUBROUTINE load_topol(xx_at2nod,&
     xx_trad2py_at,xx_trad2py_nod,&
     xx_symm_ats_start,xx_symm_ats_crits,xx_symm_ats,&
     xx_symm_nods,xx_symm_nods_num,&
     xx_num_ats,xx_num_nods,&
     xx_symm_ats_dim,xx_symm_nods_dim)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::xx_num_ats,xx_num_nods,xx_symm_nods_num
  INTEGER,INTENT(IN)::xx_symm_ats_dim,xx_symm_nods_dim
  INTEGER,DIMENSION(xx_num_ats),INTENT(IN)::xx_at2nod
  INTEGER,DIMENSION(xx_num_ats),INTENT(IN)::xx_trad2py_at
  INTEGER,DIMENSION(xx_num_nods),INTENT(IN)::xx_trad2py_nod
  INTEGER,DIMENSION(xx_num_nods),INTENT(IN)::xx_symm_ats_start
  INTEGER,DIMENSION(xx_num_nods),INTENT(IN)::xx_symm_ats_crits
  INTEGER,DIMENSION(xx_symm_ats_dim),INTENT(IN)::xx_symm_ats
  INTEGER,DIMENSION(xx_symm_nods_dim),INTENT(IN)::xx_symm_nods

  INTEGER::gg,hh,ii,jj,kk,ll,mm
  INTEGER,DIMENSION(:),ALLOCATABLE::num_ats_nod
  TYPE(p_shell),POINTER::shell


  num_ats     = xx_num_ats
  num_nods    = xx_num_nods


  IF (ALLOCATED(at2nod))           DEALLOCATE(at2nod)
  IF (ALLOCATED(trad2py_at))       DEALLOCATE(trad2py_at)
  IF (ALLOCATED(trad2py_nod))      DEALLOCATE(trad2py_nod)
  IF (ALLOCATED(list_shells))      DEALLOCATE(list_shells)
  IF (ALLOCATED(lev_supsets))      DEALLOCATE(lev_supsets)

  ALLOCATE(at2nod(num_ats),trad2py_at(num_ats),trad2py_nod(num_nods))
  ALLOCATE(lev_supsets(num_nods))
  ALLOCATE(list_shells(num_nods))

  ALLOCATE(num_ats_nod(num_nods))


  at2nod(:)      = xx_at2nod(:)
  trad2py_at(:)  = xx_trad2py_at(:)
  trad2py_nod(:) = xx_trad2py_nod(:)

  num_ats_nod(:)=0

  DO ii=1,num_ats
     jj=at2nod(ii)
     num_ats_nod(jj)=num_ats_nod(jj)+1
  END DO

  DO ii=1,num_nods
     lev_supsets(ii)=ii
  END DO


  gg=0
  DO ii=1,num_nods
     kk=num_ats_nod(ii)
     shell=>list_shells(ii)
     shell%ind=ii
     shell%nats=kk
     shell%nats2=kk*2
     shell%nnods=0
     shell%ntot=0
     ALLOCATE(shell%order(kk))
     ALLOCATE(shell%ats(kk))
     ALLOCATE(shell%wsymm(kk))
     shell%order(:)=(/(jj,jj=1,kk)/)
     shell%wsymm(:)=0
     DO jj=1,kk
        gg=gg+1
        shell%order(jj)=jj
        shell%ats(jj)%ind=gg
     END DO
     jj=xx_symm_ats_crits(ii)
     shell%symm%num_crits=jj
     kk=xx_symm_ats_start(ii)-1
     mm=0
     DO ll=1,jj
        mm=mm+1
        mm=mm+xx_symm_ats(kk+mm)
     END DO
     shell%symm%length=mm
     ALLOCATE(shell%symm%crit(mm))
     shell%symm%crit(:)=xx_symm_ats((kk+1):(kk+mm))
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
        
  ! REAJUSTO LEVELS en shells

  DO ii=1,num_nods
     list_shells(ii)%lev_supsets = lev_supsets(ii)
  END DO

  num_supsets=1
  DO ii=1,num_nods-1
     IF (lev_supsets(ii)/=(lev_supsets(ii+1))) THEN
        lev_supsets(ii)=num_supsets
        num_supsets=num_supsets+1
     ELSE
        lev_supsets(ii)=num_supsets
     END IF
  END DO
  lev_supsets(num_nods)=num_supsets

  ! TRANSLATOR

  ALLOCATE(translator%supset(num_supsets))
  ALLOCATE(translator%total_filtro(num_nods))
  ALLOCATE(translator%dict(num_nods))
  translator%num_supsets=num_supsets
  DO ii=1,num_supsets
     translator%supset(ii)%num_nods=0
  END DO

  DO ii=1,num_nods
     jj=lev_supsets(ii)
     translator%supset(jj)%num_nods=translator%supset(jj)%num_nods+1
  END DO

  jj=0
  DO ii=1,num_supsets
     translator%supset(ii)%tope=jj
     jj=jj+translator%supset(ii)%num_nods
  END DO

  NULLIFY(shell)
  DEALLOCATE(num_ats_nod)

END SUBROUTINE LOAD_TOPOL


SUBROUTINE LOAD_NET(xx_hbs,xx_bs,xx_num_Hbs_at,xx_num_Bs_at,&
     xx_Total_num_hbs,xx_Total_num_bs,xx_num_ats)
 
  IMPLICIT NONE
 
  INTEGER,INTENT(IN)::xx_num_ats,xx_Total_num_hbs,xx_Total_num_bs
  INTEGER,DIMENSION(xx_num_ats),INTENT(IN)::xx_num_hbs_at,xx_num_bs_at
  INTEGER,DIMENSION(xx_Total_num_hbs),INTENT(IN)::xx_hbs
  INTEGER,DIMENSION(xx_Total_num_bs),INTENT(IN)::xx_bs
 
  INTEGER::ii,jj,kk,ll,gghb,ggb,ggat,numhbs,numbs
  INTEGER,DIMENSION(:),ALLOCATABLE::vect_aux
  TYPE(p_shell),POINTER::shell
  TYPE(p_bonded),POINTER::bonded
  TYPE(p_at),POINTER::at
 
  gghb=0
  ggb=0
  DO ii=1,num_nods
     shell=>list_shells(ii)
     CALL SHELL_DOWN (shell)
     DO jj=1,shell%nats
        at=>shell%ats(jj)
        CALL AT_DOWN(at)
        ggat=at%ind
        numhbs=xx_num_hbs_at(ggat)
        numbs=xx_num_bs_at(ggat)
        at%bonded(1)%num=numhbs
        at%bonded(2)%num=numbs
        CALL AT_UP(at)
        !HBS
        kk=gghb+1
        gghb=gghb+numhbs
        at%bonded(1)%ats(:)=xx_hbs(kk:gghb)
        !BS
        kk=ggb+1
        ggb=ggb+numbs
        at%bonded(2)%ats(:)=xx_bs(kk:ggb)
        !! RECOLOCO
        IF ((numhbs+numbs)>0) at%empty=.FALSE.
        DO ll=1,2
           bonded=>at%bonded(ll)
           shell%nnods=shell%nnods+bonded%num
           ALLOCATE(vect_aux(bonded%num))
           vect_aux(:)                 = at2nod(bonded%ats(:))
           bonded%nods(:)              = vect_aux(:)
           bonded%lev_supsets(:)       = lev_supsets(vect_aux(:))
           DEALLOCATE(vect_aux)
        END DO
     END DO
     shell%ntot=1+shell%nats+shell%nats2+shell%nnods
     CALL BUILD_LOOPS_IN(shell)
  END DO
 
  NULLIFY(shell,bonded,at)
  
END SUBROUTINE LOAD_NET
 

!!!###################################################################################################
!!!####  BUILDING SHELLS
!!!###################################################################################################

SUBROUTINE BUILD_SHELL1ST (core)
 
  IMPLICIT NONE
 
  INTEGER,INTENT(IN)::core
  TYPE(p_shell),TARGET::shell
  TYPE(p_at),POINTER::at
  TYPE(p_bonded),POINTER::bonded
 
  INTEGER,DIMENSION(:),ALLOCATABLE::privilegios
  INTEGER::dim_privil,ii,ll
 
 
  shell=list_shells(core)
   
  !! order bonds
   
  dim_privil=1
  ALLOCATE(privilegios(1))
  privilegios(1)=core
 
  DO ii=1,shell%nats
     at=>shell%ats(ii)
     DO ll=1,2
        IF (at%bonded(ll)%symm%num_crits>0) THEN
           bonded=>at%bonded(ll)
           CALL upload_symm(bonded%symm)
           CALL order_bonded(dim_privil,privilegios,bonded)
           CALL download_symm(bonded%symm)
        END IF
     END DO
  END DO
 
  !! order ats

  IF (shell%symm%num_crits>0) THEN
     CALL upload_symm(shell%symm)
     CALL remove_triv_symm_ats(shell)
     IF (symm%num_crits>0) THEN
        CALL order_ats(dim_privil,privilegios,shell)
     END IF
     CALL download_symm(shell%symm)
  END IF
  
  DEALLOCATE(privilegios)

  CALL build_mss_ats(shell)
  CALL build_mss_nods()
  CALL build_mss_symm(shell)
  CALL build_mss()
  CALL build_trans_mss_ats()
  CALL build_trans_mss_nods()

   
  NULLIFY(at,bonded)
 
END SUBROUTINE BUILD_SHELL1ST
!############################

SUBROUTINE BUILD_SHELL2ND (core)
 
  IMPLICIT NONE
 
  INTEGER,INTENT(IN)::core
 
  TYPE(p_shell),TARGET::shell1st
  TYPE(p_shell),DIMENSION(:),ALLOCATABLE,TARGET::shell2nd
 
  TYPE(p_bonded),POINTER::bonded
  TYPE(p_at),POINTER::at
  TYPE(p_shell),POINTER::shell
 
  INTEGER::nnods,totntot,core2
  INTEGER::ntot1sh,nnods1sh,ntot2sh,nnods2sh
  INTEGER::ii,jj,kk,gg,iii,ggg,ll
  INTEGER,DIMENSION(:),ALLOCATABLE::privilegios
  INTEGER::dim_privil,dim_repe,dim_carro1,dim_carro2
  INTEGER,DIMENSION(:),ALLOCATABLE::aux_ind_ats,aux_ind_nods,aux_ord,aux_symm
  INTEGER,DIMENSION(:),ALLOCATABLE::carro,carro_aux,cant_repe1,cant_repe2
 
 
  !! Building the structure
 
  shell1st=list_shells(core)
 
  nnods=shell1st%nnods
  totntot=shell1st%ntot

  ALLOCATE(shell2nd(shell1st%nnods))
  gg=0
  DO ii=1,shell1st%nats
     DO ll=1,2
        bonded=>shell1st%ats(ii)%bonded(ll)
        DO jj=1,bonded%num
           gg=gg+1
           bonded%aux2sh(jj)=gg
           kk=bonded%nods(jj)
           shell2nd(gg)=list_shells(kk)
           totntot=totntot+shell2nd(gg)%ntot
        END DO
     END DO
  END DO
 
  !! Building loops_same 2nd

  gg=0
  DO iii=1,nnods
     gg=gg+shell2nd(iii)%nnods_norepe
  END DO
  ALLOCATE(carro(gg),carro_aux(gg),cant_repe2(gg))
  ggg=0
  DO iii=1,nnods
     gg=ggg+1
     ggg=ggg+shell2nd(iii)%nnods_norepe
     carro(gg:ggg)=shell2nd(iii)%nods_norepe(:)
     carro_aux(gg:ggg)=iii
  END DO
  dim_carro2=ggg

  CALL BUILDING_REPE_LIST2 (core,carro,dim_carro2,cant_repe2)

  DO kk=1,dim_carro2
     IF (cant_repe2(kk)>1) THEN
        gg=carro(kk)
        shell=>shell2nd(carro_aux(kk))
        DO ii=1,shell%nats
           at=>shell%ats(ii)
           DO ll=1,2
              bonded=>at%bonded(ll)
              DO jj=1,bonded%num
                 IF (bonded%nods(jj)==gg) THEN
                    bonded%loops_same(jj)=cant_repe2(kk)
                    bonded%with_loops_same=.TRUE.
                 END IF
              END DO
           END DO
        END DO
     END IF
  END DO

  !! Building loops_post 1st and loops_ant 2nd

  dim_carro1=shell1st%nnods_norepe
  ALLOCATE(cant_repe1(dim_carro1))
  CALL BUILDING_REPE_LIST3 (core,shell1st%nods_norepe,carro,dim_carro1,dim_carro2,cant_repe1,cant_repe2)

  DO kk=1,dim_carro2
     IF (cant_repe2(kk)>0) THEN
        gg=carro(kk)
        shell=>shell2nd(carro_aux(kk))
        DO ii=1,shell%nats
           at=>shell%ats(ii)
           DO ll=1,2
              bonded=>at%bonded(ll)
              DO jj=1,bonded%num
                 IF (bonded%nods(jj)==gg) THEN
                    bonded%loops_ant(jj)=cant_repe2(kk)
                    bonded%with_loops_ant=.TRUE.
                 END IF
              END DO
           END DO
        END DO
     END IF
  END DO

  DO kk=1,dim_carro1
     IF (cant_repe1(kk)>0) THEN
        gg=shell1st%nods_norepe(kk)
        DO ii=1,shell1st%nats
           at=>shell1st%ats(ii)
           DO ll=1,2
              bonded=>at%bonded(ll)
              DO jj=1,bonded%num
                 IF (bonded%nods(jj)==gg) THEN
                    bonded%loops_post(jj)=cant_repe1(kk)
                    bonded%with_loops_post=.TRUE.
                 END IF
              END DO
           END DO
        END DO
     END IF
  END DO

  DEALLOCATE(carro,cant_repe1,cant_repe2)
 
!  !! removing symm in 2ndsh
! 
!  dim_privil=2
!  ALLOCATE(privilegios(2))
!  privilegios(2)=core
! 
!  DO iii=1,nnods
! 
!     shell=>shell2nd(iii)
!     core2=shell%ind
! 
!     !list of privilegios
!     privilegios(1)=core2 
! 
!     ! order bonds 2sh
!     DO ii=1,shell%nats
!        at=>shell%ats(ii)
!        DO ll=1,2
!           IF (at%bonded(ll)%symm%num_crits>0) THEN
!              bonded=>at%bonded(ll)
!              CALL upload_symm(bonded%symm)
!              CALL order_bonded(dim_privil,privilegios,bonded)
!              IF (symm%num_crits>0) THEN
!                 CALL order_bonded_w_23order(bonded)
!              END IF
!              IF (symm%num_crits>0) THEN 
!                 !! quito falsas simetrias porque me quede sin criterios para los nodos que aparecen una sola vez
!                 CALL quito_falsos_bonded_symm(bonded)
!              END IF
!              CALL download_symm(bonded%symm)
!           END IF
!        END DO
!     END DO
! 
!     ! order ats 2sh
!     IF (shell%symm%num_crits>0) THEN
!        CALL upload_symm(shell%symm)
!        CALL remove_triv_symm_ats(shell) ! el que no tiene nada no tiene ni privilegios
!        IF (symm%num_crits>0) THEN
!           CALL order_ats(dim_privil,privilegios,shell)
!        END IF
!        IF (symm%num_crits>0) THEN
!           CALL order_ats_w_23order(shell)
!        END IF
!        IF (symm%num_crits>0) THEN
!           CALL quito_falsos_ats_symm(shell)
!        END IF
!        CALL download_symm(shell%symm)
!     END IF
! 
!  END DO
! 
!  !! removing symm in 1stsh
! 
!  DEALLOCATE(privilegios)
!  dim_privil=1
!  ALLOCATE(privilegios(1))
!  privilegios(1)=core
! 
!  ! order bonds 1sh
!  DO ii=1,shell1st%nats
!     at=>shell1st%ats(ii)
!     DO ll=1,2
!        IF (at%bonded(ll)%symm%num_crits>0) THEN
!           bonded=>at%bonded(ll)
!           CALL upload_symm(bonded%symm)
!           CALL order_bonded(dim_privil,privilegios,bonded)
!           IF (symm%num_crits>0) THEN
!              CALL order_bonded_w_2order(bonded)
!           END IF
!           IF (symm%num_crits>0) THEN
!              CALL order_bonded_w_next_shells(bonded,shell2nd,nnods)
!           END IF
!           IF (symm%num_crits>0) THEN
!              CALL quito_falsos_bonded_symm2(bonded,shell2nd,nnods)
!           END IF
!           CALL download_symm(bonded%symm)
!        END IF
!     END DO
!  END DO
!  
!  ! order ats 1sh
!  IF (shell1st%symm%num_crits>0) THEN
!     CALL upload_symm(shell1st%symm)
!     CALL remove_triv_symm_ats(shell1st) ! el que no tiene nada no tiene ni privilegios
!     IF (symm%num_crits>0) THEN 
!        CALL order_ats(dim_privil,privilegios,shell1st)
!        IF (symm%num_crits>0) THEN 
!           CALL order_ats_w_2order(shell1st)
!        END IF
!        IF (symm%num_crits>0) THEN 
!           CALL order_ats_w_next_shells(shell1st,shell2nd,nnods)
!        END IF
!        IF (symm%num_crits>0) THEN
!           CALL quito_falsos_ats_symm2(shell1st,shell2nd,nnods)
!        END IF
!     END IF
!     CALL download_symm(shell1st%symm)
!  END IF
! 
!  !! Building the msss
!  ! Poner aqui un if
!  CALL solving_last_symmetries_permuting (shell1st,shell2nd,nnods,totntot)
! 
! 
!  IF (ALLOCATED(mss_ind_ats))   DEALLOCATE(mss_ind_ats)
!  IF (ALLOCATED(mss_ind_nods))  DEALLOCATE(mss_ind_nods)
!  IF (ALLOCATED(mss_symm))      DEALLOCATE(mss_symm)
! 
!  ALLOCATE(mss_ind_ats(totntot),mss_ind_nods(totntot),mss_symm(totntot))
! 
!  gg=0
!  ntot1sh=shell1st%ntot
!  nnods1sh=shell1st%nnods
!  ALLOCATE(aux_ind_ats(ntot1sh),aux_ind_nods(ntot1sh),aux_symm(ntot1sh),aux_ord(nnods1sh))
!  CALL build_mss_wout_word(shell1st,aux_ind_ats,aux_ind_nods,aux_symm,aux_ord)
!  ggg=gg+1
!  gg=gg+ntot1sh
!  mss_ind_ats(ggg:gg)=aux_ind_ats(:)
!  mss_ind_nods(ggg:gg)=aux_ind_nods(:)
!  mss_symm(ggg:gg)=aux_symm(:)
!  DEALLOCATE(aux_ind_ats,aux_ind_nods,aux_symm)
! 
!  DO ii=1,nnods1sh
!     jj=aux_ord(ii)
!     shell=>shell2nd(jj)
!     ntot2sh =shell%ntot
!     nnods2sh=shell%nnods
!     ALLOCATE(aux_ind_ats(ntot2sh),aux_ind_nods(ntot2sh),aux_symm(ntot2sh))
!     CALL build_mss_wout(shell,aux_ind_ats,aux_ind_nods,aux_symm)
!     ggg=gg+1
!     gg=gg+ntot2sh
!     mss_ind_ats(ggg:gg)=aux_ind_ats(:)
!     mss_ind_nods(ggg:gg)=aux_ind_nods(:)
!     mss_symm(ggg:gg)=aux_symm(:)
!     DEALLOCATE(aux_ind_ats,aux_ind_nods,aux_symm)
!  END DO
! 
!  CALL remato_mss2sh()
! 
!  DEALLOCATE(aux_ord,shell2nd,privilegios)
!  NULLIFY(at,shell,bonded)
 
END SUBROUTINE BUILD_SHELL2ND



!!!###################################################################################################
!!!####  ORDER BONDED FUNCTIONS
!!!###################################################################################################

SUBROUTINE ORDER_BONDED(dim_privil,privilegios,bonded)
  
  IMPLICIT NONE
 
  INTEGER,INTENT(IN)::dim_privil
  INTEGER,DIMENSION(dim_privil),INTENT(IN)::privilegios
  TYPE(p_bonded),INTENT(INOUT)::bonded
  
  INTEGER::ii,gg,hh,priv,numb
  INTEGER,DIMENSION(:,:),ALLOCATABLE::valores
 
  numb=bonded%num

  gg=dim_privil+1
  ALLOCATE(valores(numb,gg))
  valores(:,:)=0
   
  DO ii=1,dim_privil
     priv=privilegios(ii)
     DO hh=1,numb
        IF (priv==bonded%nods(hh)) THEN
           valores(hh,ii)=1
        END IF
     END DO
  END DO
  valores(:,gg)=bonded%lev_supsets(:)

  CALL SORT_INT_MATRIX (numb,gg,bonded%order,valores)
  DEALLOCATE(valores)

  IF (symm%num_crits>0) THEN

     gg=1
     ALLOCATE(valores(numb,gg))
     valores(:,1)=bonded%loops_in(:)
     CALL SORT_INT_MATRIX (numb,gg,bonded%order,valores)
     DEALLOCATE(valores)
 
  END IF
 
END SUBROUTINE ORDER_BONDED

!!!###################################################################################################
!!!####  ORDER ATS FUNCTIONS
!!!###################################################################################################

SUBROUTINE REMOVE_TRIV_SYMM_ATS(shell)
 
  TYPE(p_shell),TARGET,INTENT(INOUT)::shell
 
  INTEGER::numats
  INTEGER::ii,gg,ll,idim,mm,nn
  LOGICAL::interruptor
  INTEGER::new_length,new_num_crits
  INTEGER,DIMENSION(:),ALLOCATABLE::new_crit
  INTEGER,DIMENSION(:,:),ALLOCATABLE::valores
 
  numats=shell%nats
 
  ALLOCATE(valores(numats,1))
 
  DO ii=1,numats
     IF (shell%ats(ii)%empty.eqv..FALSE.) THEN
        valores(ii,1)=1
     ELSE
        valores(ii,1)=0
     END IF
  END DO
 
  CALL SORT_INT_MATRIX (numats,1,shell%order,valores)
 
  IF (symm%num_crits>0) THEN
 
     ALLOCATE(new_crit(symm%length))
     new_num_crits=0
     new_length=0
     mm=0
 
     interruptor=.FALSE.
     gg=0
     DO ii=1,symm%num_crits
        gg=gg+1
        idim=symm%crit(gg)
        ll=shell%order(symm%crit(gg+1))
        IF (valores(ll,1)==0) THEN
           interruptor=.true.
           gg=gg+idim
        ELSE
           new_num_crits=new_num_crits+1
           new_length=new_length+1+idim
           mm=mm+1
           new_crit(mm)=idim
           DO nn=1,idim
              gg=gg+1
              mm=mm+1
              new_crit(mm)=symm%crit(gg)
           END DO
        END IF
     END DO
 
     IF (interruptor.eqv..TRUE.) THEN
        symm%num_crits=new_num_crits
        symm%length=new_length
        DEALLOCATE(symm%crit)
        ALLOCATE(symm%crit(new_length))
        symm%crit(:)=new_crit(1:new_length)
     END IF
 
     DEALLOCATE(new_crit)

  END IF
 
  DEALLOCATE(valores)
 
END SUBROUTINE  REMOVE_TRIV_SYMM_ATS
!############################


SUBROUTINE ORDER_ATS(dim_privil,privilegios,shell)
 
  IMPLICIT NONE
 
  INTEGER,INTENT(IN)::dim_privil
  INTEGER,DIMENSION(dim_privil),INTENT(IN)::privilegios
  TYPE(p_shell),TARGET,INTENT(INOUT)::shell
 
  INTEGER::ii,jj,kk,gg,priv,ll
  INTEGER::numats,num_aux
  TYPE(p_bonded),POINTER::bonded
  TYPE(p_at),POINTER::at
  INTEGER,DIMENSION(:,:),ALLOCATABLE::valores
  INTEGER,DIMENSION(:),ALLOCATABLE::val_aux
 
  numats=shell%nats

  ALLOCATE(valores(numats,dim_privil*2))
  valores(:,:)=0

  gg=0
  DO ii=1,dim_privil
     priv=privilegios(ii)
     DO jj=1,numats
        DO ll=1,2
           bonded=>shell%ats(jj)%bonded(ll)
           DO kk=1,bonded%num
              IF (priv==bonded%nods(kk)) THEN
                 valores(jj,gg+ll)=valores(jj,gg+ll)+1
              END IF
           END DO
        END DO
     END DO
     gg=gg+2
  END DO

  CALL SORT_INT_MATRIX (numats,gg,shell%order,valores)
  
  DEALLOCATE(valores)

  IF (symm%num_crits>0) THEN

     ALLOCATE(valores(numats,2))
     DO jj=1,numats
        at=>shell%ats(jj)
        DO ll=1,2
           valores(jj,ll)=at%bonded(ll)%num
        END DO
     END DO
 
     CALL SORT_INT_MATRIX (numats,2,shell%order,valores)

     DEALLOCATE(valores)

  END IF
 
  DO ll=1,2

     IF (symm%num_crits==0) EXIT

     gg=0
     DO ii=1,numats
        IF (gg<shell%ats(ii)%bonded(ll)%num) gg=shell%ats(ii)%bonded(ll)%num
     END DO
     ALLOCATE(valores(numats,gg))
     valores(:,:)=0
     DO ii=1,numats
        bonded=>shell%ats(ii)%bonded(ll)
        num_aux=bonded%num
        ALLOCATE(val_aux(num_aux))
        CALL ordered_version(bonded%order,bonded%lev_supsets,val_aux,num_aux)
        valores(ii,1:num_aux)=val_aux(:)
        DEALLOCATE(val_aux)
     END DO
     CALL SORT_INT_MATRIX (numats,gg,shell%order,valores)

     DEALLOCATE(valores)

     IF (symm%num_crits==0) EXIT

     ALLOCATE(valores(numats,gg))
     valores(:,:)=0

     DO ii=1,numats
        bonded=>shell%ats(ii)%bonded(ll)
        num_aux=bonded%num
        ALLOCATE(val_aux(num_aux))
        CALL ordered_version(bonded%order,bonded%loops_in,val_aux,num_aux)
        valores(ii,1:num_aux)=val_aux(:)
        DEALLOCATE(val_aux)
     END DO

     CALL SORT_INT_MATRIX (numats,gg,shell%order,valores)

     DEALLOCATE(valores)

  END DO

  NULLIFY(bonded,at)
 
END SUBROUTINE ORDER_ATS


!!!###################################################################################################
!!!####  REPETITION AND ORDER FUNCTIONS
!!!###################################################################################################

SUBROUTINE ORDERED_VERSION(order,original,ordered,dim)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::dim
  INTEGER,DIMENSION(dim),INTENT(IN)::order,original
  INTEGER,DIMENSION(dim),INTENT(OUT)::ordered
  
  INTEGER::ii

  DO ii=1,dim
     ordered(ii)=original(order(ii))
  END DO

END SUBROUTINE ORDERED_VERSION
!############################

SUBROUTINE BUILDING_REPE_LISTS (carro,dim_carro,nods_norepe,dim_norepe,nods_repe,cant_repe,dim_repe)
 
  IMPLICIT NONE
  INTEGER,INTENT(IN)::dim_carro
  INTEGER,DIMENSION(dim_carro),INTENT(IN)::carro
  INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(INOUT)::nods_norepe,nods_repe,cant_repe
  INTEGER,INTENT(OUT)::dim_norepe,dim_repe
 
  LOGICAL,DIMENSION(dim_carro)::filtro,filtro2
  INTEGER,DIMENSION(dim_carro)::box,repetitions,box2
  INTEGER::ii,jj,gg,kk
 
  filtro(:)=.FALSE.
  filtro2(:)=.FALSE.
  repetitions(:)=0
  DO ii=1,dim_carro
     IF (filtro(ii).eqv..FALSE.) THEN
        filtro(ii)=.TRUE.
        gg=carro(ii)
        kk=1
        DO jj=ii+1,dim_carro
           IF (filtro(jj).eqv..FALSE.) THEN
              IF (carro(jj)==gg) THEN
                 filtro(jj)=.TRUE.
                 kk=kk+1
              END IF
           END IF
        END DO
        filtro2(ii)=.TRUE.
        repetitions(ii)=kk
     END IF
  END DO
 
  dim_norepe=COUNT(filtro2(:),DIM=1)
 
  ALLOCATE(nods_norepe(dim_norepe))
  IF (dim_norepe<dim_carro) THEN
     jj=0
     kk=0
     DO ii=1,dim_carro
        IF (filtro2(ii).eqv..TRUE.) THEN
           jj=jj+1
           nods_norepe(jj)=carro(ii)
           IF (repetitions(ii)>1) THEN
              kk=kk+1
              box(kk)=carro(ii)
              box2(kk)=repetitions(ii)
           END IF
        END IF
     END DO
     ALLOCATE(nods_repe(kk),cant_repe(kk))
     dim_repe=kk
     nods_repe(:)=box(1:kk)
     cant_repe(:)=box2(1:kk)
  ELSE
     dim_repe=0
     nods_norepe(:)=carro(:)
     ALLOCATE(nods_repe(0),cant_repe(0))
  END IF
 
END SUBROUTINE BUILDING_REPE_LISTS
!############################

SUBROUTINE BUILDING_REPE_LIST2 (core,carro,dim_carro,cant_repe)
 
  IMPLICIT NONE
  INTEGER,INTENT(IN)::dim_carro,core
  INTEGER,DIMENSION(dim_carro),INTENT(IN)::carro
  INTEGER,DIMENSION(dim_carro),INTENT(OUT)::cant_repe
  INTEGER,DIMENSION(dim_carro)::aux
  LOGICAL,DIMENSION(dim_carro)::filtro
  INTEGER::ii,jj,gg,kk

  cant_repe(:)=0

  WHERE (carro==core)
     filtro=.TRUE.
  ELSEWHERE
     filtro=.FALSE.
  END WHERE

  DO ii=1,dim_carro
     IF (filtro(ii).eqv..FALSE.) THEN
        filtro(ii)=.TRUE.
        gg=carro(ii)
        kk=1
        aux(kk)=gg
        DO jj=ii+1,dim_carro
           IF (filtro(jj).eqv..FALSE.) THEN
              IF (carro(jj)==gg) THEN
                 filtro(jj)=.TRUE.
                 kk=kk+1
                 aux(kk)=jj
              END IF
           END IF
        END DO
        IF (kk>1) THEN
           DO jj=1,kk
              cant_repe(aux(jj))=kk
           END DO
        END IF
     END IF
  END DO
 
END SUBROUTINE BUILDING_REPE_LIST2
!############################

SUBROUTINE BUILDING_REPE_LIST3 (core,carro1,carro2,dim_carro1,dim_carro2,cant_repe1,cant_repe2)
 
  IMPLICIT NONE
  INTEGER,INTENT(IN)::core,dim_carro1,dim_carro2
  INTEGER,DIMENSION(dim_carro1),INTENT(IN)::carro1
  INTEGER,DIMENSION(dim_carro2),INTENT(IN)::carro2
  INTEGER,DIMENSION(dim_carro1),INTENT(OUT)::cant_repe1
  INTEGER,DIMENSION(dim_carro2),INTENT(OUT)::cant_repe2

  INTEGER,DIMENSION(dim_carro2)::aux
  LOGICAL,DIMENSION(dim_carro2)::filtro

  INTEGER::ii,jj,gg,kk

  cant_repe1(:)=0
  cant_repe2(:)=0

  WHERE (carro2==core)
     filtro=.TRUE.
  ELSEWHERE
     filtro=.FALSE.
  END WHERE

  DO ii=1,dim_carro1
     gg=carro1(ii)
     kk=0
     DO jj=1,dim_carro2
        IF (filtro(ii).eqv..FALSE.) THEN
           IF (carro2(jj)==gg) THEN
              filtro(jj)=.TRUE.
              kk=kk+1
              aux(kk)=jj
           END IF
        END IF
     END DO
     IF (kk>0) THEN
        DO jj=1,kk
           cant_repe2(aux(jj))=kk
        END DO
        cant_repe1(ii)=kk
     END IF
  END DO

END SUBROUTINE BUILDING_REPE_LIST3

!!!###################################################################################################
!!!####  SORTING FUNCTIONS
!!!###################################################################################################

SUBROUTINE SORT_INT_MATRIX (num_ats,dim_matrix,order,valores)
 
  IMPLICIT NONE
 
  INTEGER,INTENT(IN)::num_ats,dim_matrix
  INTEGER,DIMENSION(num_ats),INTENT(INOUT)::order
  INTEGER,DIMENSION(num_ats,dim_matrix),INTENT(IN)::valores
 
  INTEGER::pp,ii,jj,kk,gg,ll,mm,idim,tope,new_num_crits
  INTEGER,DIMENSION(:),ALLOCATABLE::val_aux,ind_aux,symm_aux
  INTEGER,DIMENSION(:),ALLOCATABLE::vals,inds,new_symm_aux,cajon
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro
  LOGICAL::interruptor
 
  pp=1
  DO WHILE ((symm%num_crits>0).AND.(pp<=dim_matrix))
     new_num_crits=0
     gg=0
     DO ii=1,symm%num_crits
        gg=gg+1
        idim=symm%crit(gg)
        ALLOCATE(val_aux(idim),ind_aux(idim),symm_aux(idim))
        DO jj=1,idim
           gg=gg+1
           kk=symm%crit(gg)
           ll=order(kk)
           val_aux(jj)=valores(ll,pp)
           symm_aux(jj)=kk
           ind_aux(jj)=ll
        END DO
        IF (ALL(val_aux(2:idim).eq.val_aux(1)).eqv..FALSE.) THEN
           ALLOCATE(filtro(idim),vals(idim),inds(idim))
           filtro(1:idim)=.TRUE.
           DO jj=1,idim
              kk=MAXLOC(val_aux(1:idim),DIM=1,MASK=filtro(1:idim))
              filtro(kk)=.FALSE.
              ll=symm_aux(jj)
              order(ll)=ind_aux(kk)
              inds(jj)=ll
              vals(jj)=val_aux(kk)
           END DO
           interruptor=.FALSE.
           mm=0
           DO jj=2,idim
              IF (vals(jj-1)==vals(jj)) THEN
                 IF (interruptor.eqv..FALSE.) THEN
                    ll=2
                    interruptor=.TRUE.
                    mm=jj-1
                 ELSE
                    ll=ll+1
                 END IF
              ELSE
                 IF (interruptor.eqv..True.) THEN
                    interruptor=.FALSE.
                    new_num_crits=new_num_crits+1
                    IF (new_num_crits==1) THEN
                       tope=ll+1
                       ALLOCATE(new_symm_aux(tope))
                       new_symm_aux(1)=ll
                       new_symm_aux(2:tope)=inds(mm:(mm+ll-1))
                    ELSE
                       ALLOCATE(cajon(tope))
                       cajon(1:tope)=new_symm_aux(1:tope)
                       DEALLOCATE(new_symm_aux)
                       ALLOCATE(new_symm_aux(tope+1+ll))
                       new_symm_aux(1:tope)=cajon(1:tope)
                       tope=tope+1
                       new_symm_aux(tope)=ll
                       new_symm_aux((tope+1):(tope+ll))=inds(mm:(mm+ll-1))
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
                 new_symm_aux(2:tope)=inds(mm:(mm+ll-1))
              ELSE
                 ALLOCATE(cajon(tope))
                 cajon(1:tope)=new_symm_aux(1:tope)
                 DEALLOCATE(new_symm_aux)
                 ALLOCATE(new_symm_aux(tope+1+ll))
                 new_symm_aux(1:tope)=cajon(1:tope)
                 tope=tope+1
                 new_symm_aux(tope)=ll
                 new_symm_aux((tope+1):(tope+ll))=inds(mm:(mm+ll-1))
                 tope=tope+ll
                 DEALLOCATE(cajon)
              END IF
           END IF
           DEALLOCATE(filtro,vals,inds)
        ELSE
           new_num_crits=new_num_crits+1
           IF (new_num_crits==1) THEN
              tope=idim+1
              ALLOCATE(new_symm_aux(tope))
              new_symm_aux(1)=idim
              new_symm_aux(2:tope)=symm_aux(1:idim)
           ELSE
              ALLOCATE(cajon(tope))
              cajon(1:tope)=new_symm_aux(1:tope)
              DEALLOCATE(new_symm_aux)
              ALLOCATE(new_symm_aux(tope+1+idim))
              new_symm_aux(1:tope)=cajon(1:tope)
              tope=tope+1
              new_symm_aux(tope)=idim
              new_symm_aux((tope+1):(tope+idim))=symm_aux(1:idim)
              tope=tope+idim
              DEALLOCATE(cajon)
           END IF
        END IF
        DEALLOCATE(val_aux,ind_aux,symm_aux)
     END DO
     symm%num_crits=new_num_crits
     IF (symm%num_crits>0) THEN
        DEALLOCATE(symm%crit)
        ALLOCATE(symm%crit(tope))
        symm%length=tope
        symm%crit(:)=new_symm_aux(:)
        DEALLOCATE(new_symm_aux)
     END IF
     pp=pp+1
  END DO
 
END SUBROUTINE SORT_INT_MATRIX


!!!###################################################################################################
!!!####  AUXILIAR FUNCTIONS TO BUILD MSS
!!!###################################################################################################


SUBROUTINE BUILD_MSS_ATS(shell)

  IMPLICIT NONE
 
  TYPE(p_shell),TARGET,INTENT(INOUT)::shell

  INTEGER::ntot,nats,nats2,nnods
  INTEGER::ii,jj,kk,ll,num_aux
  INTEGER,DIMENSION(:),ALLOCATABLE::oo,val_aux

  TYPE(p_at),POINTER::at
  TYPE(p_bonded),POINTER::bonded

  IF (ALLOCATED(mss_ind_ats)) THEN
     DEALLOCATE(mss_ind_ats,mss_template)
  END IF

  ntot = shell%ntot
  nats = shell%nats
  nats2 = shell%nats2
  nnods = shell%nnods
 
  ALLOCATE(mss_ind_ats(ntot),mss_template(ntot))
  ALLOCATE(oo(nats))

  mss_template(:)=.FALSE.

  oo(:)=shell%order(:)

  jj=1
  mss_ind_ats(1)  = nats
  ii=jj+1
  jj=jj+nats
  mss_ind_ats(ii:jj)  = (/(shell%ats(oo(kk))%ind,kk=1,nats)/)
  mss_template(ii:jj)  = .TRUE.
  ii=jj+1
  jj=jj+nats2
  mss_ind_ats(ii:jj)  = (/((/shell%ats(oo(kk))%bonded(1)%num,shell%ats(oo(kk))%bonded(2)%num/),kk=1,nats)/)
  DO kk=1,nats
     at=>shell%ats(oo(kk))
     IF (at%empty.eqv..FALSE.) THEN
        DO ll=1,2
           IF (at%bonded(ll)%num>0) THEN
              bonded=>at%bonded(ll)
              num_aux=bonded%num
              ii=jj+1
              jj=jj+num_aux
              ALLOCATE(val_aux(num_aux))
              CALL ordered_version(bonded%order,bonded%ats,val_aux,num_aux)
              mss_ind_ats(ii:jj) = (/val_aux(:)/)
              mss_template(ii:jj) = .TRUE.
              DEALLOCATE(val_aux)
           END IF
        END DO
     END IF
  END DO

  DEALLOCATE(oo)

  NULLIFY(at,bonded)

END SUBROUTINE BUILD_MSS_ATS
!############################

SUBROUTINE BUILD_MSS_NODS()

  IMPLICIT NONE
 
  IF (ALLOCATED(mss_ind_nods))   DEALLOCATE(mss_ind_nods)
  ALLOCATE(mss_ind_nods(SIZE(mss_ind_ats)))

  WHERE (mss_template.eqv..TRUE.)
     mss_ind_nods=at2nod(mss_ind_ats)
  ELSEWHERE
     mss_ind_nods=mss_ind_ats
  END WHERE

END SUBROUTINE BUILD_MSS_NODS
!############################

SUBROUTINE BUILD_MSS_SYMM(shell)

  IMPLICIT NONE
 
  TYPE(p_shell),TARGET,INTENT(INOUT)::shell
 
  TYPE(p_at),POINTER::at
  TYPE(p_bonded),POINTER::bonded
 
  INTEGER::ntot,nats,nats2,nnods
  INTEGER::ii,jj,kk,ll,mm,nn,pp
  INTEGER,DIMENSION(:),ALLOCATABLE::oo

  IF (ALLOCATED(mss_symm))      DEALLOCATE(mss_symm)

  ntot = shell%ntot
  nats = shell%nats
  nats2 = shell%nats2
  nnods = shell%nnods
  ALLOCATE(oo(nats))
  oo(:)=shell%order(:)
  ALLOCATE(mss_symm(ntot))

  DO ii=1,nats
     at=>shell%ats(ii)
     DO ll=1,2
        IF (at%bonded(ll)%symm%num_crits>0) THEN
           bonded=>at%bonded(ll)
           mm=0
           DO nn=1,bonded%symm%num_crits
              mm=mm+1
              DO pp=1,bonded%symm%crit(mm)
                 mm=mm+1
                 bonded%wsymm(bonded%symm%crit(mm))=1
              END DO
           END DO
        END IF
     END DO
  END DO

  IF (shell%symm%num_crits>0) THEN
     mm=0
     DO nn=1,shell%symm%num_crits
        mm=mm+1
        DO ll=1,shell%symm%crit(mm)
           mm=mm+1
           shell%wsymm(shell%symm%crit(mm))=1
        END DO
     END DO
  END IF
 
  jj=1
  mss_symm(1)         = 0
  ii=jj+1
  jj=jj+nats
  mss_symm(ii:jj)     = (/(shell%wsymm(kk),kk=1,nats)/)
  ii=jj+1
  jj=jj+nats2
  mss_symm(ii:jj)     = 0
  ii=jj+1
  jj=jj+nnods
  mss_symm(ii:jj)     = (/((/shell%ats(oo(kk))%bonded(1)%wsymm(:),shell%ats(oo(kk))%bonded(2)%wsymm(:)/),kk=1,nats)/)


  DEALLOCATE(oo)
  NULLIFY(at,bonded)

END SUBROUTINE BUILD_MSS_SYMM
!############################

SUBROUTINE BUILD_MSS ()

  IMPLICIT NONE

  INTEGER::ii,dim

  CALL RESET_TRANSLATOR()

  dim=SIZE(mss_ind_nods)

  IF (ALLOCATED(mss))  DEALLOCATE(mss)
  ALLOCATE(mss(dim))
 
  DO ii=1,dim
     IF (mss_template(ii).eqv..TRUE.) THEN
        mss(ii)=TRANS_WO_SUPSETS(mss_ind_nods(ii))
     ELSE
        mss(ii)= mss_ind_nods(ii)
     END IF
  END DO

END SUBROUTINE BUILD_MSS
!############################

SUBROUTINE BUILD_TRANS_MSS_ATS()
 
  WHERE (mss_template.eqv..TRUE.)
     mss_ind_ats=trad2py_at(mss_ind_ats)
  ELSEWHERE
     mss_ind_ats=mss_ind_ats
  END WHERE

END SUBROUTINE BUILD_TRANS_MSS_ATS
!############################

SUBROUTINE BUILD_TRANS_MSS_NODS()

  WHERE (mss_template.eqv..TRUE.)
     mss_ind_nods=trad2py_nod(mss_ind_nods)
  ELSEWHERE
     mss_ind_nods=mss_ind_nods
  END WHERE

END SUBROUTINE BUILD_TRANS_MSS_NODS
!############################

INTEGER FUNCTION TRANS_WO_SUPSETS(in_node)
 
  IMPLICIT NONE
 
  INTEGER,INTENT(IN)::in_node
 
  INTEGER::ii_supset,aa,jj
  TYPE(p_trans_set),POINTER::supset_aux

  IF (translator%total_filtro(in_node).eqv..FALSE.) THEN

     translator%total_filtro(in_node)=.TRUE.
 
     ii_supset=lev_supsets(in_node)
     supset_aux=>translator%supset(ii_supset)
     
     aa=supset_aux%count_sets+1
     supset_aux%count_sets=aa
     jj=supset_aux%tope+aa
     translator%dict(in_node)=jj
     TRANS_WO_SUPSETS=jj

     NULLIFY(supset_aux)
 
  ELSE

     TRANS_WO_SUPSETS=translator%dict(in_node)

  END IF

END FUNCTION TRANS_WO_SUPSETS
 

END MODULE GLOB


