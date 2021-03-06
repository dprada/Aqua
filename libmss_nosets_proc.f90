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
   INTEGER,DIMENSION(:),ALLOCATABLE::cant_loops_1order,cant_loops_2order,cant_loops_3order
   LOGICAL::with_loops_1order,with_loops_2order,with_loops_3order
   INTEGER::num
   LOGICAL::unsolved_loops
   INTEGER::num_permutations
   INTEGER,DIMENSION(:,:),ALLOCATABLE::permutations
   TYPE(p_symm)::symm
END TYPE p_bonded

TYPE p_at
   TYPE(p_bonded),DIMENSION(2)::bonded
   INTEGER::num_hbs_bs
   INTEGER::ind
END TYPE p_at

TYPE p_shell
   TYPE(p_at),DIMENSION(:),ALLOCATABLE::ats
   INTEGER,DIMENSION(:),ALLOCATABLE::ats_ord
   INTEGER::nats,nats2,nnods,ntot
   INTEGER::ind
   INTEGER,DIMENSION(:),ALLOCATABLE::wsymm
   TYPE(p_symm)::symm
   INTEGER::lev_supsets
   INTEGER::nnods_norepe,nnods_repe
   LOGICAL::unsolved_loops
   INTEGER::num_permutations
   INTEGER,DIMENSION(:,:),ALLOCATABLE::permutations
   INTEGER,DIMENSION(:),ALLOCATABLE::list_nods,list_nods_norepe,list_nods_repe,list_cant_repe
END TYPE p_shell

TYPE p_trans_set
   INTEGER::tope,num_nods
   INTEGER::count_sets
END type p_trans_set

TYPE p_translation
   TYPE(p_trans_set),DIMENSION(:),ALLOCATABLE::supset
   INTEGER::num_supsets
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

!f2py intent(hide)::precalc_2, precalc2_num
!f2py intent(hide)::precalc_3, precalc3_num
!f2py intent(hide)::precalc_4, precalc4_num
INTEGER,PARAMETER::precalc2_num=2
INTEGER,PARAMETER::precalc3_num=6
INTEGER,PARAMETER::precalc4_num=24
INTEGER,DIMENSION(2,2),PARAMETER::&
       precalc2=(/(/1,2/),(/2,1/)/)
INTEGER,DIMENSION(3,6),PARAMETER::&
       precalc3=(/(/1,2,3/),(/1,3,2/),(/2,1,3/),&
       (/2,3,1/),(/3,1,2/),(/3,2,1/)/)
INTEGER,DIMENSION(4,24),PARAMETER::&
       precalc4=(/(/1,2,3,4/),(/1,2,4,3/),(/1,3,2,4/),(/1,3,4,2/),(/1,4,2,3/),(/1,4,3,2/),&
       (/2,1,3,4/),(/2,1,4,3/),(/2,3,1,4/),(/2,3,4,1/),(/2,4,1,3/),(/2,4,3,1/),&
       (/3,1,2,4/),(/3,1,4,2/),(/3,2,1,4/),(/3,2,4,1/),(/3,4,1,2/),(/3,4,2,1/),&
       (/4,1,2,3/),(/4,1,3,2/),(/4,2,1,3/),(/4,2,3,1/),(/4,3,1,2/),(/4,3,2,1/)/)


!!! OUTPUT

INTEGER,DIMENSION(:),ALLOCATABLE::mss_ind_ats,mss_ind_nods,mss,mss_symm


CONTAINS

!!!###################################################################################################
!!!####  AUXILIAR TOPOL AND NET SUBROUTINES
!!!###################################################################################################

SUBROUTINE ATS_UP(at)
 
  IMPLICIT NONE
 
  TYPE(p_at),TARGET,INTENT(INOUT)::at
  TYPE(p_bonded),POINTER::bonded
  INTEGER::ii,kk,num

  DO ii=1,2
     bonded=>at%bonded(ii)
     num=bonded%num
     ALLOCATE(bonded%ats(num))
     ALLOCATE(bonded%nods(num))
     ALLOCATE(bonded%order(num))
     ALLOCATE(bonded%lev_supsets(num))
     ALLOCATE(bonded%cant_loops_1order(num))
     ALLOCATE(bonded%cant_loops_2order(num))
     ALLOCATE(bonded%cant_loops_3order(num))
     ALLOCATE(bonded%aux2sh(num)) 
     ALLOCATE(bonded%wsymm(num)) 
     ALLOCATE(bonded%symm%crit(num+1)) !!!
     bonded%order                =(/(kk,kk=1,num)/)
     bonded%cant_loops_1order(:) = 0
     bonded%cant_loops_2order(:) = 0
     bonded%cant_loops_3order(:) = 0
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
     bonded%unsolved_loops       =.FALSE.
     bonded%with_loops_1order    =.FALSE.
     bonded%with_loops_2order    =.FALSE.
     bonded%with_loops_3order    =.FALSE.
  END DO

  NULLIFY(bonded)

END SUBROUTINE ATS_UP

SUBROUTINE ATS_DOWN(at)
 
  IMPLICIT NONE
 
  TYPE(p_at),TARGET,INTENT(INOUT)::at
  TYPE(p_bonded),POINTER::bonded
  INTEGER::ii

  DO ii=1,2
     bonded=>at%bonded(ii)
     DEALLOCATE(bonded%ats)
     DEALLOCATE(bonded%nods)
     DEALLOCATE(bonded%order)
     DEALLOCATE(bonded%lev_supsets)
     DEALLOCATE(bonded%cant_loops_1order)
     DEALLOCATE(bonded%cant_loops_2order)
     DEALLOCATE(bonded%cant_loops_3order)
     DEALLOCATE(bonded%aux2sh)
     DEALLOCATE(bonded%wsymm)
     DEALLOCATE(bonded%symm%crit)
  END DO

  NULLIFY(bonded)

END SUBROUTINE ATS_DOWN

SUBROUTINE UPLOAD_SYMM(loc_symm)
 
  IMPLICIT NONE
 
  TYPE(p_symm),INTENT(IN)::loc_symm
 
  symm%num_crits=loc_symm%num_crits
  symm%length=loc_symm%length
 
  ALLOCATE(symm%crit(symm%length))
  symm%crit(:)=loc_symm%crit(:)
 
END SUBROUTINE UPLOAD_SYMM

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


SUBROUTINE RESET_TRANSLATOR()
 
  IMPLICIT NONE
 
  INTEGER::ii
 
  DO ii=1,translator%num_supsets
     translator%supset(ii)%count_sets=0
  END DO
 
END SUBROUTINE RESET_TRANSLATOR


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
  TYPE(p_shell),POINTER::shell_aux


  num_ats     = xx_num_ats
  num_nods    = xx_num_nods


  IF (ALLOCATED(at2nod))           DEALLOCATE(at2nod)
  IF (ALLOCATED(trad2py_at))       DEALLOCATE(trad2py_at)
  IF (ALLOCATED(trad2py_nod))      DEALLOCATE(trad2py_nod)
  IF (ALLOCATED(list_shells))      DEALLOCATE(list_shells) !! cuidado aqui

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
     shell_aux=>list_shells(ii)
     shell_aux%ind=ii
     shell_aux%nats=kk
     shell_aux%nats2=kk*2
     shell_aux%nnods=0
     shell_aux%ntot=0
     ALLOCATE(shell_aux%ats_ord(kk))
     ALLOCATE(shell_aux%ats(kk))
     ALLOCATE(shell_aux%wsymm(kk))
     shell_aux%ats_ord(:)=(/(jj,jj=1,kk)/)
     shell_aux%wsymm(:)=0
     DO jj=1,kk
        gg=gg+1
        shell_aux%ats_ord(jj)=jj
        shell_aux%ats(jj)%ind=gg
     END DO
     jj=xx_symm_ats_crits(ii)
     shell_aux%symm%num_crits=jj
     kk=xx_symm_ats_start(ii)-1
     mm=0
     DO ll=1,jj
        mm=mm+1
        mm=mm+xx_symm_ats(kk+mm)
     END DO
     shell_aux%symm%length=mm
     ALLOCATE(shell_aux%symm%crit(mm))
     shell_aux%symm%crit(:)=xx_symm_ats((kk+1):(kk+mm))
     ALLOCATE(shell_aux%list_nods_norepe(0),shell_aux%list_nods_repe(0),shell_aux%list_cant_repe(0),&
          shell_aux%list_nods(0))
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

  NULLIFY(shell_aux)
  DEALLOCATE(num_ats_nod)

END SUBROUTINE LOAD_TOPOL





SUBROUTINE load_net(xx_hbs,xx_bs,xx_num_Hbs_at,xx_num_Bs_at,&
     xx_Total_num_hbs,xx_Total_num_bs,xx_num_ats)
 
  IMPLICIT NONE
 
  INTEGER,INTENT(IN)::xx_num_ats,xx_Total_num_hbs,xx_Total_num_bs
  INTEGER,DIMENSION(xx_num_ats),INTENT(IN)::xx_num_hbs_at,xx_num_bs_at
  INTEGER,DIMENSION(xx_Total_num_hbs),INTENT(IN)::xx_hbs
  INTEGER,DIMENSION(xx_Total_num_bs),INTENT(IN)::xx_bs
 
  INTEGER::ii,jj,kk,gghb,ggb,ggat,numhbs,numbs,totnum
  INTEGER,DIMENSION(:),ALLOCATABLE::vect_aux_hbs,vect_aux_bs
  TYPE(p_shell),POINTER::shell
  TYPE(p_bonded),POINTER::bonded
  TYPE(p_at),POINTER::at
 
  gghb=0
  ggb=0
  DO ii=1,num_nods
     shell=>list_shells(ii)
     shell%unsolved_loops=.FALSE.
     shell%nnods=0
     DEALLOCATE(shell%list_nods_norepe)
     DEALLOCATE(shell%list_nods_repe)
     DEALLOCATE(shell%list_cant_repe)
     DEALLOCATE(shell%list_nods)
     DO jj=1,shell%nats
        at=>shell%ats(jj)
        CALL ATS_DOWN(at)
        ggat=at%ind
        numhbs=xx_num_hbs_at(ggat)
        numbs=xx_num_bs_at(ggat)
        at%bonded(1)%num=numhbs
        at%bonded(2)%num=numbs
        CALL ATS_UP(at)
        !HBS
        kk=gghb+1
        gghb=gghb+num
        at%bonded(1)%ats(:)=xx_hbs(kk:gghb)
        !BS
        kk=ggb+1
        ggb=ggb+num
        at%bonded(2)%ats(:)=xx_bs(kk:ggb)
        !! RECOLOCO
        DO ll=1,2
           bonded=>at%bonded(ll)
           at%num_hbs_bs=at%num_hbs_bs+bonded%num
           shell%nnods=shell%nnods+bonded%num
           ALLOCATE(vect_aux(bonded%num))
           vect_aux(:)                 = at2nod(bonded%ats(:))
           bonded%nods(:)              = vect_aux(:)
           bonded%lev_supsets(:)       = lev_supsets(vect_aux(:))
           DEALLOCATE(vect_aux)
        END DO
     END DO
     ALLOCATE(shell%list_nods(shell%nnods))
     shell%list_nods(:)=(/((/shell%ats(kk)%bonded(1)%nods(:),shell%ats(kk)%bonded(2)%nods(:)/),kk=1,shell%nats)/)
     CALL BUILDING_REPE_LIST (shell%list_nods,shell%nnods,shell%list_nods_norepe,&
          shell%nnods_norepe,shell%list_nods_repe,shell%list_cant_repe,shell%nnods_repe)
     shell%ntot=1+shell%nats+shell%nats2+shell%nnods
     IF (shell%nnods_repe>0) THEN
        DO jj=1,shell%nats
           at=>shell%ats(jj)
           DO ll=1,2
              IF (at%bonded(ll)%num>0) THEN
                 CALL COUNT_NODE_REPE (shell%list_nods_repe,shell%list_cant_repe,shell%nnods_repe,&
                      at%bonded(ll)%nods,at%bonded(ll)%num,at%bonded(ll)%cant_loops_1order)
                 IF (SUM(at%bonded(ll)%cant_loops_1order(:))>0) at%bonded(ll)%with_loops_1order=.TRUE.
              END IF
           END DO
        END DO
     END IF
  END DO
 
  NULLIFY(shell,bonded,at)
  
END SUBROUTINE load_net
 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE build_shell1st (core)
 
  IMPLICIT NONE
 
  INTEGER,INTENT(IN)::core
  TYPE(p_shell),TARGET::shell
  TYPE(p_at),POINTER::at
  TYPE(p_bonded),POINTER::bonded
 
  INTEGER,DIMENSION(:),ALLOCATABLE::privilegios
  INTEGER::dim_privil,ii
 
 
  shell1st=list_shells(core)
   
  !! order bonds
   
  dim_privil=1
  ALLOCATE(privilegios(1))
  privilegios(1)=core
 
  DO ii=1,shell1st%nats
     at=>shell1st%ats(ii)
     DO jj=1,2
        IF (at%bonded(jj)%symm%num_crits>0) THEN
           bonded=>at%bonded(jj)
           CALL upload_symm(bonded%symm)
           CALL order_bonded(dim_privil,privilegios,bonded)
           CALL download_symm(bonded%symm)
        END IF
     END DO
  END DO
 
  !! order ats

  IF (shell1st%symm%num_crits>0) THEN
     CALL upload_symm(shell1st%symm)
     CALL remove_triv_symm_ats(shell1st)
     IF (symm%num_crits>0) THEN
        CALL order_ats(dim_privil,privilegios,shell1st)
     END IF
     CALL download_symm(shell1st%symm)
  END IF
  
  DEALLOCATE(privilegios)

  CALL build_mss(shell1st)
  CALL remato_mss()
   
  NULLIFY(at,bonded,shell1st)
 
END SUBROUTINE build_shell1st

 
SUBROUTINE build_shell2nd (core)
 
  IMPLICIT NONE
 
  INTEGER,INTENT(IN)::core
 
  TYPE(p_shell),TARGET::shell1st
  TYPE(p_shell),DIMENSION(:),ALLOCATABLE,TARGET::shell2nd
 
  TYPE(p_bonded),POINTER::bs_aux
  TYPE(p_at),POINTER::at_aux
  TYPE(p_shell),POINTER::shell_aux
 
  INTEGER::nnods,totntot,core2
  INTEGER::ntot1sh,nnods1sh,ntot2sh,nnods2sh
  INTEGER::ii,jj,kk,gg,iii,ggg,ll
  INTEGER,DIMENSION(:),ALLOCATABLE::privilegios
  INTEGER::dim_privil,dim_repe
  INTEGER,DIMENSION(:),ALLOCATABLE::aux_ind_ats,aux_ind_nods,aux_ord,aux_symm
  INTEGER,DIMENSION(:),ALLOCATABLE::carro_loops3,list_repe,cant_repe
 
 
  !! Building the structure
 
  shell1st=list_shells(core)
 
  nnods=shell1st%nnods
  totntot=shell1st%ntot
 
  ALLOCATE(shell2nd(shell1st%nnods))
  gg=0
  DO ii=1,shell1st%nats
     DO ll=1,2
        DO jj=1,shell1st%ats(ii)%bonded(ll)%num
           gg=gg+1
           shell1st%ats(ii)%bonded(ll)%aux2sh(jj)=gg
           kk=shell1st%ats(ii)%bonded(ll)%bonded_nods(jj)
           shell2nd(gg)=list_shells(kk)
           totntot=totntot+shell2nd(gg)%ntot
        END DO
     END DO
  END DO
 
  !! Building cant_loops_2order(:)
 
  DO iii=1,nnods
     shell_aux=>shell2nd(iii)
     DO ii=1,shell_aux%nats
        at_aux=>shell_aux%ats(ii)
        DO ll=1,2
           bonded=>at%bonded(ll)
           DO jj=1,bonded%num
              gg=COUNT((shell1st%list_nods(:)==bonded%nods(jj)),DIM=1)
              bonded%cant_loops_2order(jj)=gg
              IF (gg>1) bonded%with_loops_2order=.TRUE.
           END DO
        END DO
     END DO
  END DO


  !! Building cant_loops_3order(:) !! hay que quitar el core
 
  gg=0
  DO iii=1,nnods
     gg=gg+shell2nd(iii)%nnods_norepe
  END DO
  ALLOCATE(carro_loops3(gg))
  ggg=0
  DO iii=1,nnods
     gg=ggg+1
     ggg=ggg+shell2nd(iii)%nnods_norepe
     carro_loops3(gg:ggg)=shell2nd(iii)%list_nods_norepe(:)
  END DO
 
  CALL BUILDING_REPE_LIST2 (core,carro_loops3,ggg,list_repe,cant_repe,dim_repe) !quito el core

  IF (dim_repe>0) THEN
     DO iii=1,nnods
        shell_aux=>shell2nd(iii)
        DO ii=1,shell_aux%nats
           at=>shell_aux%ats(ii)
           DO ll=1,2
              bonded=>at%bonded(ll)
              DO jj=1,bonded%num
                 DO ggg=1,dim_repe
                    IF (bonded%nods(jj)==list_repe(ggg)) THEN
                       bonded%cant_loops_3order(jj)=cant_repe(ggg)
                       IF (cant_repe(ggg)>0) bonded%with_loops_3order=.TRUE.
                    END IF
                 END DO
              END DO
           END DO
        END DO
     END DO
     shell_aux=>shell1st
     DO ii=1,shell_aux%nats
        at=>shell_aux%ats(ii)
        DO ll=1,2
           bonded=>at%bonded(ll)
           DO jj=1,bonded%num
              DO ggg=1,dim_repe
                 IF (bonded%nods(jj)==list_repe(ggg)) THEN
                    bonded%cant_loops_2order(jj)=cant_repe(ggg)
                    IF (cant_repe(ggg)>0) bonded%with_loops_2order=.TRUE.
                 END IF
              END DO
           END DO
        END DO
     END DO
  END IF
 
  DEALLOCATE(list_repe,cant_repe)
 
  !! removing symm in 2ndsh
 
  dim_privil=2
  ALLOCATE(privilegios(2))
  privilegios(2)=core
 
  DO iii=1,nnods
 
     shell_aux=>shell2nd(iii)
     core2=shell_aux%ind
 
     !list of privilegios
     privilegios(1)=core2 
 
     ! order bonds 2sh
     DO ii=1,shell_aux%nats
        at=>shell_aux%ats(ii)
        DO ll=1,2
           IF (at%bonded(ll)%symm%num_crits>0) THEN
              bonded=>at%bonded(ll)
              CALL upload_symm(bonded%symm)
              CALL order_bonded(dim_privil,privilegios,bonded)
              IF (symm%num_crits>0) THEN
                 CALL order_bonded_w_23order(bonded)
              END IF
              IF (symm%num_crits>0) THEN 
                 !! quito falsas simetrias porque me quede sin criterios para los nodos que aparecen una sola vez
                 CALL quito_falsos_bonded_symm(bonded)
              END IF
              CALL download_symm(bonded%symm)
           END IF
        END DO
     END DO
 
     ! order ats 2sh
     IF (shell_aux%symm%num_crits>0) THEN
        CALL upload_symm(shell_aux%symm)
        CALL remove_triv_symm_ats(shell_aux) ! el que no tiene nada no tiene ni privilegios
        IF (symm%num_crits>0) THEN
           CALL order_ats(dim_privil,privilegios,shell_aux)
        END IF
        IF (symm%num_crits>0) THEN
           CALL order_ats_w_23order(shell_aux)
        END IF
        IF (symm%num_crits>0) THEN
           CALL quito_falsos_ats_symm(shell_aux)
        END IF
        CALL download_symm(shell_aux%symm)
     END IF
 
  END DO
 
  !! removing symm in 1stsh
 
  DEALLOCATE(privilegios)
  dim_privil=1
  ALLOCATE(privilegios(1))
  privilegios(1)=core
 
  ! order bonds 1sh
  DO ii=1,shell1st%nats
     at=>shell1st%ats(ii)
     DO ll=1,2
        IF (at%bonded(ll)%symm%num_crits>0) THEN
           bonded=>at%bonded(ll)
           CALL upload_symm(bonded%symm)
           CALL order_bonded(dim_privil,privilegios,bonded)
           IF (symm%num_crits>0) THEN
              CALL order_bonded_w_2order(bonded)
           END IF
           IF (symm%num_crits>0) THEN
              CALL order_bonded_w_next_shells(bonded,shell2nd,nnods)
           END IF
           IF (symm%num_crits>0) THEN
              CALL quito_falsos_bonded_symm2(bonded,shell2nd,nnods)
           END IF
           CALL download_symm(bonded%symm)
        END IF
     END DO
  END DO
  
  ! order ats 1sh
  IF (shell1st%symm%num_crits>0) THEN
     CALL upload_symm(shell1st%symm)
     CALL remove_triv_symm_ats(shell1st) ! el que no tiene nada no tiene ni privilegios
     IF (symm%num_crits>0) THEN 
        CALL order_ats(dim_privil,privilegios,shell1st)
        IF (symm%num_crits>0) THEN 
           CALL order_ats_w_2order(shell1st)
        END IF
        IF (symm%num_crits>0) THEN 
           CALL order_ats_w_next_shells(shell1st,shell2nd,nnods)
        END IF
        IF (symm%num_crits>0) THEN
           CALL quito_falsos_ats_symm2(shell1st,shell2nd,nnods)
        END IF
     END IF
     CALL download_symm(shell1st%symm)
  END IF
 
  !! Building the msss
  ! Poner aqui un if
  CALL solving_last_symmetries_permuting (shell1st,shell2nd,nnods,totntot)


  IF (ALLOCATED(mss_ind_ats))   DEALLOCATE(mss_ind_ats)
  IF (ALLOCATED(mss_ind_nods))  DEALLOCATE(mss_ind_nods)
  IF (ALLOCATED(mss_symm))      DEALLOCATE(mss_symm)
 
  ALLOCATE(mss_ind_ats(totntot),mss_ind_nods(totntot),mss_symm(totntot))
 
  gg=0
  ntot1sh=shell1st%ntot
  nnods1sh=shell1st%nnods
  ALLOCATE(aux_ind_ats(ntot1sh),aux_ind_nods(ntot1sh),aux_symm(ntot1sh),aux_ord(nnods1sh))
  CALL build_mss_wout_word(shell1st,aux_ind_ats,aux_ind_nods,aux_symm,aux_ord)
  ggg=gg+1
  gg=gg+ntot1sh
  mss_ind_ats(ggg:gg)=aux_ind_ats(:)
  mss_ind_nods(ggg:gg)=aux_ind_nods(:)
  mss_symm(ggg:gg)=aux_symm(:)
  DEALLOCATE(aux_ind_ats,aux_ind_nods,aux_symm)
 
  DO ii=1,nnods1sh
     jj=aux_ord(ii)
     shell_aux=>shell2nd(jj)
     ntot2sh =shell_aux%ntot
     nnods2sh=shell_aux%nnods
     ALLOCATE(aux_ind_ats(ntot2sh),aux_ind_nods(ntot2sh),aux_symm(ntot2sh))
     CALL build_mss_wout(shell_aux,aux_ind_ats,aux_ind_nods,aux_symm)
     ggg=gg+1
     gg=gg+ntot2sh
     mss_ind_ats(ggg:gg)=aux_ind_ats(:)
     mss_ind_nods(ggg:gg)=aux_ind_nods(:)
     mss_symm(ggg:gg)=aux_symm(:)
     DEALLOCATE(aux_ind_ats,aux_ind_nods,aux_symm)
  END DO
 
  CALL remato_mss2sh()
 
  DEALLOCATE(aux_ord,shell2nd,privilegios)
  NULLIFY(at_aux,shell_aux,bs_aux)
 
END SUBROUTINE build_shell2nd
 
 
SUBROUTINE remove_triv_symm_ats(shell)
 
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
     IF (shell%ats(ii)%num_hbs_bs>0) THEN
        valores(ii,1)=1
     ELSE
        valores(ii,1)=0
     END IF
  END DO
 
  CALL SORT_INT_MATRIX (numats,1,shell%ats_ord,valores)
 
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
        ll=shell%ats_ord(symm%crit(gg+1))
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
 
END SUBROUTINE remove_triv_symm_ats
 
 
SUBROUTINE order_ats(dim_privil,privilegios,shell)
 
  IMPLICIT NONE
 
  INTEGER,INTENT(IN)::dim_privil
  INTEGER,DIMENSION(dim_privil),INTENT(IN)::privilegios
  TYPE(p_shell),TARGET,INTENT(INOUT)::shell
 
  INTEGER::ii,jj,kk,gg,priv,ll
  INTEGER::numats
  TYPE(p_bonded),POINTER::bonded
  TYPE(p_at),POINTER::at_aux
  INTEGER,DIMENSION(:,:),ALLOCATABLE::valores
 
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
              IF (priv==bonded%bonded_nods(kk)) THEN
                 valores(jj,gg+ll)=valores(jj,gg+ll)+1
              END IF
           END DO
        END DO
     END DO
     gg=gg+2
  END DO

  CALL SORT_INT_MATRIX (numats,gg,shell%ats_ord,valores)
  
  DEALLOCATE(valores)

  IF (symm%num_crits>0) THEN

     ALLOCATE(valores(numats,2))
     DO jj=1,numats
        at_aux=>shell%ats(jj)
        DO ll=1,2
           valores(jj,ll)=at_aux%bonded(ll)%num
        END DO
     END DO
 
     CALL SORT_INT_MATRIX (numats,2,shell%ats_ord,valores)

     DEALLOCATE(valores)

  END IF
 
  DO ll=1,2

     IF (symm%num_crits>0) THEN

        gg=0
        DO ii=1,numats
           IF (gg<shell%ats(ii)%bonded(ll)%num) gg=shell%ats(ii)%bonded(ll)%num
        END DO
        ALLOCATE(valores(numats,gg))
        valores(:,:)=0
        DO ii=1,numats
           valores(ii,1:shell%ats(ii)%bonded(ll)%num)=shell%ats(ii)%bonded(ll)%lev_supsets(:)
        END DO
        CALL SORT_INT_MATRIX (numats,gg,shell%ats_ord,valores)

        IF (symm%num_crits>0) THEN
           valores(:,:)=0
           DO ii=1,numats
              valores(ii,1:shell%ats(ii)%bonded(ll)%num)=shell%ats(ii)%bonded(ll)%cant_loops_1order(:)
           END DO
           ii=symm%num_crits
           CALL SORT_INT_MATRIX (numats,gg,shell%ats_ord,valores)
        END IF
        DEALLOCATE(valores)

     END IF

  END DO
 
  NULLIFY(bonded,at_aux)
 
END SUBROUTINE order_ats
 
SUBROUTINE order_ats_w_23order(shell)
 
  IMPLICIT NONE
 
  TYPE(p_shell),TARGET,INTENT(INOUT)::shell
 
  INTEGER::ii,gg,ll
  INTEGER::numats
  INTEGER,DIMENSION(:,:),ALLOCATABLE::valores
 
  numats=shell%nats
 
  DO ll=1,2

     IF (symm%num_crits>0) THEN
        
        gg=0
        DO ii=1,numats
           IF (gg<shell%ats(ii)%bonded(ll)%num) gg=shell%ats(ii)%bonded(ll)%num
        END DO
        ALLOCATE(valores(numats,gg))
        valores(:,:)=0
        
        DO ii=1,numats
           valores(ii,1:shell%ats(ii)%bonded(ll)%num)=shell%ats(ii)%bonded(ll)%cant_loops_2order(:)
        END DO
        
        CALL SORT_INT_MATRIX (numats,gg,shell%ats_ord,valores)
        
        IF (symm%num_crits>0) THEN
           DO ii=1,numats
              valores(ii,1:shell%ats(ii)%bonded(ll)%num)=shell%ats(ii)%bonded(ll)%cant_loops_3order(:)
           END DO
           CALL SORT_INT_MATRIX (numats,gg,shell%ats_ord,valores)
        END IF
     
        DEALLOCATE(valores)

     END IF
 
  END DO
 
END SUBROUTINE order_ats_w_23order

SUBROUTINE order_ats_w_2order(shell)
 
  IMPLICIT NONE
 
  TYPE(p_shell),TARGET,INTENT(INOUT)::shell
 
  INTEGER::ii,gg,ll
  INTEGER::numats
  INTEGER,DIMENSION(:,:),ALLOCATABLE::valores
 
  numats=shell%nats
 
  DO ll=1,2
     IF (symm%num_crits>0) THEN
        gg=0
        DO ii=1,numats
           IF (gg<shell%ats(ii)%bonded(ll)%num) gg=shell%ats(ii)%bonded(ll)%num
        END DO
        ALLOCATE(valores(numats,gg))
        valores(:,:)=0
        
        DO ii=1,numats
           valores(ii,1:shell%ats(ii)%bonded(ll)%num)=shell%ats(ii)%bonded(ll)%cant_loops_2order(:)
        END DO
        
        CALL SORT_INT_MATRIX (numats,gg,shell%ats_ord,valores)
        
        DEALLOCATE(valores)

     END IF
  END DO
 
END SUBROUTINE order_ats_w_2order



SUBROUTINE order_ats_w_next_shells(shell1st,shell2nd,nnods)
 
  IMPLICIT NONE
 
  INTEGER,INTENT(IN)::nnods
  TYPE(p_shell),TARGET,INTENT(INOUT)::shell1st
  TYPE(p_shell),DIMENSION(nnods),TARGET,INTENT(IN)::shell2nd
 
  ! tienen igual: privilegios, num hbs y num bs, y los mismos bonded en hbs y bs en su shell
  ! falta: comparar lo que sucede en las siguientes shell
 
  INTEGER::numats
  INTEGER::ii,jj,kk,gg,ll,mm,nn,qq,aa,aaa,dim_matrix
  INTEGER,DIMENSION(:),ALLOCATABLE::aux_ord
  INTEGER,DIMENSION(:,:),ALLOCATABLE::valores
  TYPE(p_at),POINTER::at_aux
  TYPE(p_shell),POINTER::shell_aux
 
  dim_matrix=0
  numats=shell1st%nats
   
  ALLOCATE(aux_ord(numats))
  aux_ord=shell1st%ats_ord(:)
 
  gg=0
  DO ii=1,symm%num_crits
     gg=gg+1
     DO jj=1,symm%crit(gg)
        gg=gg+1
        kk=aux_ord(symm%crit(gg))
        at_aux=>shell1st%ats(kk)
        aa=0
        DO ll=1,at_aux%hbs%num
           mm=at_aux%hbs%aux2sh(ll)
           aa=aa+shell2nd(mm)%nats2
        END DO
        DO ll=1,at_aux%bs%num
           mm=at_aux%bs%aux2sh(ll)
           aa=aa+shell2nd(mm)%nats2
        END DO
        IF (dim_matrix<aa) dim_matrix=aa
     END DO
  END DO
 
  ALLOCATE(valores(numats,dim_matrix))
  valores(:,:)=0
 
  gg=0
  DO ii=1,symm%num_crits
     gg=gg+1
     DO jj=1,symm%crit(gg)
        gg=gg+1
        kk=aux_ord(symm%crit(gg))
        at_aux=>shell1st%ats(kk)
        aa=0
        DO ll=1,at_aux%hbs%num
           mm=at_aux%hbs%aux2sh(ll)
           shell_aux=>shell2nd(mm)
           DO nn=1,shell_aux%nats
              qq=shell_aux%ats_ord(nn)
              aa=aa+1
              valores(kk,aa)=shell_aux%ats(qq)%hbs%num
              aa=aa+1
              valores(kk,aa)=shell_aux%ats(qq)%bs%num
           END DO
        END DO
        DO ll=1,at_aux%bs%num
           mm=at_aux%bs%aux2sh(ll)
           shell_aux=>shell2nd(mm)
           DO nn=1,shell_aux%nats
              qq=shell_aux%ats_ord(nn)
              aa=aa+1
              valores(kk,aa)=shell_aux%ats(qq)%hbs%num
              aa=aa+1
              valores(kk,aa)=shell_aux%ats(qq)%bs%num
           END DO
        END DO
     END DO
  END DO
 
  CALL SORT_INT_MATRIX (numats,dim_matrix,shell1st%ats_ord,valores)
  DEALLOCATE(valores,aux_ord)
 
 
  IF (symm%num_crits>0) THEN
 
     dim_matrix=0
     ALLOCATE(aux_ord(numats))
     aux_ord=shell1st%ats_ord(:)
 
     gg=0
     DO ii=1,symm%num_crits
        gg=gg+1
        DO jj=1,symm%crit(gg)
           gg=gg+1
           kk=aux_ord(symm%crit(gg))
           at_aux=>shell1st%ats(kk)
           aa=0
           DO ll=1,at_aux%hbs%num
              mm=at_aux%hbs%aux2sh(ll)
              aa=aa+shell2nd(mm)%nnods
           END DO
           DO ll=1,at_aux%bs%num
              mm=at_aux%bs%aux2sh(ll)
              aa=aa+shell2nd(mm)%nnods
           END DO
           IF (dim_matrix<aa) dim_matrix=aa
        END DO
     END DO
 
     IF (dim_matrix>0) THEN
 
        ALLOCATE(valores(numats,dim_matrix))
        valores(:,:)=0
 
        gg=0
        DO ii=1,symm%num_crits
           gg=gg+1
           DO jj=1,symm%crit(gg)
              gg=gg+1
              kk=aux_ord(symm%crit(gg))
              at_aux=>shell1st%ats(kk)
              aa=0
              DO ll=1,at_aux%hbs%num
                 mm=at_aux%hbs%aux2sh(ll)
                 shell_aux=>shell2nd(mm)
                 DO nn=1,shell_aux%nats
                    qq=shell_aux%ats_ord(nn)
                    aaa=aa+1
                    aa=aa+shell_aux%ats(qq)%hbs%num
                    valores(kk,aaa:aa)=shell_aux%ats(qq)%hbs%lev_supsets(:)
                    aaa=aa+1
                    aa=aa+shell_aux%ats(qq)%bs%num
                    valores(kk,aaa:aa)=shell_aux%ats(qq)%bs%lev_supsets(:)
                 END DO
              END DO
              DO ll=1,at_aux%bs%num
                 mm=at_aux%bs%aux2sh(ll)
                 shell_aux=>shell2nd(mm)
                 DO nn=1,shell_aux%nats
                    qq=shell_aux%ats_ord(nn)
                    aaa=aa+1
                    aa=aa+shell_aux%ats(qq)%hbs%num
                    valores(kk,aaa:aa)=shell_aux%ats(qq)%hbs%lev_supsets(:)
                    aaa=aa+1
                    aa=aa+shell_aux%ats(qq)%bs%num
                    valores(kk,aaa:aa)=shell_aux%ats(qq)%bs%lev_supsets(:)
                 END DO
              END DO
           END DO
        END DO
 
        CALL SORT_INT_MATRIX (numats,dim_matrix,shell1st%ats_ord,valores)
 
        DEALLOCATE(valores)
 
     END IF
 
     DEALLOCATE(aux_ord)
 
  END IF
 
  IF (symm%num_crits>0) THEN
 
     dim_matrix=0
     ALLOCATE(aux_ord(numats))
     aux_ord=shell1st%ats_ord(:)
 
     gg=0
     DO ii=1,symm%num_crits
        gg=gg+1
        DO jj=1,symm%crit(gg)
           gg=gg+1
           kk=aux_ord(symm%crit(gg))
           at_aux=>shell1st%ats(kk)
           aa=0
           DO ll=1,at_aux%hbs%num
              mm=at_aux%hbs%aux2sh(ll)
              aa=aa+shell2nd(mm)%nnods
           END DO
           DO ll=1,at_aux%bs%num
              mm=at_aux%bs%aux2sh(ll)
              aa=aa+shell2nd(mm)%nnods
           END DO
           IF (dim_matrix<aa) dim_matrix=aa
        END DO
     END DO
 
     IF (dim_matrix>0) THEN
 
        ALLOCATE(valores(numats,dim_matrix))
        valores(:,:)=0
 
        gg=0
        DO ii=1,symm%num_crits
           gg=gg+1
           DO jj=1,symm%crit(gg)
              gg=gg+1
              kk=aux_ord(symm%crit(gg))
              at_aux=>shell1st%ats(kk)
              aa=0
              DO ll=1,at_aux%hbs%num
                 mm=at_aux%hbs%aux2sh(ll)
                 shell_aux=>shell2nd(mm)
                 DO nn=1,shell_aux%nats
                    qq=shell_aux%ats_ord(nn)
                    aaa=aa+1
                    aa=aa+shell_aux%ats(qq)%hbs%num
                    valores(kk,aaa:aa)=shell_aux%ats(qq)%hbs%cant_loops_1order(:)
                    aaa=aa+1
                    aa=aa+shell_aux%ats(qq)%bs%num
                    valores(kk,aaa:aa)=shell_aux%ats(qq)%bs%cant_loops_1order(:)
                 END DO
              END DO
              DO ll=1,at_aux%bs%num
                 mm=at_aux%bs%aux2sh(ll)
                 shell_aux=>shell2nd(mm)
                 DO nn=1,shell_aux%nats
                    qq=shell_aux%ats_ord(nn)
                    aaa=aa+1
                    aa=aa+shell_aux%ats(qq)%hbs%num
                    valores(kk,aaa:aa)=shell_aux%ats(qq)%hbs%cant_loops_1order(:)
                    aaa=aa+1
                    aa=aa+shell_aux%ats(qq)%bs%num
                    valores(kk,aaa:aa)=shell_aux%ats(qq)%bs%cant_loops_1order(:)
                 END DO
              END DO
           END DO
        END DO
 
        CALL SORT_INT_MATRIX (numats,dim_matrix,shell1st%ats_ord,valores)
 
        DEALLOCATE(valores)
 
      END IF
 
     DEALLOCATE(aux_ord)
 
  END IF
 
  IF (symm%num_crits>0) THEN
 
     dim_matrix=0
     ALLOCATE(aux_ord(numats))
     aux_ord=shell1st%ats_ord(:)
 
     gg=0
     DO ii=1,symm%num_crits
        gg=gg+1
        DO jj=1,symm%crit(gg)
           gg=gg+1
           kk=aux_ord(symm%crit(gg))
           at_aux=>shell1st%ats(kk)
           aa=0
           DO ll=1,at_aux%hbs%num
              mm=at_aux%hbs%aux2sh(ll)
              aa=aa+shell2nd(mm)%nnods
           END DO
           DO ll=1,at_aux%bs%num
              mm=at_aux%bs%aux2sh(ll)
              aa=aa+shell2nd(mm)%nnods
           END DO
           IF (dim_matrix<aa) dim_matrix=aa
        END DO
     END DO
 
     IF (dim_matrix>0) THEN
 
        ALLOCATE(valores(numats,dim_matrix))
        valores(:,:)=0
 
        gg=0
        DO ii=1,symm%num_crits
           gg=gg+1
           DO jj=1,symm%crit(gg)
              gg=gg+1
              kk=aux_ord(symm%crit(gg))
              at_aux=>shell1st%ats(kk)
              aa=0
              DO ll=1,at_aux%hbs%num
                 mm=at_aux%hbs%aux2sh(ll)
                 shell_aux=>shell2nd(mm)
                 DO nn=1,shell_aux%nats
                    qq=shell_aux%ats_ord(nn)
                    aaa=aa+1
                    aa=aa+shell_aux%ats(qq)%hbs%num
                    valores(kk,aaa:aa)=shell_aux%ats(qq)%hbs%cant_loops_2order(:)
                    aaa=aa+1
                    aa=aa+shell_aux%ats(qq)%bs%num
                    valores(kk,aaa:aa)=shell_aux%ats(qq)%bs%cant_loops_2order(:)
                 END DO
              END DO
              DO ll=1,at_aux%bs%num
                 mm=at_aux%bs%aux2sh(ll)
                 shell_aux=>shell2nd(mm)
                 DO nn=1,shell_aux%nats
                    qq=shell_aux%ats_ord(nn)
                    aaa=aa+1
                    aa=aa+shell_aux%ats(qq)%hbs%num
                    valores(kk,aaa:aa)=shell_aux%ats(qq)%hbs%cant_loops_2order(:)
                    aaa=aa+1
                    aa=aa+shell_aux%ats(qq)%bs%num
                    valores(kk,aaa:aa)=shell_aux%ats(qq)%bs%cant_loops_2order(:)
                 END DO
              END DO
           END DO
        END DO
 
        CALL SORT_INT_MATRIX (numats,dim_matrix,shell1st%ats_ord,valores)
 
        DEALLOCATE(valores)
 
      END IF
 
     DEALLOCATE(aux_ord)
 
  END IF
 
  IF (symm%num_crits>0) THEN
 
     dim_matrix=0
     ALLOCATE(aux_ord(numats))
     aux_ord=shell1st%ats_ord(:)
 
     gg=0
     DO ii=1,symm%num_crits
        gg=gg+1
        DO jj=1,symm%crit(gg)
           gg=gg+1
           kk=aux_ord(symm%crit(gg))
           at_aux=>shell1st%ats(kk)
           aa=0
           DO ll=1,at_aux%hbs%num
              mm=at_aux%hbs%aux2sh(ll)
              aa=aa+shell2nd(mm)%nnods
           END DO
           DO ll=1,at_aux%bs%num
              mm=at_aux%bs%aux2sh(ll)
              aa=aa+shell2nd(mm)%nnods
           END DO
           IF (dim_matrix<aa) dim_matrix=aa
        END DO
     END DO
 
     IF (dim_matrix>0) THEN
 
        ALLOCATE(valores(numats,dim_matrix))
        valores(:,:)=0
 
        gg=0
        DO ii=1,symm%num_crits
           gg=gg+1
           DO jj=1,symm%crit(gg)
              gg=gg+1
              kk=aux_ord(symm%crit(gg))
              at_aux=>shell1st%ats(kk)
              aa=0
              DO ll=1,at_aux%hbs%num
                 mm=at_aux%hbs%aux2sh(ll)
                 shell_aux=>shell2nd(mm)
                 DO nn=1,shell_aux%nats
                    qq=shell_aux%ats_ord(nn)
                    aaa=aa+1
                    aa=aa+shell_aux%ats(qq)%hbs%num
                    valores(kk,aaa:aa)=shell_aux%ats(qq)%hbs%cant_loops_3order(:)
                    aaa=aa+1
                    aa=aa+shell_aux%ats(qq)%bs%num
                    valores(kk,aaa:aa)=shell_aux%ats(qq)%bs%cant_loops_3order(:)
                 END DO
              END DO
              DO ll=1,at_aux%bs%num
                 mm=at_aux%bs%aux2sh(ll)
                 shell_aux=>shell2nd(mm)
                 DO nn=1,shell_aux%nats
                    qq=shell_aux%ats_ord(nn)
                    aaa=aa+1
                    aa=aa+shell_aux%ats(qq)%hbs%num
                    valores(kk,aaa:aa)=shell_aux%ats(qq)%hbs%cant_loops_3order(:)
                    aaa=aa+1
                    aa=aa+shell_aux%ats(qq)%bs%num
                    valores(kk,aaa:aa)=shell_aux%ats(qq)%bs%cant_loops_3order(:)
                 END DO
              END DO
           END DO
        END DO
 
        CALL SORT_INT_MATRIX (numats,dim_matrix,shell1st%ats_ord,valores)
 
        DEALLOCATE(valores)
 
      END IF
 
     DEALLOCATE(aux_ord)
 
  END IF
  
END SUBROUTINE order_ats_w_next_shells
 
SUBROUTINE order_bonded(dim_privil,privilegios,bonded)
  
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
        IF (priv==bonded%bonded_nods(hh)) THEN
           valores(hh,ii)=1
        END IF
     END DO
  END DO
  valores(:,gg)=bonded%lev_supsets(:)
   
  CALL SORT_INT_MATRIX (numb,gg,bonded%order,valores)
  DEALLOCATE(valores)

  IF (symm%num_crits>0) THEN  !Miro los loops en 1er orden

     gg=1
     ALLOCATE(valores(numb,gg))
     valores(:,1)=bonded%cant_loops_1order(:)
     CALL SORT_INT_MATRIX (numb,gg,bonded%order,valores)
     DEALLOCATE(valores)
 
  END IF
 
  CALL REORDER_BONDED (bonded)
 
END SUBROUTINE order_bonded
 
SUBROUTINE order_bonded_w_23order(bonded)
  
  IMPLICIT NONE
 
  TYPE(p_bonded),INTENT(INOUT)::bonded
  
  INTEGER::numb
  INTEGER,DIMENSION(:,:),ALLOCATABLE::valores
 
  numb=bonded%num
   
   
  ALLOCATE(valores(numb,2))
  valores(:,1)=bonded%cant_loops_2order(:)
  valores(:,2)=bonded%cant_loops_3order(:)
   
  CALL SORT_INT_MATRIX (numb,2,bonded%order,valores)
 
  CALL REORDER_BONDED (bonded)
   
  DEALLOCATE(valores)
 
END SUBROUTINE order_bonded_w_23order
 
SUBROUTINE order_bonded_w_2order(bonded)
  
  IMPLICIT NONE
 
  TYPE(p_bonded),INTENT(INOUT)::bonded
  
  INTEGER::numb
  INTEGER,DIMENSION(:,:),ALLOCATABLE::valores
  INTEGER::ii
 
  numb=bonded%num
   
  ii=symm%num_crits
  ALLOCATE(valores(numb,1))
  valores(:,1)=bonded%cant_loops_2order(:)
   
  CALL SORT_INT_MATRIX (numb,1,bonded%order,valores)
 
  CALL REORDER_BONDED (bonded)

  DEALLOCATE(valores)
 
  IF (symm%num_crits<ii) THEN
     print*,'AQUIIIIIIIIIIIIIIIIIIIIIIII'
  END IF

END SUBROUTINE order_bonded_w_2order
 
  
 
SUBROUTINE order_bonded_w_next_shells(bonded,shell2nd,nnods)

  IMPLICIT NONE
 
  INTEGER,INTENT(IN)::nnods
  TYPE(p_bonded),INTENT(INOUT)::bonded
  TYPE(p_shell),DIMENSION(nnods),TARGET,INTENT(IN)::shell2nd
 
  INTEGER::numb,dim_matrix
  INTEGER::ii,jj,gg,mm,kk,aa,qq,aaa,nn
  INTEGER,DIMENSION(:,:),ALLOCATABLE::valores
  TYPE(p_shell),POINTER::shell_aux
 
  dim_matrix=0
  numb=bonded%num
 
  gg=0
  DO ii=1,symm%num_crits
     gg=gg+1
     DO jj=1,symm%crit(gg)
        gg=gg+1
        kk=bonded%order(symm%crit(gg))
        mm=bonded%aux2sh(kk)
        IF (dim_matrix<shell2nd(mm)%nats2) dim_matrix=shell2nd(mm)%nats2
     END DO
  END DO
 
  ALLOCATE(valores(numb,dim_matrix))
  valores(:,:)=0
 
  gg=0
  DO ii=1,symm%num_crits
     gg=gg+1
     DO jj=1,symm%crit(gg)
        gg=gg+1
        kk=bonded%order(symm%crit(gg))
        mm=bonded%aux2sh(kk)
        shell_aux=>shell2nd(mm)
        aa=0
        DO nn=1,shell_aux%nats
           qq=shell_aux%ats_ord(nn)
           DO ll=1,2
              aa=aa+1
              valores(kk,aa)=shell_aux%ats(qq)%bonded(ll)%num
           END DO
        END DO
     END DO
  END DO
 
  CALL SORT_INT_MATRIX (numb,dim_matrix,bonded%order,valores)
  DEALLOCATE(valores)

 
  IF (symm%num_crits>0) THEN
 
     dim_matrix=0
     numb=bonded%num
     
     gg=0
     DO ii=1,symm%num_crits
        gg=gg+1
        DO jj=1,symm%crit(gg)
           gg=gg+1
           kk=bonded%order(symm%crit(gg))
           mm=bonded%aux2sh(kk)
           IF (dim_matrix<shell2nd(mm)%nnods) dim_matrix=shell2nd(mm)%nnods
        END DO
     END DO
     
     IF (dim_matrix>0) THEN
 
        ALLOCATE(valores(numb,dim_matrix))
        valores(:,:)=0
 
        gg=0
        DO ii=1,symm%num_crits
           gg=gg+1
           DO jj=1,symm%crit(gg)
              gg=gg+1
              kk=bonded%order(symm%crit(gg))
              mm=bonded%aux2sh(kk)
              shell_aux=>shell2nd(mm)
              aa=0
              DO nn=1,shell_aux%nats
                 qq=shell_aux%ats_ord(nn)
                 DO ll=1,2
                    aaa=aa+1
                    aa=aa+shell_aux%ats(qq)%bonds(ll)%num
                    valores(kk,aaa:aa)=shell_aux%ats(qq)%bonds(ll)%lev_supsets(:)
                 END DO
              END DO
           END DO
        END DO
 
        CALL SORT_INT_MATRIX (numb,dim_matrix,bonded%order,valores)
        DEALLOCATE(valores)
 
     END IF
  END IF
 
 
  IF (symm%num_crits>0) THEN
 
     dim_matrix=0
     numb=bonded%num
     
     gg=0
     DO ii=1,symm%num_crits
        gg=gg+1
        DO jj=1,symm%crit(gg)
           gg=gg+1
           kk=bonded%order(symm%crit(gg))
           mm=bonded%aux2sh(kk)
           IF (dim_matrix<shell2nd(mm)%nnods) dim_matrix=shell2nd(mm)%nnods
        END DO
     END DO
     
     IF (dim_matrix>0) THEN
 
        ALLOCATE(valores(numb,dim_matrix))
        valores(:,:)=0
 
        gg=0
        DO ii=1,symm%num_crits
           gg=gg+1
           DO jj=1,symm%crit(gg)
              gg=gg+1
              kk=bonded%order(symm%crit(gg))
              mm=bonded%aux2sh(kk)
              shell_aux=>shell2nd(mm)
              aa=0
              DO nn=1,shell_aux%nats
                 qq=shell_aux%ats_ord(nn)
                 DO ll=1,2
                    aaa=aa+1
                    aa=aa+shell_aux%ats(qq)%bonds(ll)%num
                    valores(kk,aaa:aa)=shell_aux%ats(qq)%bonds(ll)%cant_loops_1order(:)
                 END DO
              END DO
           END DO
        END DO
 
        CALL SORT_INT_MATRIX (numb,dim_matrix,bonded%order,valores)
        DEALLOCATE(valores)

     END IF
  END IF
 
  IF (symm%num_crits>0) THEN
 
     dim_matrix=0
     numb=bonded%num
     
     gg=0
     DO ii=1,symm%num_crits
        gg=gg+1
        DO jj=1,symm%crit(gg)
           gg=gg+1
           kk=bonded%order(symm%crit(gg))
           mm=bonded%aux2sh(kk)
           IF (dim_matrix<shell2nd(mm)%nnods) dim_matrix=shell2nd(mm)%nnods
        END DO
     END DO
     
     IF (dim_matrix>0) THEN
 
        ALLOCATE(valores(numb,dim_matrix))
        valores(:,:)=0
 
        gg=0
        DO ii=1,symm%num_crits
           gg=gg+1
           DO jj=1,symm%crit(gg)
              gg=gg+1
              kk=bonded%order(symm%crit(gg))
              mm=bonded%aux2sh(kk)
              shell_aux=>shell2nd(mm)
              aa=0
              DO nn=1,shell_aux%nats
                 qq=shell_aux%ats_ord(nn)
                 DO ll=1,2
                    aaa=aa+1
                    aa=aa+shell_aux%ats(qq)%bonds(ll)%num
                    valores(kk,aaa:aa)=shell_aux%ats(qq)%bonds(ll)%cant_loops_2order(:)
                 END DO
              END DO
           END DO
        END DO
 
        CALL SORT_INT_MATRIX (numb,dim_matrix,bonded%order,valores)
        DEALLOCATE(valores)

     END IF
  END IF
 
  IF (symm%num_crits>0) THEN
 
     dim_matrix=0
     numb=bonded%num
     gg=0
     DO ii=1,symm%num_crits
        gg=gg+1
        DO jj=1,symm%crit(gg)
           gg=gg+1
           kk=bonded%order(symm%crit(gg))
           mm=bonded%aux2sh(kk)
           IF (dim_matrix<shell2nd(mm)%nnods) dim_matrix=shell2nd(mm)%nnods
        END DO
     END DO
     
     IF (dim_matrix>0) THEN
 
        ALLOCATE(valores(numb,dim_matrix))
        valores(:,:)=0
 
        gg=0
        DO ii=1,symm%num_crits
           gg=gg+1
           DO jj=1,symm%crit(gg)
              gg=gg+1
              kk=bonded%order(symm%crit(gg))
              mm=bonded%aux2sh(kk)
              shell_aux=>shell2nd(mm)
              aa=0
              DO nn=1,shell_aux%nats
                 qq=shell_aux%ats_ord(nn)
                 DO ll=1,2
                    aaa=aa+1
                    aa=aa+shell_aux%ats(qq)%bonds(ll)%num
                    valores(kk,aaa:aa)=shell_aux%ats(qq)%bonds(ll)%cant_loops_3order(:)
                 END DO
              END DO
           END DO
        END DO
 
        CALL SORT_INT_MATRIX (numb,dim_matrix,bonded%order,valores)
        DEALLOCATE(valores)

     END IF
  END IF
 
  CALL REORDER_BONDED (bonded)
 
END SUBROUTINE order_bonded_w_next_shells
 
SUBROUTINE REORDER_BONDED (bonded)
 
  IMPLICIT NONE
 
  TYPE(p_bonded),INTENT(INOUT)::bonded
  INTEGER,DIMENSION(:),ALLOCATABLE::box
  INTEGER::ii,numb
 
  ALLOCATE(box(bonded%num))
  box(:)=bonded%bonded_ats(:)
  bonded%bonded_ats(:)=box(bonded%order(:))
  box(:)=bonded%bonded_nods(:)
  bonded%bonded_nods(:)=box(bonded%order(:))
  box(:)=bonded%lev_supsets(:)
  bonded%lev_supsets(:)=box(bonded%order(:))
  box(:)=bonded%aux2sh(:)
  bonded%aux2sh(:)=box(bonded%order(:))
  box(:)=bonded%cant_loops_1order(:)
  bonded%cant_loops_1order(:)=box(bonded%order(:))
  box(:)=bonded%cant_loops_2order(:)
  bonded%cant_loops_2order(:)=box(bonded%order(:))
  box(:)=bonded%cant_loops_3order(:)
  bonded%cant_loops_3order(:)=box(bonded%order(:))
  DEALLOCATE(box)
 
  numb=bonded%num
  bonded%order(:)=(/(ii,ii=1,numb)/)
 
END SUBROUTINE REORDER_BONDED
 
SUBROUTINE COMPLETE_WSYMM(shell)
 
  TYPE(p_shell),TARGET,INTENT(INOUT)::shell
  TYPE(p_at),POINTER::at_aux
  TYPE(p_bonded),POINTER::bs_aux
 
  INTEGER::mm,nn,ll,ii
 
  DO ii=1,shell%nats
     at_aux=>shell%ats(ii)
     DO ll=1,2
        IF (at_aux%bonds(ll)%symm%num_crits>0) THEN
           bs_aux=>at_aux%bonds(ll)
           mm=0
           DO nn=1,bs_aux%symm%num_crits
              mm=mm+1
              DO ll=1,bs_aux%symm%crit(mm)
                 mm=mm+1
                 bs_aux%wsymm(bs_aux%symm%crit(mm))=1
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
  
  NULLIFY(at_aux,bs_aux)
 
END SUBROUTINE COMPLETE_WSYMM
 
SUBROUTINE build_mss(shell)
 
  IMPLICIT NONE
 
  TYPE(p_shell),INTENT(INOUT)::shell
 
  INTEGER::ntot,nats,nats2,nnods
  INTEGER::ii,jj,kk
  INTEGER,DIMENSION(:),ALLOCATABLE::oo
  
 
  IF (ALLOCATED(mss_ind_ats))   DEALLOCATE(mss_ind_ats)
  IF (ALLOCATED(mss_ind_nods))  DEALLOCATE(mss_ind_nods)
  IF (ALLOCATED(mss_symm))      DEALLOCATE(mss_symm)
 
  ntot = shell%ntot
  nats = shell%nats
  nats2 = shell%nats2
  nnods = shell%nnods

  ALLOCATE(oo(nats))
  ALLOCATE(mss_ind_ats(ntot),mss_ind_nods(ntot),mss_symm(ntot))
 
  CALL COMPLETE_WSYMM(shell)
 
  oo(:)=shell%ats_ord(:)
  mss_symm(:)=0
 
  jj=1
  mss_ind_ats(1)  = nats
  mss_ind_nods(1) = nats
  ii=jj+1
  jj=jj+nats
  mss_symm(ii:jj)     = (/(shell%wsymm(kk),kk=1,nats)/)
  mss_ind_ats(ii:jj)  = (/(shell%ats(oo(kk))%ind,kk=1,nats)/) !! ojo
  mss_ind_nods(ii:jj) = at2nod(mss_ind_ats(ii:jj))
  ii=jj+1
  jj=jj+nats2
  mss_ind_ats(ii:jj)  = (/((/shell%ats(oo(kk))%bonded(1)%num,shell%ats(oo(kk))%bonded(2)%num/),kk=1,nats)/)
  mss_ind_nods(ii:jj) = mss_ind_ats(ii:jj)
  ii=jj+1
  jj=jj+nnods
  mss_symm(ii:jj)     = (/((/shell%ats(oo(kk))%bonded(1)%wsymm(:),shell%ats(oo(kk))%bonded(2)%wsymm(:)/),kk=1,nats)/)
  mss_ind_ats(ii:jj)  = (/((/shell%ats(oo(kk))%bonded(1)%ats(:),shell%ats(oo(kk))%bonded(2)%ats(:)/),kk=1,nats)/)
  mss_ind_nods(ii:jj) = at2nod(mss_ind_ats(ii:jj))
 
  DEALLOCATE(oo)
 
END SUBROUTINE build_mss
 
 
SUBROUTINE build_mss_wout(shell2nd,aux_ind_ats,aux_ind_nods,aux_symm)
 
  IMPLICIT NONE
 
  TYPE(p_shell),INTENT(INOUT)::shell2nd
  INTEGER,DIMENSION(shell2nd%ntot),INTENT(OUT)::aux_ind_ats,aux_ind_nods,aux_symm
 
  INTEGER::ntot,nats,nats2,nnods
  INTEGER::ii,jj,kk
  INTEGER,DIMENSION(:),ALLOCATABLE::oo
 
  ntot = shell2nd%ntot
  nats = shell2nd%nats
  nats2 = shell2nd%nats2
  nnods = shell2nd%nnods
 
  ALLOCATE(oo(nats))
 
  oo(:)=shell2nd%ats_ord(:)
 
  CALL COMPLETE_WSYMM(shell2nd)
 
  aux_symm(:)=0
 
  jj=1
  aux_ind_ats(1)  = nats
  aux_ind_nods(1) = nats
  ii=jj+1
  jj=jj+nats
  aux_symm(ii:jj)     = (/(shell2nd%wsymm(kk),kk=1,nats)/)
  aux_ind_ats(ii:jj)  = (/(shell2nd%ats(oo(kk))%ind,kk=1,nats)/) !! ojo
  aux_ind_nods(ii:jj) = at2nod(aux_ind_ats(ii:jj))
  ii=jj+1
  jj=jj+nats2
  aux_ind_ats(ii:jj)  = (/((/shell2nd%ats(oo(kk))%bonded(1)%num,shell2nd%ats(oo(kk))%bonded(2)%num/),kk=1,nats)/)
  aux_ind_nods(ii:jj) = aux_ind_ats(ii:jj)
  ii=jj+1
  jj=jj+nnods
  aux_symm(ii:jj)     = (/((/shell2nd%ats(oo(kk))%bonded(1)%wsymm(:),shell2nd%ats(oo(kk))%bonded(2)%wsymm(:)/),kk=1,nats)/)
  aux_ind_ats(ii:jj)  = (/((/shell2nd%ats(oo(kk))%bonded(1)%ats(:),shell2nd%ats(oo(kk))%bonded(2)%ats(:)/),kk=1,nats)/)
  aux_ind_nods(ii:jj) = at2nod(aux_ind_ats(ii:jj))
 
  DEALLOCATE(oo)
 
END SUBROUTINE build_mss_wout
 
SUBROUTINE build_mss_wout_word(shell2nd,aux_ind_ats,aux_ind_nods,aux_symm,aux_ord)
  
  IMPLICIT NONE
 
  TYPE(p_shell),INTENT(INOUT)::shell2nd
  INTEGER,DIMENSION(shell2nd%ntot),INTENT(OUT)::aux_ind_ats,aux_ind_nods,aux_symm
  INTEGER,DIMENSION(shell2nd%nnods),INTENT(OUT)::aux_ord
 
  INTEGER::ntot,nats,nats2,nnods
  INTEGER::ii,jj,kk
  INTEGER,DIMENSION(:),ALLOCATABLE::oo
 
  ntot = shell2nd%ntot
  nats = shell2nd%nats
  nats2 = shell2nd%nats2
  nnods = shell2nd%nnods
 
  ALLOCATE(oo(nats))
 
  oo(:)=shell2nd%ats_ord(:)
 
  CALL COMPLETE_WSYMM(shell2nd)
  aux_symm(:)=0
 
  jj=1
  aux_ind_ats(1)  = nats
  aux_ind_nods(1) = nats
  ii=jj+1
  jj=jj+nats
  aux_symm(ii:jj)     = (/(shell2nd%wsymm(kk),kk=1,nats)/)
  aux_ind_ats(ii:jj)  = (/(shell2nd%ats(oo(kk))%ind,kk=1,nats)/) !! ojo
  aux_ind_nods(ii:jj) = at2nod(aux_ind_ats(ii:jj))
  ii=jj+1
  jj=jj+nats2
  aux_ind_ats(ii:jj)  = (/((/shell2nd%ats(oo(kk))%bonded(1)%num,shell2nd%ats(oo(kk))%bonded(2)%num/),kk=1,nats)/)
  aux_ind_nods(ii:jj) = aux_ind_ats(ii:jj)
  ii=jj+1
  jj=jj+nnods
  aux_symm(ii:jj)     = (/((/shell2nd%ats(oo(kk))%bonded(1)%wsymm(:),shell2nd%ats(oo(kk))%bonded(2)%wsymm(:)/),kk=1,nats)/)
  aux_ord(:)          = (/((/shell2nd%ats(oo(kk))%bonded(1)%ats(:),shell2nd%ats(oo(kk))%bonded(2)%ats(:)/),kk=1,nats)/)
  aux_ind_ats(ii:jj)  = aux_ord(:)
  aux_ord(:)          = at2nod(aux_ord(:))
  aux_ind_nods(ii:jj) = aux_ord(:)
  aux_ord(:)          = (/((/shell2nd%ats(oo(kk))%bonded(1)%aux2sh(:),shell2nd%ats(oo(kk))%bonded(2)%aux2sh(:)/),kk=1,nats)/)
 
  DEALLOCATE(oo)
 
END SUBROUTINE build_mss_wout_word
 
 
SUBROUTINE remato_mss()
 
  IMPLICIT NONE
 
  INTEGER::ii,in_node,out_node,flaga,flagz,aa,bb,tnod,tat
  INTEGER,DIMENSION(num_nods)::translation
  LOGICAL,DIMENSION(num_nods)::filtro
 
  filtro=.FALSE.
 
  IF (ALLOCATED(mss))  DEALLOCATE(mss)
  ALLOCATE(mss(SIZE(mss_ind_nods)))
 
  CALL RESET_TRANSLATOR()
 
  flaga=1
  flagz=2
  aa=mss_ind_nods(flaga)
  in_node=mss_ind_nods(flagz)
  tnod=trad2py_nod(in_node)
  tat=trad2py_at(mss_ind_ats(flagz))
  IF (filtro(in_node).eqv..FALSE.) THEN
     filtro(in_node)=.TRUE.
     CALL TRANSLATION_WO_SUPSETS(in_node,out_node)
     translation(in_node)=out_node
  ELSE
     out_node=translation(in_node)
  END IF
  mss(flaga)=aa
  flaga=flagz
  flagz=flagz+aa-1
  DO ii=flaga,flagz
     mss(ii)=out_node
     mss_ind_nods(ii)=tnod
     mss_ind_ats(ii)=tat
  END DO
  flaga=flagz+1
  flagz=flagz+2*aa
  bb=0
  DO ii=flaga,flagz
     bb=bb+mss_ind_nods(ii)
     mss(ii)=mss_ind_nods(ii)
  END DO
  flaga=flagz+1
  flagz=flagz+bb
  DO ii=flaga,flagz
     in_node=mss_ind_nods(ii)
     IF (filtro(in_node).eqv..FALSE.) THEN
        filtro(in_node)=.TRUE.
        CALL TRANSLATION_WO_SUPSETS(in_node,out_node)
        translation(in_node)=out_node
     ELSE
        out_node=translation(in_node)
     END IF
     mss(ii)=out_node
     mss_ind_nods(ii)=trad2py_nod(mss_ind_nods(ii))
     mss_ind_ats(ii)=trad2py_at(mss_ind_ats(ii))
  END DO
 
END SUBROUTINE remato_mss
 
SUBROUTINE remato_mss2sh()
 
  IMPLICIT NONE
 
  INTEGER::ii,jj,in_node,out_node,flaga,flagz,aa,bb,cc,tnod,tat
  INTEGER,DIMENSION(num_nods)::translation
  LOGICAL,DIMENSION(num_nods)::filtro
 
  filtro=.FALSE.
 
  IF (ALLOCATED(mss))  DEALLOCATE(mss)
  ALLOCATE(mss(SIZE(mss_ind_nods)))
 
  CALL RESET_TRANSLATOR()

  flaga=1
  flagz=2
  aa=mss_ind_nods(flaga)
  in_node=mss_ind_nods(flagz)
  tnod=trad2py_nod(in_node)
  tat=trad2py_at(mss_ind_ats(flagz))
  IF (filtro(in_node).eqv..FALSE.) THEN
     filtro(in_node)=.TRUE.
     CALL TRANSLATION_WO_SUPSETS(in_node,out_node)
     translation(in_node)=out_node
  ELSE
     out_node=translation(in_node)
  END IF
  mss(flaga)=aa
  flaga=flagz
  flagz=flagz+aa-1
  DO ii=flaga,flagz
     mss(ii)=out_node
     mss_ind_nods(ii)=tnod
     mss_ind_ats(ii)=tat
  END DO
  flaga=flagz+1
  flagz=flagz+2*aa
  bb=0
  DO ii=flaga,flagz
     bb=bb+mss_ind_nods(ii)
     mss(ii)=mss_ind_nods(ii)
  END DO
  flaga=flagz+1
  flagz=flagz+bb
  DO ii=flaga,flagz
     in_node=mss_ind_nods(ii)
     IF (filtro(in_node).eqv..FALSE.) THEN
        filtro(in_node)=.TRUE.
        CALL TRANSLATION_WO_SUPSETS(in_node,out_node)
        translation(in_node)=out_node
     ELSE
        out_node=translation(in_node)
     END IF
     mss(ii)=out_node
     mss_ind_nods(ii)=trad2py_nod(mss_ind_nods(ii))
     mss_ind_ats(ii)=trad2py_at(mss_ind_ats(ii))
  END DO
  DO jj=1,bb
     flaga=flagz+1
     flagz=flagz+2
     aa=mss_ind_nods(flaga)
     in_node=mss_ind_nods(flagz)
     tnod=trad2py_nod(in_node)
     tat=trad2py_at(mss_ind_ats(flagz))
     IF (filtro(in_node).eqv..FALSE.) THEN
        filtro(in_node)=.TRUE.
        CALL TRANSLATION_WO_SUPSETS(in_node,out_node)
        translation(in_node)=out_node
     ELSE
        out_node=translation(in_node)
     END IF
     mss(flaga)=aa
     flaga=flagz
     flagz=flagz+aa-1
     DO ii=flaga,flagz
        mss(ii)=out_node
        mss_ind_nods(ii)=tnod
        mss_ind_ats(ii)=tat
     END DO
     flaga=flagz+1
     flagz=flagz+2*aa
     cc=0
     DO ii=flaga,flagz
        cc=cc+mss_ind_nods(ii)
        mss(ii)=mss_ind_nods(ii)
     END DO
     flaga=flagz+1
     flagz=flagz+cc
     DO ii=flaga,flagz
        in_node=mss_ind_nods(ii)
        IF (filtro(in_node).eqv..FALSE.) THEN
           filtro(in_node)=.TRUE.
           CALL TRANSLATION_WO_SUPSETS(in_node,out_node)
           translation(in_node)=out_node
        ELSE
           out_node=translation(in_node)
        END IF
        mss(ii)=out_node
        mss_ind_nods(ii)=trad2py_nod(mss_ind_nods(ii))
        mss_ind_ats(ii)=trad2py_at(mss_ind_ats(ii))
     END DO
  END DO

 
END SUBROUTINE remato_mss2sh

 
SUBROUTINE TRANSLATION_WO_SUPSETS(in_node,out_node)
 
  IMPLICIT NONE
 
  INTEGER,INTENT(IN)::in_node
  INTEGER,INTENT(OUT)::out_node
 
  INTEGER::ii_supset,aa
  TYPE(p_trans_set),POINTER::supset_aux
 
  ii_supset=lev_supsets(in_node)
  supset_aux=>translator%supset(ii_supset)
 
  aa=supset_aux%count_sets+1
  supset_aux%count_sets=aa
  out_node=supset_aux%tope+aa
 
  NULLIFY(supset_aux)
 
END SUBROUTINE TRANSLATION_WO_SUPSETS
 
 
 
 
 
!!!!!!!!!!!!!!!!!!!!!##########
!!!!#################################
!!!!####  GENERAL FUNCTIONS TO SORT  
!!!!#################################
 
 
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

SUBROUTINE BUILDING_REPE_LIST (carro,dim_carro,list_norepe,dim_norepe,list_repe,cant_repe,dim_repe)
 
  IMPLICIT NONE
  INTEGER,INTENT(IN)::dim_carro
  INTEGER,DIMENSION(dim_carro),INTENT(IN)::carro
  INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(INOUT)::list_norepe,list_repe,cant_repe
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
 
  ALLOCATE(list_norepe(dim_norepe))
  IF (dim_norepe<dim_carro) THEN
     jj=0
     kk=0
     DO ii=1,dim_carro
        IF (filtro2(ii).eqv..TRUE.) THEN
           jj=jj+1
           list_norepe(jj)=carro(ii)
           IF (repetitions(ii)>1) THEN
              kk=kk+1
              box(kk)=carro(ii)
              box2(kk)=repetitions(ii)
           END IF
        END IF
     END DO
     ALLOCATE(list_repe(kk),cant_repe(kk))
     dim_repe=kk
     list_repe(:)=box(1:kk)
     cant_repe(:)=box2(1:kk)
  ELSE
     dim_repe=0
     list_norepe(:)=carro(:)
     ALLOCATE(list_repe(0),cant_repe(0))
  END IF
 
END SUBROUTINE BUILDING_REPE_LIST
 
SUBROUTINE BUILDING_REPE_LIST2 (core,carro,dim_carro,list_repe,cant_repe,dim_repe)
 
  IMPLICIT NONE
  INTEGER,INTENT(IN)::dim_carro,core
  INTEGER,DIMENSION(dim_carro),INTENT(IN)::carro
  INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(INOUT)::list_repe,cant_repe
  INTEGER,INTENT(OUT)::dim_repe
 
  LOGICAL,DIMENSION(dim_carro)::filtro,filtro2
  INTEGER,DIMENSION(dim_carro)::box,repetitions,box2
  INTEGER::ii,jj,gg,kk,dim_norepe
 
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
 
  IF (dim_norepe<dim_carro) THEN
     jj=0
     kk=0
     DO ii=1,dim_carro
        IF (filtro2(ii).eqv..TRUE.) THEN
           jj=jj+1
           IF (repetitions(ii)>1) THEN
              IF (carro(ii)/=core) THEN
                 kk=kk+1
                 box(kk)=carro(ii)
                 box2(kk)=repetitions(ii)
              END IF
           END IF
        END IF
     END DO
     ALLOCATE(list_repe(kk),cant_repe(kk))
     dim_repe=kk
     list_repe(:)=box(1:kk)
     cant_repe(:)=box2(1:kk)
  ELSE
     dim_repe=0
     ALLOCATE(list_repe(0),cant_repe(0))
  END IF
 
 
END SUBROUTINE BUILDING_REPE_LIST2

SUBROUTINE COUNT_NODE_REPE (list_repe,cant_repe,dim_repe,nodes,nnods,salida)
 
  IMPLICIT NONE
 
  INTEGER,INTENT(IN)::dim_repe,nnods
  INTEGER,DIMENSION(dim_repe),INTENT(IN)::list_repe,cant_repe
  INTEGER,DIMENSION(nnods),INTENT(IN)::nodes
  INTEGER,DIMENSION(nnods),INTENT(OUT)::salida
 
  INTEGER::ii,jj,kk
 
  DO ii=1,dim_repe
     jj=list_repe(ii)
     DO kk=1,nnods
        IF (jj==nodes(kk)) THEN
           salida(kk)=cant_repe(ii)
        END IF
     END DO
  END DO
 
END SUBROUTINE COUNT_NODE_REPE
 

SUBROUTINE quito_falsos_bonded_symm(bonded)

  IMPLICIT NONE
  
  TYPE(p_bonded),INTENT(INOUT)::bonded

  INTEGER::ii,jj,gg,kk,ll
  INTEGER::new_num_crits
  INTEGER,DIMENSION(symm%length)::new_crit
  LOGICAL::interruptor,interruptor2

  new_num_crits=0
  ll=0

  interruptor2=.FALSE.
  gg=0
  DO ii=1,symm%num_crits
     gg=gg+1
     kk=symm%crit(gg+1)
     interruptor=.FALSE.
     IF (bonded%cant_loops_1order(kk)==0) THEN
        IF (bonded%cant_loops_2order(kk)==0) THEN
           IF (bonded%cant_loops_3order(kk)==0) THEN
              interruptor=.TRUE.
           END IF
        END IF
     END IF
     IF (interruptor.eqv..TRUE.) THEN
        interruptor2=.TRUE.
        gg=gg+symm%crit(gg)
     ELSE
        new_num_crits=new_num_crits+1
        ll=ll+1
        new_crit(ll)=symm%crit(gg)
        DO jj=1,symm%crit(gg)
           ll=ll+1
           gg=gg+1
           new_crit(ll)=symm%crit(gg)
        END DO
     END IF
  END DO


  IF (new_num_crits==0) THEN
     symm%num_crits=0
  ELSE
     bonded%unsolved_loops=.TRUE.
     IF (interruptor2.eqv..TRUE.) THEN
        symm%length=ll
        symm%num_crits=new_num_crits
        DEALLOCATE(symm%crit)
        ALLOCATE(symm%crit(ll))
        symm%crit(:)=new_crit(1:ll)
     END IF
  END IF

END SUBROUTINE quito_falsos_bonded_symm

SUBROUTINE quito_falsos_ats_symm(shell)

  IMPLICIT NONE
  
  TYPE(p_shell),INTENT(INOUT)::shell
 
  INTEGER::ii,jj,gg,ggg,kk,ll
  INTEGER::new_num_crits
  INTEGER,DIMENSION(symm%length)::new_crit
  LOGICAL::interruptor,interruptor2

  new_num_crits=0
  ll=0

  interruptor2=.FALSE.
  gg=0
  DO ii=1,symm%num_crits
     gg=gg+1
     interruptor=.FALSE.
     ggg=gg
     DO jj=1,symm%crit(gg)
        ggg=ggg+1
        kk=symm%crit(gg+1)
        kk=shell%ats_ord(kk)
        interruptor=interruptor.OR.shell%ats(kk)%hbs%with_loops_2order
        interruptor=interruptor.OR.shell%ats(kk)%bs%with_loops_2order
        interruptor=interruptor.OR.shell%ats(kk)%hbs%with_loops_3order
        interruptor=interruptor.OR.shell%ats(kk)%bs%with_loops_3order
     END DO
     IF (interruptor.eqv..FALSE.) THEN
        interruptor2=.TRUE.
        gg=gg+symm%crit(gg)
     ELSE
        new_num_crits=new_num_crits+1
        ll=ll+1
        new_crit(ll)=symm%crit(gg)
        DO jj=1,symm%crit(gg)
           ll=ll+1
           gg=gg+1
           new_crit(ll)=symm%crit(gg)
        END DO
     END IF
  END DO

  IF (new_num_crits==0) THEN
     symm%num_crits=0
  ELSE
     shell%unsolved_loops=.TRUE.
     IF (interruptor2.eqv..TRUE.) THEN
        symm%length=ll
        symm%num_crits=new_num_crits
        DEALLOCATE(symm%crit)
        ALLOCATE(symm%crit(ll))
        symm%crit(:)=new_crit(1:ll)
     END IF
  END IF

END SUBROUTINE quito_falsos_ats_symm


SUBROUTINE quito_falsos_bonded_symm2(bonded,shell2nd,nnods)

  INTEGER,INTENT(IN)::nnods
  TYPE(p_bonded),INTENT(INOUT)::bonded
  TYPE(p_shell),DIMENSION(nnods),TARGET,INTENT(IN)::shell2nd

  INTEGER::ii,jj,gg,kk,ll,mm,pp,ggg
  INTEGER::new_num_crits
  INTEGER,DIMENSION(symm%length)::new_crit
  LOGICAL::interruptor,interruptor2,interruptor3

  new_num_crits=0
  ll=0
  interruptor2=.FALSE.
  gg=0
  DO ii=1,symm%num_crits
     gg=gg+1
     kk=symm%crit(gg+1)
     interruptor=.FALSE.
     IF (bonded%cant_loops_1order(kk)==0) THEN
        IF (bonded%cant_loops_2order(kk)==0) THEN
           interruptor=.TRUE.
        END IF
     END IF
     IF (interruptor.eqv..TRUE.) THEN
        interruptor3=.FALSE.
        ggg=gg
        DO jj=1,symm%crit(gg)
           ggg=ggg+1
           kk=bonded%order(symm%crit(ggg))
           mm=bonded%aux2sh(kk)
           interruptor3=interruptor3.OR.shell2nd(mm)%unsolved_loops
           DO pp=1,shell2nd(mm)%nats
              interruptor3=interruptor3.OR.shell2nd(mm)%ats(pp)%hbs%with_loops_3order
              interruptor3=interruptor3.OR.shell2nd(mm)%ats(pp)%bs%with_loops_3order
              interruptor3=interruptor3.OR.shell2nd(mm)%ats(pp)%hbs%with_loops_2order
              interruptor3=interruptor3.OR.shell2nd(mm)%ats(pp)%bs%with_loops_2order
           END DO
        END DO
        IF (interruptor3.eqv..TRUE.) THEN
           interruptor=.FALSE.
        END IF
        IF (interruptor.eqv..TRUE.) THEN
           interruptor2=.TRUE.
           gg=gg+symm%crit(gg)
        ELSE
           new_num_crits=new_num_crits+1
           ll=ll+1
           new_crit(ll)=symm%crit(gg)
           DO jj=1,symm%crit(gg)
              ll=ll+1
              gg=gg+1
              new_crit(ll)=symm%crit(gg)
           END DO
        END IF
     END IF
  END DO

  IF (new_num_crits==0) THEN
     symm%num_crits=0
  ELSE
     bonded%unsolved_loops=.TRUE.
     IF (interruptor2.eqv..TRUE.) THEN
        print*,'AQUIIIIIIIIIIIIIIIII 2'
        symm%length=ll
        symm%num_crits=new_num_crits
        DEALLOCATE(symm%crit)
        ALLOCATE(symm%crit(ll))
        symm%crit(:)=new_crit(1:ll)
     END IF
  END IF

END SUBROUTINE quito_falsos_bonded_symm2


SUBROUTINE quito_falsos_ats_symm2(shell1st,shell2nd,nnods)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::nnods
  TYPE(p_shell),INTENT(INOUT)::shell1st
  TYPE(p_shell),DIMENSION(nnods),TARGET,INTENT(IN)::shell2nd

 
  INTEGER::ii,jj,gg,ggg,kk,ll,pp,qq,mm
  INTEGER::new_num_crits
  INTEGER,DIMENSION(symm%length)::new_crit
  LOGICAL::interruptor,interruptor2

  new_num_crits=0
  ll=0

  interruptor2=.FALSE.
  gg=0
  DO ii=1,symm%num_crits
     gg=gg+1
     interruptor=.FALSE.
     ggg=gg
     DO jj=1,symm%crit(gg)
        ggg=ggg+1
        kk=symm%crit(gg+1)
        kk=shell1st%ats_ord(kk)
        interruptor=interruptor.OR.shell1st%ats(kk)%hbs%with_loops_2order
        DO pp=1,shell1st%ats(kk)%hbs%num
           mm=shell1st%ats(kk)%hbs%aux2sh(pp)
           interruptor=interruptor.OR.shell2nd(mm)%unsolved_loops
           DO qq=1,shell2nd(mm)%nats
              interruptor=interruptor.OR.shell2nd(mm)%ats(qq)%hbs%with_loops_3order
              interruptor=interruptor.OR.shell2nd(mm)%ats(qq)%bs%with_loops_3order
              interruptor=interruptor.OR.shell2nd(mm)%ats(qq)%hbs%with_loops_2order
              interruptor=interruptor.OR.shell2nd(mm)%ats(qq)%bs%with_loops_2order
           END DO
        END DO
        interruptor=interruptor.OR.shell1st%ats(kk)%bs%with_loops_2order
        DO pp=1,shell1st%ats(kk)%bs%num
           mm=shell1st%ats(kk)%bs%aux2sh(pp)
           interruptor=interruptor.OR.shell2nd(mm)%unsolved_loops
           DO qq=1,shell2nd(mm)%nats
              interruptor=interruptor.OR.shell2nd(mm)%ats(qq)%hbs%with_loops_3order
              interruptor=interruptor.OR.shell2nd(mm)%ats(qq)%bs%with_loops_3order
              interruptor=interruptor.OR.shell2nd(mm)%ats(qq)%hbs%with_loops_2order
              interruptor=interruptor.OR.shell2nd(mm)%ats(qq)%bs%with_loops_2order
           END DO
        END DO
     END DO
     IF (interruptor.eqv..FALSE.) THEN
        interruptor2=.TRUE.
        gg=gg+symm%crit(gg)
     ELSE
        new_num_crits=new_num_crits+1
        ll=ll+1
        new_crit(ll)=symm%crit(gg)
        DO jj=1,symm%crit(gg)
           ll=ll+1
           gg=gg+1
           new_crit(ll)=symm%crit(gg)
        END DO
     END IF
  END DO

  IF (new_num_crits==0) THEN
     symm%num_crits=0
  ELSE
     shell1st%unsolved_loops=.TRUE.
     IF (interruptor2.eqv..TRUE.) THEN
        symm%length=ll
        symm%num_crits=new_num_crits
        DEALLOCATE(symm%crit)
        ALLOCATE(symm%crit(ll))
        symm%crit(:)=new_crit(1:ll)
     END IF
  END IF

END SUBROUTINE quito_falsos_ats_symm2


SUBROUTINE SOLVING_LAST_SYMMETRIES_PERMUTING (shell1st,shell2nd,nnods,totntot)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::nnods,totntot
  TYPE(p_shell),INTENT(INOUT)::shell1st
  TYPE(p_shell),DIMENSION(nnods),TARGET,INTENT(INOUT)::shell2nd

  INTEGER::total_num_permutations
  INTEGER::ii,jj,iii
  INTEGER::gg,kk,ll,hh,extra_dim,kkk,ggg

  INTEGER,DIMENSION(:),ALLOCATABLE::num_permutaciones
  INTEGER,DIMENSION(:,:),ALLOCATABLE::dale

  INTEGER::ntot1sh,nnods1sh,ntot2sh,nnods2sh
  INTEGER,DIMENSION(:),ALLOCATABLE::aux_ind_ats,aux_ind_nods,aux_ord,aux_symm

  TYPE(p_shell),POINTER::shell_aux


  !1st shell

  total_num_permutations=1
  extra_dim=shell1st%nats2+1
  DO ii=1,nnods
     extra_dim=extra_dim+shell2nd(ii)%nats2+1
  END DO
  ALLOCATE(num_permutaciones(extra_dim))
  hh=0

  IF (shell1st%unsolved_loops.eqv..TRUE.) THEN
     CALL upload_symm(shell1st%symm)
     CALL permuto (shell1st%ats_ord,shell1st%nats,shell1st%permutations,shell1st%num_permutations)
     CALL download_symm(shell1st%symm)
  ELSE
     shell1st%num_permutations=1
  END IF
  total_num_permutations=total_num_permutations*shell1st%num_permutations
  hh=hh+1
  num_permutaciones(hh)=shell1st%num_permutations

  DO ii=1,shell1st%nats
     IF (shell1st%ats(ii)%hbs%unsolved_loops) THEN
        CALL upload_symm(shell1st%ats(ii)%hbs%symm)
        CALL permuto (shell1st%ats(ii)%hbs%order,shell1st%ats(ii)%hbs%num,shell1st%ats(ii)%hbs%permutations,shell1st%ats(ii)%hbs%num_permutations)
        CALL download_symm(shell1st%ats(ii)%hbs%symm)
     ELSE
        shell1st%ats(ii)%hbs%num_permutations=1
     END IF
     total_num_permutations=total_num_permutations*shell1st%ats(ii)%hbs%num_permutations
     hh=hh+1
     num_permutaciones(hh)=shell1st%ats(ii)%hbs%num_permutations
     IF (shell1st%ats(ii)%bs%unsolved_loops) THEN
        CALL upload_symm(shell1st%ats(ii)%bs%symm)
        CALL permuto (shell1st%ats(ii)%bs%order,shell1st%ats(ii)%bs%num,shell1st%ats(ii)%bs%permutations,shell1st%ats(ii)%bs%num_permutations)
        CALL download_symm(shell1st%ats(ii)%bs%symm)
     ELSE
        shell1st%ats(ii)%bs%num_permutations=1
     END IF
     total_num_permutations=total_num_permutations*shell1st%ats(ii)%bs%num_permutations
     hh=hh+1
     num_permutaciones(hh)=shell1st%ats(ii)%bs%num_permutations
  END DO
  
  DO jj=1,nnods
     IF (shell2nd(jj)%unsolved_loops.eqv..TRUE.) THEN
        CALL upload_symm(shell2nd(jj)%symm)
        CALL permuto (shell2nd(jj)%ats_ord,shell2nd(jj)%nats,shell2nd(jj)%permutations,shell2nd(jj)%num_permutations)
        CALL download_symm(shell2nd(jj)%symm)
     ELSE
        shell2nd(jj)%num_permutations=1
     END IF
     total_num_permutations=total_num_permutations*shell2nd(jj)%num_permutations
     hh=hh+1
     num_permutaciones(hh)=shell2nd(jj)%num_permutations
     DO ii=1,shell2nd(jj)%nats
        IF (shell2nd(jj)%ats(ii)%hbs%unsolved_loops) THEN
           CALL upload_symm(shell2nd(jj)%ats(ii)%hbs%symm)
           CALL permuto (shell2nd(jj)%ats(ii)%hbs%order,shell2nd(jj)%ats(ii)%hbs%num,shell2nd(jj)%ats(ii)%hbs%permutations,shell2nd(jj)%ats(ii)%hbs%num_permutations)
           CALL download_symm(shell2nd(jj)%ats(ii)%hbs%symm)
        ELSE
           shell2nd(jj)%ats(ii)%hbs%num_permutations=1
        END IF
        total_num_permutations=total_num_permutations*shell2nd(jj)%ats(ii)%hbs%num_permutations
        hh=hh+1
        num_permutaciones(hh)=shell2nd(jj)%ats(ii)%hbs%num_permutations
        IF (shell2nd(jj)%ats(ii)%bs%unsolved_loops) THEN
           CALL upload_symm(shell2nd(jj)%ats(ii)%bs%symm)
           CALL permuto (shell2nd(jj)%ats(ii)%bs%order,shell2nd(jj)%ats(ii)%bs%num,shell2nd(jj)%ats(ii)%bs%permutations,shell2nd(jj)%ats(ii)%bs%num_permutations)
           CALL download_symm(shell2nd(jj)%ats(ii)%bs%symm)
        ELSE
           shell2nd(jj)%ats(ii)%bs%num_permutations=1
        END IF
        total_num_permutations=total_num_permutations*shell2nd(jj)%ats(ii)%bs%num_permutations
        hh=hh+1
        num_permutaciones(hh)=shell2nd(jj)%ats(ii)%bs%num_permutations
     END DO
  END DO


  IF (total_num_permutations>1) THEN

     ALLOCATE(dale(total_num_permutations,extra_dim))
     DO ii=1,extra_dim
        kkk=1
        DO kk=ii+1,extra_dim
           kkk=kkk*num_permutaciones(kk)
        END DO
        ll=1
        DO WHILE (ll<=total_num_permutations)
           DO jj=1,num_permutaciones(ii)
              DO kk=1,kkk
                 dale(ll,ii)=jj
                 ll=ll+1
              END DO
           END DO
        END DO
     END DO

     DO iii=1,total_num_permutations
        PRINT*,dale(iii,:)
        print*,'--'
     END DO

     DO iii=1,total_num_permutations

        hh=0
        hh=hh+1
        IF (shell1st%unsolved_loops.eqv..TRUE.) THEN
           shell1st%ats_ord(:)=shell1st%permutations(dale(iii,hh),:)
        END IF
        
        DO ii=1,shell1st%nats
           hh=hh+1
           IF (shell1st%ats(ii)%hbs%unsolved_loops) THEN
              print*,'si',dale(iii,hh)
              print*,shell1st%ats(ii)%hbs%permutations(dale(iii,hh),:)
              shell1st%ats(ii)%hbs%order(:)=shell1st%ats(ii)%hbs%permutations(dale(iii,hh),:)
              CALL REORDER_BONDED (shell1st%ats(ii)%hbs)
           END IF
           hh=hh+1
           IF (shell1st%ats(ii)%bs%unsolved_loops) THEN !! cagada
              shell1st%ats(ii)%bs%order(:)=shell1st%ats(ii)%bs%permutations(dale(iii,hh),:)
              CALL REORDER_BONDED (shell1st%ats(ii)%bs)
           END IF
        END DO

        DO jj=1,nnods
           hh=hh+1
           IF (shell2nd(jj)%unsolved_loops.eqv..TRUE.) THEN
              shell2nd(jj)%ats_ord(:)=shell2nd(jj)%permutations(dale(iii,hh),:)
           END IF
           DO ii=1,shell2nd(jj)%nats
              hh=hh+1
              IF (shell2nd(jj)%ats(ii)%hbs%unsolved_loops) THEN
                 shell2nd(jj)%ats(ii)%hbs%order(:)=shell2nd(jj)%ats(ii)%hbs%permutations(dale(iii,hh),:)
                 CALL REORDER_BONDED (shell2nd(jj)%ats(ii)%hbs)
              END IF
              hh=hh+1
              IF (shell2nd(jj)%ats(ii)%bs%unsolved_loops) THEN
                 shell2nd(jj)%ats(ii)%bs%order(:)=shell2nd(jj)%ats(ii)%bs%permutations(dale(iii,hh),:)
                 CALL REORDER_BONDED (shell2nd(jj)%ats(ii)%bs)
              END IF
           END DO
        END DO

        IF (ALLOCATED(mss_ind_ats))   DEALLOCATE(mss_ind_ats)
        IF (ALLOCATED(mss_ind_nods))  DEALLOCATE(mss_ind_nods)
        IF (ALLOCATED(mss_symm))      DEALLOCATE(mss_symm)
        

        ALLOCATE(mss_ind_ats(totntot),mss_ind_nods(totntot),mss_symm(totntot))
        
        gg=0
        ntot1sh=shell1st%ntot
        nnods1sh=shell1st%nnods
        ALLOCATE(aux_ind_ats(ntot1sh),aux_ind_nods(ntot1sh),aux_symm(ntot1sh),aux_ord(nnods1sh))
        CALL build_mss_wout_word(shell1st,aux_ind_ats,aux_ind_nods,aux_symm,aux_ord)
        ggg=gg+1
        gg=gg+ntot1sh
        mss_ind_ats(ggg:gg)=aux_ind_ats(:)
        mss_ind_nods(ggg:gg)=aux_ind_nods(:)
        mss_symm(ggg:gg)=aux_symm(:)
        DEALLOCATE(aux_ind_ats,aux_ind_nods,aux_symm)
        
        DO ii=1,nnods1sh
           jj=aux_ord(ii)
           shell_aux=>shell2nd(jj)
           ntot2sh =shell_aux%ntot
           nnods2sh=shell_aux%nnods
           ALLOCATE(aux_ind_ats(ntot2sh),aux_ind_nods(ntot2sh),aux_symm(ntot2sh))
           CALL build_mss_wout(shell_aux,aux_ind_ats,aux_ind_nods,aux_symm)
           ggg=gg+1
           gg=gg+ntot2sh
           mss_ind_ats(ggg:gg)=aux_ind_ats(:)
           mss_ind_nods(ggg:gg)=aux_ind_nods(:)
           mss_symm(ggg:gg)=aux_symm(:)
           DEALLOCATE(aux_ind_ats,aux_ind_nods,aux_symm)
        END DO
        DEALLOCATE(aux_ord)

        CALL remato_mss2sh()
        print*,'>>'
        print'(85I3)',mss(:)
 
     END DO
     NULLIFY(shell_aux)
     print*,'&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
  END IF


END SUBROUTINE SOLVING_LAST_SYMMETRIES_PERMUTING


SUBROUTINE PERMUTO (orden_original,dim_orden,permutations,num_permutations)

  IMPLICIT NONE

  TYPE int_2d_pointer
     INTEGER,DIMENSION(:,:),ALLOCATABLE::p2d
  END type int_2d_pointer

  INTEGER,INTENT(IN)::dim_orden
  INTEGER,DIMENSION(dim_orden),INTENT(IN)::orden_original
  INTEGER,INTENT(INOUT)::num_permutations
  INTEGER,DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT)::permutations

  TYPE(int_2d_pointer),DIMENSION(symm%num_crits)::perms
  INTEGER,DIMENSION(symm%num_crits)::num_perms,num_permats
  INTEGER,ALLOCATABLE,DIMENSION(:,:)::dale

  INTEGER::ii,jj,kk,gg,kkk,ll,cantidad,pos,ind,aa
  INTEGER,DIMENSION(:),ALLOCATABLE::desord


  gg=0
  num_permutations=1
  DO ii=1,symm%num_crits
     gg=gg+1
     cantidad=symm%crit(gg)
     num_perms(ii)=cantidad
     ALLOCATE(desord(cantidad))
     desord(:)=symm%crit((gg+1):(gg+cantidad))
     SELECT CASE (cantidad)
     CASE (2)
        num_perms(ii)=2
        num_permats(ii)=2
        ALLOCATE(perms(ii)%p2d(2,2))
        DO aa=1,2
           perms(ii)%p2d(:,aa)=desord(precalc2(:,aa))
        END DO
        num_permutations=num_permutations*2
     CASE (3)
        num_perms(ii)=6
        num_permats(ii)=3
        ALLOCATE(perms(ii)%p2d(3,6))
        DO aa=1,6
           perms(ii)%p2d(:,aa)=desord(precalc3(:,aa))
        END DO
        num_permutations=num_permutations*6
     CASE (4)
        num_perms(ii)=24
        num_permats(ii)=4
        ALLOCATE(perms(ii)%p2d(4,24))
        DO aa=1,24
           perms(ii)%p2d(:,aa)=desord(precalc4(:,aa))
        END DO
        num_permutations=num_permutations*24
     CASE DEFAULT
        print*,'# Error: there is no precalculated list for', cantidad,' elements permutations.'
        exit
     END SELECT
     DEALLOCATE(desord)
     gg=gg+cantidad
  END DO

  ALLOCATE(permutations(num_permutations,dim_orden))
  ALLOCATE(dale(num_permutations,symm%num_crits))

  DO ii=1,symm%num_crits
     kkk=1
     DO kk=ii+1,symm%num_crits
        kkk=kkk*num_perms(kk)
     END DO
     ll=1
     DO WHILE (ll<=num_permutations)
        DO jj=1,num_perms(ii)
           DO kk=1,kkk
              dale(ll,ii)=jj
              ll=ll+1
           END DO
        END DO
     END DO
  END DO

  DO ii=1,num_permutations
     permutations(ii,:)=orden_original(:)
     DO jj=1,symm%num_crits
        ll=dale(ii,jj)
        DO kk=1,num_permats(jj)
           pos=perms(jj)%p2d(1,kk)
           ind=orden_original(perms(jj)%p2d(ll,kk))
           permutations(ii,pos)=ind
        END DO
     END DO
  END DO


END SUBROUTINE PERMUTO


END MODULE GLOB


