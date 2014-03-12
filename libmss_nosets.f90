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
   INTEGER,DIMENSION(:),ALLOCATABLE::loops_in,loops_ant,loops_post
   LOGICAL::with_loops_in,with_loops_ant,with_loops_post

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
   INTEGER,DIMENSION(:),ALLOCATABLE::ats_ord
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
     ALLOCATE(bonded%loops_post(num))
     ALLOCATE(bonded%loops_ant(num))
     ALLOCATE(bonded%aux2sh(num)) 
     ALLOCATE(bonded%wsymm(num)) 
     ALLOCATE(bonded%symm%crit(num+1)) !!!
     bonded%order                =(/(kk,kk=1,num)/)
     bonded%loops_in(:)          = 0
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
     ALLOCATE(shell%ats_ord(kk))
     ALLOCATE(shell%ats(kk))
     ALLOCATE(shell%wsymm(kk))
     shell%ats_ord(:)=(/(jj,jj=1,kk)/)
     shell%wsymm(:)=0
     DO jj=1,kk
        gg=gg+1
        shell%ats_ord(jj)=jj
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



!!!###################################################################################################
!!!####  REPETITION FUNCTIONS
!!!###################################################################################################



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


END MODULE GLOB


