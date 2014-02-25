!!!#################################
!!!####  COMMON VARIABLES AND
!!!####  SUBROUTINES TO UPLOAD THEM
!!!#################################

MODULE GLOB

  IMPLICIT NONE

TYPE p_bonded
   INTEGER,DIMENSION(:),ALLOCATABLE::bonded_ats,bonded_nods,aux2sh
   INTEGER,DIMENSION(:),ALLOCATABLE::order,wsymm
   INTEGER,DIMENSION(:),ALLOCATABLE::lev_supsets,lev_sets,lev_nods,lev_cantsets
   LOGICAL::filtro_supsets
   INTEGER::num
END TYPE p_bonded

TYPE p_at
   TYPE(p_bonded)::hbs,bs
   INTEGER::num_hbs_bs
   INTEGER::ind
END TYPE p_at

TYPE p_shell
   TYPE(p_at),DIMENSION(:),ALLOCATABLE::ats
   INTEGER,DIMENSION(:),ALLOCATABLE::ats_ord
   INTEGER::nats,nats2,nnods,ntot
   INTEGER::ind
   INTEGER::symm_num_crits,symm_length
   INTEGER,DIMENSION(:),ALLOCATABLE::symm,wsymm
   LOGICAL::filtro_supsets
   INTEGER::lev_supsets,lev_sets,lev_nods
END TYPE p_shell

!f2py   intent(hide)::list_shells
TYPE(p_shell),DIMENSION(:),ALLOCATABLE,TARGET::list_shells

!! TOPOLOGY

!f2py   intent(hide)::num_ats,num_nods,num_sets,num_supsets
!f2py   intent(hide)::at2nod,lev_nods,lev_sets,lev_supsets
!f2py   intent(hide)::topes_supsets,filtro_supsets,superfiltro_supsets
!f2py   intent(hide)::trad2py_nod,trad2py_at
INTEGER::num_ats,num_nods,num_sets,num_supsets
INTEGER,DIMENSION(:),ALLOCATABLE::at2nod
INTEGER,DIMENSION(:),ALLOCATABLE::lev_nods,lev_sets,lev_supsets
INTEGER,DIMENSION(:),ALLOCATABLE::topes_supsets
LOGICAL,DIMENSION(:),ALLOCATABLE::filtro_supsets
INTEGER,DIMENSION(:),ALLOCATABLE::trad2py_nod,trad2py_at
LOGICAL::superfiltro_supsets


!!! AUXILIAR

!f2py   intent(hide)::symm_crit
INTEGER,DIMENSION(:),ALLOCATABLE::symm_crit

!!! OUTPUT

INTEGER,DIMENSION(:),ALLOCATABLE::mss_ind_ats,mss_ind_nods,mss,mss_symm


CONTAINS

!!!!!!!!!!!!!!!!!!!!##########

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

  INTEGER::gg,hh,ii,jj,kk,ll,mm,nn
  INTEGER,DIMENSION(:),ALLOCATABLE::num_ats_nod
  TYPE(p_shell),POINTER::shell_aux


  num_ats     = xx_num_ats
  num_nods    = xx_num_nods
  num_sets    = xx_num_nods
  num_supsets = xx_num_nods



  IF (ALLOCATED(at2nod))           DEALLOCATE(at2nod)
  IF (ALLOCATED(trad2py_at))       DEALLOCATE(trad2py_at)
  IF (ALLOCATED(trad2py_nod))      DEALLOCATE(trad2py_nod)
  IF (ALLOCATED(list_shells))     DEALLOCATE(list_shells) !! cuidado aqui

  IF (ALLOCATED(lev_nods))         DEALLOCATE(lev_nods)
  IF (ALLOCATED(lev_sets))         DEALLOCATE(lev_sets)
  IF (ALLOCATED(lev_supsets))      DEALLOCATE(lev_supsets)
  IF (ALLOCATED(topes_supsets))    DEALLOCATE(topes_supsets)
  IF (ALLOCATED(filtro_supsets))   DEALLOCATE(filtro_supsets)


  ALLOCATE(at2nod(num_ats),trad2py_at(num_ats),trad2py_nod(num_nods))
  ALLOCATE(lev_nods(num_nods),lev_sets(num_nods),lev_supsets(num_nods))
  ALLOCATE(filtro_supsets(num_nods))
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
     lev_nods(ii)=1
     lev_sets(ii)=1
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
     shell_aux%symm_num_crits=jj
     kk=xx_symm_ats_start(ii)-1
     mm=0
     DO ll=1,jj
        mm=mm+1
        mm=mm+xx_symm_ats(kk+mm)
     END DO
     shell_aux%symm_length=mm
     ALLOCATE(shell_aux%symm(mm))
     shell_aux%symm(:)=xx_symm_ats((kk+1):(kk+mm))
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

  superfiltro_supsets=.FALSE.
  filtro_supsets(:)=.FALSE.
  gg=1
  DO ii=1,xx_symm_sets_num
     superfiltro_supsets=.TRUE.
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

  ! REAJUSTO LEVELS en shells

  DO ii=1,num_nods
     shell_aux=>list_shells(ii)
     shell_aux%filtro_supsets = filtro_supsets(ii)
     shell_aux%lev_supsets    = lev_supsets(ii)
     shell_aux%lev_sets       = lev_sets(ii)
     shell_aux%lev_nods       = lev_nods(ii)
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

  NULLIFY(shell_aux)
  DEALLOCATE(num_ats_nod)

END SUBROUTINE load_topol

SUBROUTINE reset_at(at_aux)

  IMPLICIT NONE

  TYPE(p_at),INTENT(INOUT)::at_aux

  IF (ALLOCATED(at_aux%hbs%bonded_ats)) THEN
     DEALLOCATE(at_aux%hbs%bonded_ats,at_aux%hbs%bonded_nods,at_aux%hbs%order)
     DEALLOCATE(at_aux%hbs%lev_supsets,at_aux%hbs%lev_sets)
     DEALLOCATE(at_aux%hbs%lev_cantsets,at_aux%hbs%lev_nods)
     DEALLOCATE(at_aux%hbs%aux2sh,at_aux%hbs%wsymm)
  END IF

  IF (ALLOCATED(at_aux%bs%bonded_ats)) THEN
     DEALLOCATE(at_aux%bs%bonded_ats,at_aux%bs%bonded_nods,at_aux%bs%order)
     DEALLOCATE(at_aux%bs%lev_supsets,at_aux%bs%lev_sets)
     DEALLOCATE(at_aux%bs%lev_cantsets,at_aux%bs%lev_nods)
     DEALLOCATE(at_aux%bs%aux2sh,at_aux%bs%wsymm)
  END IF

END SUBROUTINE reset_at

SUBROUTINE load_net(xx_hbs,xx_bs,xx_num_Hbs_at,xx_num_Bs_at,&
     xx_Total_num_hbs,xx_Total_num_bs,xx_num_ats)
 
  IMPLICIT NONE
 
  INTEGER,INTENT(IN)::xx_num_ats,xx_Total_num_hbs,xx_Total_num_bs
  INTEGER,DIMENSION(xx_num_ats),INTENT(IN)::xx_num_hbs_at,xx_num_bs_at
  INTEGER,DIMENSION(xx_Total_num_hbs),INTENT(IN)::xx_hbs
  INTEGER,DIMENSION(xx_Total_num_bs),INTENT(IN)::xx_bs
 
  INTEGER::ii,jj,kk,gghb,ggb,ggat,numhbs,numbs,totnum
  INTEGER,DIMENSION(:),ALLOCATABLE::vect_aux_hbs,vect_aux_bs
  TYPE(p_shell),POINTER::shell_aux
  TYPE(p_at),POINTER::at_aux

  gghb=0
  ggb=0
  DO ii=1,num_nods
     shell_aux=>list_shells(ii)
     shell_aux%nnods=0
     DO jj=1,shell_aux%nats
        at_aux=>shell_aux%ats(jj)
        ggat=at_aux%ind
        numhbs=xx_num_hbs_at(ggat)
        numbs=xx_num_bs_at(ggat)
        totnum=numhbs+numbs
        at_aux%num_hbs_bs=totnum
        at_aux%hbs%num=numhbs
        at_aux%bs%num=numbs
        CALL RESET_AT(at_aux)
        ALLOCATE(at_aux%hbs%bonded_ats(numhbs),at_aux%hbs%bonded_nods(numhbs),at_aux%hbs%order(numhbs))
        ALLOCATE(at_aux%hbs%lev_supsets(numhbs),at_aux%hbs%lev_sets(numhbs))
        ALLOCATE(at_aux%hbs%lev_cantsets(numhbs),at_aux%hbs%lev_nods(numhbs))
        ALLOCATE(at_aux%hbs%aux2sh(numhbs),at_aux%hbs%wsymm(numhbs))
        ALLOCATE(at_aux%bs%bonded_ats(numbs),at_aux%bs%bonded_nods(numbs),at_aux%bs%order(numbs))
        ALLOCATE(at_aux%bs%lev_supsets(numbs),at_aux%bs%lev_sets(numbs))
        ALLOCATE(at_aux%bs%lev_cantsets(numbs),at_aux%bs%lev_nods(numbs))
        ALLOCATE(at_aux%bs%aux2sh(numbs),at_aux%bs%wsymm(numbs))
        ALLOCATE(vect_aux_hbs(numhbs),vect_aux_bs(numbs))
        kk=gghb+1
        gghb=gghb+numhbs
        vect_aux_hbs(:)=xx_hbs(kk:gghb)
        kk=ggb+1
        ggb=ggb+numbs
        vect_aux_bs(:)=xx_bs(kk:ggb)
        at_aux%hbs%bonded_ats(:)   = vect_aux_hbs(:)
        at_aux%bs%bonded_ats(:)    = vect_aux_bs(:)
        vect_aux_hbs(:)            = at2nod(vect_aux_hbs(:))
        vect_aux_bs(:)             = at2nod(vect_aux_bs(:))
        at_aux%hbs%bonded_nods(:)  = vect_aux_hbs(:)
        at_aux%bs%bonded_nods(:)   = vect_aux_bs(:)
        at_aux%hbs%lev_supsets(:)  = lev_supsets(vect_aux_hbs(:))
        at_aux%hbs%lev_sets(:)     = lev_sets(vect_aux_hbs(:))
        at_aux%hbs%lev_nods(:)     = lev_nods(vect_aux_hbs(:))
        at_aux%bs%lev_supsets(:)   = lev_supsets(vect_aux_bs(:))
        at_aux%bs%lev_sets(:)      = lev_sets(vect_aux_bs(:))
        at_aux%bs%lev_nods(:)      = lev_nods(vect_aux_bs(:))
        at_aux%hbs%order(:)        = (/(kk,kk=1,numhbs)/)
        at_aux%bs%order(:)         = (/(kk,kk=1,numbs)/)
        at_aux%hbs%lev_cantsets(:) = 0
        at_aux%bs%lev_cantsets(:)  = 0
        at_aux%hbs%aux2sh(:)       = 0
        at_aux%bs%aux2sh(:)        = 0
        at_aux%hbs%wsymm(:)        = 0
        at_aux%bs%wsymm(:)         = 0
        at_aux%hbs%filtro_supsets  = .FALSE.
        at_aux%bs%filtro_supsets   = .FALSE.
        IF (superfiltro_supsets.eqv..TRUE.) THEN
           IF (ANY(filtro_supsets(vect_aux_hbs(:))).eqv..TRUE.) THEN
              at_aux%hbs%filtro_supsets=.TRUE.
              IF (at_aux%hbs%num>1) THEN
                 CALL count_cantsets (at_aux%hbs%num,at_aux%hbs%lev_supsets,at_aux%hbs%lev_sets,at_aux%hbs%lev_cantsets)
              END IF
           END IF
           IF (ANY(filtro_supsets(vect_aux_bs(:))).eqv..TRUE.) THEN
              at_aux%bs%filtro_supsets=.TRUE.
              IF (at_aux%bs%num>1) THEN
                 CALL count_cantsets (at_aux%bs%num,at_aux%bs%lev_supsets,at_aux%bs%lev_sets,at_aux%bs%lev_cantsets)
              END IF
           END IF
        END IF
        DEALLOCATE(vect_aux_hbs,vect_aux_bs)
        shell_aux%nnods=shell_aux%nnods+totnum
     END DO
     shell_aux%ntot=1+shell_aux%nats+shell_aux%nats2+shell_aux%nnods
  END DO

  NULLIFY(shell_aux,at_aux)
  
END SUBROUTINE load_net


SUBROUTINE count_cantsets(dim,vector1,vector2,cant)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::dim
  INTEGER,DIMENSION(dim),INTENT(IN)::vector1,vector2
  INTEGER,DIMENSION(dim),INTENT(OUT)::cant

  INTEGER::ii,jj,kk,gg,ggg,aa,bb
  LOGICAL,DIMENSION(dim)::filtro
  INTEGER,DIMENSION(dim)::cajon

  cant(:)=0
  filtro(:)=.TRUE.

  ii=0
  ggg=0

  DO WHILE (ggg<dim)
     ii=ii+1
     IF (filtro(ii).eqv..TRUE.) THEN
        filtro(ii)=.FALSE.
        gg=1
        cajon(gg)=ii
        aa=vector1(ii)
        bb=vector2(ii)
        DO jj=ii+1,dim
           IF (filtro(jj).eqv..TRUE.) THEN
              IF ((aa==vector1(jj)).AND.(bb==vector2(jj))) THEN
                 filtro(jj)=.FALSE.
                 gg=gg+1
                 cajon(gg)=jj
              END IF
           END IF
        END DO
        DO kk=1,gg
           cant(cajon(kk))=gg
        END DO
        ggg=ggg+gg
     END IF
  END DO

END SUBROUTINE count_cantsets


SUBROUTINE build_shell1st (core)
 
  IMPLICIT NONE
 
  INTEGER,INTENT(IN)::core
  TYPE(p_shell),TARGET::shell1st
  TYPE(p_at),POINTER::at_aux

  INTEGER,DIMENSION(:),ALLOCATABLE::privilegios
  INTEGER::dim_privil,aux_num,ii,num_crits

  shell1st=list_shells(core)
   
  !! order bonds
   
  dim_privil=1
  ALLOCATE(privilegios(1))
  privilegios(1)=core
   
  DO ii=1,shell1st%nats
     at_aux=>shell1st%ats(ii)
     IF (at_aux%hbs%num>1) THEN
        aux_num=at_aux%hbs%num
        ALLOCATE(symm_crit(aux_num+1))
        num_crits=1
        symm_crit(:)=(/aux_num,at_aux%hbs%order(:)/)
        CALL order_bonded(dim_privil,privilegios,at_aux%hbs,num_crits)
        DEALLOCATE(symm_crit)
     END IF
     IF (at_aux%bs%num>1) THEN
        aux_num=at_aux%bs%num
        ALLOCATE(symm_crit(aux_num+1))
        num_crits=1
        symm_crit(:)=(/aux_num,at_aux%bs%order(:)/)
        CALL order_bonded(dim_privil,privilegios,at_aux%bs,num_crits)
        DEALLOCATE(symm_crit)
     END IF
  END DO
   
  !! order ats
   
  IF (shell1st%symm_num_crits>0) THEN
     num_crits=shell1st%symm_num_crits
     ALLOCATE(symm_crit(shell1st%symm_length))
     symm_crit(:)=shell1st%symm(:)
     CALL order_ats(dim_privil,privilegios,shell1st,num_crits)
     shell1st%symm_num_crits=num_crits
     DEALLOCATE(shell1st%symm)
     ii=SIZE(symm_crit)
     ALLOCATE(shell1st%symm(ii))
     shell1st%symm(:)=symm_crit(:)
     shell1st%symm_length=ii
     DEALLOCATE(symm_crit)
  END IF
  
  DEALLOCATE(privilegios)

  CALL build_mss(shell1st)
  CALL encode_mss()
   
  NULLIFY(at_aux)

END SUBROUTINE build_shell1st



SUBROUTINE build_shell2nd (core)
 
  IMPLICIT NONE
 
  INTEGER,INTENT(IN)::core

  TYPE(p_shell),TARGET::shell1st
  TYPE(p_shell),DIMENSION(:),ALLOCATABLE,TARGET::shell2nd

  TYPE(p_at),POINTER::at_aux
  TYPE(p_shell),POINTER::shell_aux

  INTEGER::nnods,totntot,core2
  INTEGER::ntot1sh,nnods1sh,ntot2sh,nnods2sh
  INTEGER::ii,jj,kk,gg,iii,ggg,nn,ll,mm
  INTEGER,DIMENSION(:),ALLOCATABLE::privilegios
  INTEGER::dim_privil,aux_num,num_crits
  INTEGER,DIMENSION(:),ALLOCATABLE::aux_ind_ats,aux_ind_nods,aux_ord,aux_symm


  !! Building the structure

  shell1st=list_shells(core)

  nnods=shell1st%nnods
  totntot=shell1st%ntot

  ALLOCATE(shell2nd(shell1st%nnods))
  gg=0
  DO ii=1,shell1st%nats
     DO jj=1,shell1st%ats(ii)%hbs%num
        gg=gg+1
        shell1st%ats(ii)%hbs%aux2sh(jj)=gg
        kk=shell1st%ats(ii)%hbs%bonded_nods(jj)
        shell2nd(gg)=list_shells(kk)
        totntot=totntot+shell2nd(gg)%ntot
     END DO
     DO jj=1,shell1st%ats(ii)%bs%num
        gg=gg+1
        shell1st%ats(ii)%bs%aux2sh(jj)=gg
        kk=shell1st%ats(ii)%bs%bonded_nods(jj)
        shell2nd(gg)=list_shells(kk)
        totntot=totntot+shell2nd(gg)%ntot
     END DO
  END DO

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
        at_aux=>shell_aux%ats(ii)
        IF (at_aux%hbs%num>1) THEN
           aux_num=at_aux%hbs%num
           ALLOCATE(symm_crit(aux_num+1))
           num_crits=1
           symm_crit(:)=(/aux_num,at_aux%hbs%order(:)/)
           CALL order_bonded(dim_privil,privilegios,at_aux%hbs,num_crits)
           IF (num_crits>0) THEN !! PARA QUITAR
              mm=0
              DO nn=1,num_crits
                 mm=mm+1
                 DO ll=1,symm_crit(mm)
                    mm=mm+1
                    at_aux%hbs%wsymm(symm_crit(mm))=1
                 END DO
              END DO
           END IF
           DEALLOCATE(symm_crit)
        END IF
        IF (at_aux%bs%num>1) THEN
           aux_num=at_aux%bs%num
           ALLOCATE(symm_crit(aux_num+1))
           num_crits=1
           symm_crit(:)=(/aux_num,at_aux%bs%order(:)/)
           CALL order_bonded(dim_privil,privilegios,at_aux%bs,num_crits)
           IF (num_crits>0) THEN !! PARA QUITAR
              mm=0
              DO nn=1,num_crits
                 mm=mm+1
                 DO ll=1,symm_crit(mm)
                    mm=mm+1
                    at_aux%bs%wsymm(symm_crit(mm))=1
                 END DO
              END DO
           END IF
           DEALLOCATE(symm_crit)
        END IF
     END DO

     ! order ats 2sh
     IF (shell_aux%symm_num_crits>0) THEN
        num_crits=shell_aux%symm_num_crits
        ALLOCATE(symm_crit(shell_aux%symm_length))
        symm_crit(:)=shell_aux%symm(:)
        CALL order_ats(dim_privil,privilegios,shell_aux,num_crits)
        shell_aux%symm_num_crits=num_crits
        DEALLOCATE(shell_aux%symm)
        ii=SIZE(symm_crit)
        ALLOCATE(shell_aux%symm(ii))
        shell_aux%symm(:)=symm_crit(:)
        shell_aux%symm_length=ii
        IF (num_crits>0) THEN !! PARA QUITAR
           mm=0
           DO nn=1,num_crits
              mm=mm+1
              DO ll=1,symm_crit(mm)
                 mm=mm+1
                 shell_aux%wsymm(symm_crit(mm))=1
              END DO
           END DO
        END IF
        DEALLOCATE(symm_crit)
     END IF

  END DO

  !! removing symm in 1stsh

  DEALLOCATE(privilegios)
  dim_privil=1
  ALLOCATE(privilegios(1))
  privilegios(1)=core

  ! order bonds 1sh
  DO ii=1,shell1st%nats
     at_aux=>shell1st%ats(ii)
     IF (at_aux%hbs%num>1) THEN
        aux_num=at_aux%hbs%num
        ALLOCATE(symm_crit(aux_num+1))
        num_crits=1
        symm_crit(:)=(/aux_num,at_aux%hbs%order(:)/)
        CALL order_bonded(dim_privil,privilegios,at_aux%hbs,num_crits)
        IF (num_crits>0) THEN
           CALL order_bonded_w_next_shells(at_aux%hbs,shell2nd,nnods,num_crits)
        END IF
        IF (num_crits>0) THEN !! PARA QUITAR
           mm=0
           DO nn=1,num_crits
              mm=mm+1
              DO ll=1,symm_crit(mm)
                 mm=mm+1
                 at_aux%hbs%wsymm(symm_crit(mm))=1
              END DO
           END DO
        END IF
        DEALLOCATE(symm_crit)
     END IF
     IF (at_aux%bs%num>1) THEN
        aux_num=at_aux%bs%num
        ALLOCATE(symm_crit(aux_num+1))
        num_crits=1
        symm_crit(:)=(/aux_num,at_aux%bs%order(:)/)
        CALL order_bonded(dim_privil,privilegios,at_aux%bs,num_crits)
        IF (num_crits>0) THEN
           CALL order_bonded_w_next_shells(at_aux%bs,shell2nd,nnods,num_crits)
        END IF
        IF (num_crits>0) THEN !! PARA QUITAR
           mm=0
           DO nn=1,num_crits
              mm=mm+1
              DO ll=1,symm_crit(mm)
                 mm=mm+1
                 at_aux%bs%wsymm(symm_crit(mm))=1
              END DO
           END DO
        END IF
        DEALLOCATE(symm_crit)
     END IF
  END DO


  ! order ats 1sh
  IF (shell1st%symm_num_crits>0) THEN
     num_crits=shell1st%symm_num_crits
     ALLOCATE(symm_crit(shell1st%symm_length))
     symm_crit(:)=shell1st%symm(:)
     CALL order_ats(dim_privil,privilegios,shell1st,num_crits)
     IF (num_crits>0) THEN 
        CALL order_ats_w_next_shells(shell1st,shell2nd,nnods,num_crits)
     END IF
     shell1st%symm_num_crits=num_crits
     DEALLOCATE(shell1st%symm)
     ii=SIZE(symm_crit)
     ALLOCATE(shell1st%symm(ii))
     shell1st%symm(:)=symm_crit(:)
     shell1st%symm_length=ii
     IF (num_crits>0) THEN !! PARA QUITAR
        mm=0
        DO nn=1,num_crits
           mm=mm+1
           DO ll=1,symm_crit(mm)
              mm=mm+1
              shell1st%wsymm(symm_crit(mm))=1
           END DO
        END DO
     END IF
     DEALLOCATE(symm_crit)     
  END IF

  !! Building the msss
  
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

  !CALL encode_mss()

  DEALLOCATE(aux_ord,shell2nd,privilegios)
  NULLIFY(at_aux,shell_aux)

END SUBROUTINE build_shell2nd





SUBROUTINE order_ats(dim_privil,privilegios,shell,num_crits)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::dim_privil
  INTEGER,DIMENSION(dim_privil),INTENT(IN)::privilegios
  TYPE(p_shell),TARGET,INTENT(INOUT)::shell
  INTEGER,INTENT(INOUT)::num_crits

  INTEGER::ii,jj,kk,gg,ggg,gggg,aa,bb,priv
  INTEGER::numats
  TYPE(p_bonded),POINTER::bonded
  TYPE(p_at),POINTER::at_aux
  INTEGER,DIMENSION(:,:),ALLOCATABLE::valores_aux,valores
  LOGICAL::sihay1,sihay2

  numats=shell%nats

  IF (superfiltro_supsets.eqv..TRUE.) THEN
     ALLOCATE(valores_aux(numats,dim_privil*4))
  ELSE
     ALLOCATE(valores_aux(numats,dim_privil*2))
  END IF

  valores_aux(:,:)=0
  gg=0

  IF (superfiltro_supsets.eqv..TRUE.) THEN !! si hay supersets en el sistema
     DO ii=1,dim_privil
        priv=privilegios(ii)
        IF (filtro_supsets(priv).eqv..TRUE.) THEN
           aa=lev_supsets(priv)
           bb=lev_sets(priv)
           sihay1=.FALSE. !! con hbs
           sihay2=.FALSE.
           ggg=gg+1
           gggg=ggg+1
           DO jj=1,numats
              IF (shell%ats(jj)%hbs%filtro_supsets.eqv..TRUE.) THEN
                 bonded=>shell%ats(jj)%hbs
                 DO kk=1,bonded%num
                    IF ((bonded%lev_supsets(kk)==aa).AND.(bonded%lev_sets(kk)==bb)) THEN
                       sihay1=.TRUE.
                       valores_aux(jj,ggg)=valores_aux(jj,ggg)+1
                       IF (priv==bonded%bonded_nods(kk)) THEN
                          sihay2=.TRUE.
                          valores_aux(jj,gggg)=valores_aux(jj,gggg)+1
                       END IF
                    END IF
                 END DO
              END IF
           END DO
           IF (sihay2.eqv..TRUE.) THEN
              gg=gggg
           ELSE
              IF (sihay1.eqv..TRUE.) THEN
                 gg=ggg
              END IF
           END IF
           sihay1=.FALSE. !! con bs
           sihay2=.FALSE.
           ggg=gg+1
           gggg=ggg+1
           DO jj=1,numats
              IF (shell%ats(jj)%bs%filtro_supsets.eqv..TRUE.) THEN
                 bonded=>shell%ats(jj)%bs
                 DO kk=1,bonded%num
                    IF ((bonded%lev_supsets(kk)==aa).AND.(bonded%lev_sets(kk)==bb)) THEN
                       sihay1=.TRUE.
                       valores_aux(jj,ggg)=valores_aux(jj,ggg)+1
                       IF (priv==bonded%bonded_nods(kk)) THEN
                          sihay2=.TRUE.
                          valores_aux(jj,gggg)=valores_aux(jj,gggg)+1
                       END IF
                    END IF
                 END DO
              END IF
           END DO
           IF (sihay2.eqv..TRUE.) THEN
              gg=gggg
           ELSE
              IF (sihay1.eqv..TRUE.) THEN
                 gg=ggg
              END IF
           END IF
        ELSE
           sihay1=.FALSE. ! con hbs
           ggg=gg+1
           DO jj=1,numats
              bonded=>shell%ats(jj)%hbs
              DO kk=1,bonded%num
                 IF (priv==bonded%bonded_nods(kk)) THEN
                    sihay1=.TRUE.
                    valores_aux(jj,ggg)=valores_aux(jj,ggg)+1
                 END IF
              END DO
           END DO
           IF (sihay1.eqv..TRUE.) THEN
              gg=ggg
           END IF
           sihay1=.FALSE. ! con hbs
           ggg=gg+1
           DO jj=1,numats
              bonded=>shell%ats(jj)%bs
              DO kk=1,bonded%num
                 IF (priv==bonded%bonded_nods(kk)) THEN
                    sihay1=.TRUE.
                    valores_aux(jj,ggg)=valores_aux(jj,ggg)+1
                 END IF
              END DO
           END DO
           IF (sihay1.eqv..TRUE.) THEN
              gg=ggg
           END IF
        END IF
     END DO
  ELSE
     DO ii=1,dim_privil
        priv=privilegios(ii)
        sihay1=.FALSE. ! con hbs
        ggg=gg+1
        DO jj=1,numats
           bonded=>shell%ats(jj)%hbs
           DO kk=1,bonded%num
              IF (priv==bonded%bonded_nods(kk)) THEN
                 sihay1=.TRUE.
                 valores_aux(jj,ggg)=valores_aux(jj,ggg)+1
              END IF
           END DO
        END DO
        IF (sihay1.eqv..TRUE.) THEN
           gg=ggg
        END IF
        sihay1=.FALSE. ! con hbs
        ggg=gg+1
        DO jj=1,numats
           bonded=>shell%ats(jj)%bs
           DO kk=1,bonded%num
              IF (priv==bonded%bonded_nods(kk)) THEN
                 sihay1=.TRUE.
                 valores_aux(jj,ggg)=valores_aux(jj,ggg)+1
              END IF
           END DO
        END DO
        IF (sihay1.eqv..TRUE.) THEN
           gg=ggg
        END IF
     END DO
  END IF

  ALLOCATE(valores(numats,gg+2))
  valores(:,1:gg)=valores_aux(:,1:gg)
  DO jj=1,numats
     at_aux=>shell%ats(jj)
     valores(jj,gg+1)=at_aux%hbs%num
     valores(jj,gg+2)=at_aux%bs%num
  END DO
  gg=gg+2

  DEALLOCATE(valores_aux)

  CALL SORT_INT_MATRIX (numats,gg,shell%ats_ord,valores,num_crits)

  DEALLOCATE(valores)

  IF (num_crits>0) THEN !los ligados a hbs
     gg=0
     DO ii=1,numats
        IF (gg<shell%ats(ii)%hbs%num) gg=shell%ats(ii)%hbs%num
     END DO
     ALLOCATE(valores(numats,gg))
     valores(:,:)=0
     DO ii=1,numats
        valores(ii,1:shell%ats(ii)%hbs%num)=shell%ats(ii)%hbs%lev_supsets(:)
     END DO
     CALL SORT_INT_MATRIX (numats,gg,shell%ats_ord,valores,num_crits)
     IF (superfiltro_supsets.eqv..TRUE.) THEN
        IF (num_crits>0) THEN
           valores(:,:)=0
           DO ii=1,numats
              valores(ii,1:shell%ats(ii)%hbs%num)=shell%ats(ii)%hbs%lev_cantsets(:)
           END DO
           CALL SORT_INT_MATRIX (numats,gg,shell%ats_ord,valores,num_crits)
           IF (num_crits>0) THEN
              valores(:,:)=0
              DO ii=1,numats
                 valores(ii,1:shell%ats(ii)%hbs%num)=shell%ats(ii)%hbs%lev_nods(:)
              END DO
              CALL SORT_INT_MATRIX (numats,gg,shell%ats_ord,valores,num_crits)
           END IF
        END IF
     END IF
     DEALLOCATE(valores)
  END IF

  IF (num_crits>0) THEN !los ligados a bs
     gg=0
     DO ii=1,numats
        IF (gg<shell%ats(ii)%bs%num) gg=shell%ats(ii)%bs%num
     END DO
     ALLOCATE(valores(numats,gg))
     valores(:,:)=0
     DO ii=1,numats
        valores(ii,1:shell%ats(ii)%bs%num)=shell%ats(ii)%bs%lev_supsets(:)
     END DO
     CALL SORT_INT_MATRIX (numats,gg,shell%ats_ord,valores,num_crits)
     IF (superfiltro_supsets.eqv..TRUE.) THEN
        IF (num_crits>0) THEN
           valores(:,:)=0
           DO ii=1,numats
              valores(ii,1:shell%ats(ii)%bs%num)=shell%ats(ii)%bs%lev_cantsets(:)
           END DO
           CALL SORT_INT_MATRIX (numats,gg,shell%ats_ord,valores,num_crits)
           IF (num_crits>0) THEN
              valores(:,:)=0
              DO ii=1,numats
                 valores(ii,1:shell%ats(ii)%bs%num)=shell%ats(ii)%bs%lev_nods(:)
              END DO
              CALL SORT_INT_MATRIX (numats,gg,shell%ats_ord,valores,num_crits)
           END IF
        END IF
     END IF
     DEALLOCATE(valores)
  END IF

  NULLIFY(bonded,at_aux)

END SUBROUTINE order_ats


SUBROUTINE order_ats_w_next_shells(shell1st,shell2nd,nnods,num_crits)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::nnods
  TYPE(p_shell),TARGET,INTENT(INOUT)::shell1st
  TYPE(p_shell),DIMENSION(nnods),TARGET,INTENT(IN)::shell2nd
  INTEGER,INTENT(INOUT)::num_crits

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
  DO ii=1,num_crits
     gg=gg+1
     DO jj=1,symm_crit(gg)
        gg=gg+1
        kk=aux_ord(symm_crit(gg))
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
  DO ii=1,num_crits
     gg=gg+1
     DO jj=1,symm_crit(gg)
        gg=gg+1
        kk=aux_ord(symm_crit(gg))
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

  CALL SORT_INT_MATRIX (numats,dim_matrix,shell1st%ats_ord,valores,num_crits)
  DEALLOCATE(valores,aux_ord)


  IF (num_crits>0) THEN

     dim_matrix=0
     ALLOCATE(aux_ord(numats))
     aux_ord=shell1st%ats_ord(:)

     gg=0
     DO ii=1,num_crits
        gg=gg+1
        DO jj=1,symm_crit(gg)
           gg=gg+1
           kk=aux_ord(symm_crit(gg))
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
        DO ii=1,num_crits
           gg=gg+1
           DO jj=1,symm_crit(gg)
              gg=gg+1
              kk=aux_ord(symm_crit(gg))
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

        CALL SORT_INT_MATRIX (numats,dim_matrix,shell1st%ats_ord,valores,num_crits)

        DEALLOCATE(valores)

        IF ((superfiltro_supsets.eqv..TRUE.).and.(num_crits>0)) THEN

           
           ALLOCATE(valores(numats,dim_matrix*2))
           valores(:,:)=0
           aux_ord=shell1st%ats_ord(:)

           gg=0
           DO ii=1,num_crits
              gg=gg+1
              DO jj=1,symm_crit(gg)
                 gg=gg+1
                 kk=aux_ord(symm_crit(gg))
                 at_aux=>shell1st%ats(kk)
                 aa=0
                 DO ll=1,at_aux%hbs%num
                    mm=at_aux%hbs%aux2sh(ll)
                    shell_aux=>shell2nd(mm)
                    DO nn=1,shell_aux%nats
                       qq=shell_aux%ats_ord(nn)
                       aaa=aa+1
                       aa=aa+shell_aux%ats(qq)%hbs%num
                       valores(kk,aaa:aa)=shell_aux%ats(qq)%hbs%lev_cantsets(:)
                       aaa=aa+1
                       aa=aa+shell_aux%ats(qq)%bs%num
                       valores(kk,aaa:aa)=shell_aux%ats(qq)%bs%lev_cantsets(:)
                    END DO
                 END DO
                 DO ll=1,at_aux%bs%num
                    mm=at_aux%bs%aux2sh(ll)
                    shell_aux=>shell2nd(mm)
                    DO nn=1,shell_aux%nats
                       qq=shell_aux%ats_ord(nn)
                       aaa=aa+1
                       aa=aa+shell_aux%ats(qq)%hbs%num
                       valores(kk,aaa:aa)=shell_aux%ats(qq)%hbs%lev_cantsets(:)
                       aaa=aa+1
                       aa=aa+shell_aux%ats(qq)%bs%num
                       valores(kk,aaa:aa)=shell_aux%ats(qq)%bs%lev_cantsets(:)
                    END DO
                 END DO
                 DO ll=1,at_aux%hbs%num
                    mm=at_aux%hbs%aux2sh(ll)
                    shell_aux=>shell2nd(mm)
                    DO nn=1,shell_aux%nats
                       qq=shell_aux%ats_ord(nn)
                       aaa=aa+1
                       aa=aa+shell_aux%ats(qq)%hbs%num
                       valores(kk,aaa:aa)=shell_aux%ats(qq)%hbs%lev_nods(:)
                       aaa=aa+1
                       aa=aa+shell_aux%ats(qq)%bs%num
                       valores(kk,aaa:aa)=shell_aux%ats(qq)%bs%lev_nods(:)
                    END DO
                 END DO
                 DO ll=1,at_aux%bs%num
                    mm=at_aux%bs%aux2sh(ll)
                    shell_aux=>shell2nd(mm)
                    DO nn=1,shell_aux%nats
                       qq=shell_aux%ats_ord(nn)
                       aaa=aa+1
                       aa=aa+shell_aux%ats(qq)%hbs%num
                       valores(kk,aaa:aa)=shell_aux%ats(qq)%hbs%lev_nods(:)
                       aaa=aa+1
                       aa=aa+shell_aux%ats(qq)%bs%num
                       valores(kk,aaa:aa)=shell_aux%ats(qq)%bs%lev_nods(:)
                    END DO
                 END DO
              END DO
           END DO
           CALL SORT_INT_MATRIX (numats,dim_matrix*2,shell1st%ats_ord,valores,num_crits)
           
           DEALLOCATE(valores)

        END IF

     END IF

     DEALLOCATE(aux_ord)

  END IF


END SUBROUTINE order_ats_w_next_shells

SUBROUTINE order_bonded(dim_privil,privilegios,bonded,num_crits)
  
  IMPLICIT NONE

  INTEGER,INTENT(IN)::dim_privil
  INTEGER,DIMENSION(dim_privil),INTENT(IN)::privilegios
  TYPE(p_bonded),INTENT(INOUT)::bonded
  INTEGER,INTENT(INOUT)::num_crits
  
  INTEGER::ii,gg,hh,priv,aa,bb,numb,ggg,gggg
  INTEGER,DIMENSION(:,:),ALLOCATABLE::valores_aux,valores
  LOGICAL::sihay1,sihay2

  numb=bonded%num
   
  IF (superfiltro_supsets.eqv..TRUE.) THEN
     ALLOCATE(valores_aux(numb,dim_privil*2))
  ELSE
     ALLOCATE(valores_aux(numb,dim_privil))
  END IF
   
   
  valores_aux(:,:)=0
  gg=0
   
  IF (superfiltro_supsets.eqv..TRUE.) THEN !! si hay supersets en el sistema
     IF (bonded%filtro_supsets.eqv..TRUE.) THEN
        DO ii=1,dim_privil
           priv=privilegios(ii)
           IF (filtro_supsets(priv).eqv..TRUE.) THEN !! si hay ssets en bonded y priv es ssets
              aa=lev_supsets(priv)
              bb=lev_sets(priv)
              sihay1=.FALSE.
              sihay2=.FALSE.
              ggg=gg+1
              gggg=ggg+1
              DO hh=1,numb
                 IF ((bonded%lev_supsets(hh)==aa).AND.(bonded%lev_sets(hh)==bb)) THEN
                    sihay1=.TRUE.
                    valores_aux(hh,ggg)=1
                    IF (priv==bonded%bonded_nods(hh)) THEN
                       sihay2=.TRUE.
                       valores_aux(hh,gggg)=1
                    END IF
                 END IF
              END DO
              IF (sihay2.eqv..TRUE.) THEN
                 gg=gggg
              ELSE
                 IF (sihay1.eqv..TRUE.) THEN
                    gg=ggg
                 END IF
              END IF
           ELSE !! si hay ssets en bonded y priv no es ssets
              ggg=gg+1
              sihay1=.FALSE.
              DO hh=1,numb
                 IF (priv==bonded%bonded_nods(hh)) THEN
                    sihay1=.TRUE.
                    valores_aux(hh,ggg)=1
                 END IF
              END DO
              IF (sihay1.eqv..TRUE.) THEN
                 gg=ggg
              END IF
           END IF
        END DO
     ELSE        ! si no hay ssets en bonded 
        DO ii=1,dim_privil
           priv=privilegios(ii)
           IF (filtro_supsets(priv).eqv..FALSE.) THEN
              ggg=gg+1
              sihay1=.FALSE.
              DO hh=1,numb
                 IF (priv==bonded%bonded_nods(hh)) THEN
                    sihay1=.TRUE.
                    valores_aux(hh,ggg)=1
                 END IF
              END DO
              IF (sihay1.eqv..TRUE.) THEN
                 gg=ggg
              END IF
           END IF
        END DO
     END IF
  ELSE
     !! si no hay supersets en el sistema
     DO ii=1,dim_privil
        priv=privilegios(ii)
        ggg=gg+1
        sihay1=.FALSE.
        DO hh=1,numb
           IF (priv==bonded%bonded_nods(hh)) THEN
              sihay1=.TRUE.
              valores_aux(hh,ggg)=1
           END IF
        END DO
        IF (sihay1.eqv..TRUE.) THEN
           gg=ggg
        END IF
     END DO
  END IF
   
  IF (bonded%filtro_supsets.eqv..TRUE.) THEN
     ALLOCATE(valores(numb,gg+3))
     valores(:,1:gg)=valores_aux(:,1:gg)
     valores(:,gg+1)=bonded%lev_supsets(:)
     valores(:,gg+2)=bonded%lev_cantsets(:)
     valores(:,gg+3)=bonded%lev_nods(:)
     gg=gg+3
  ELSE
     ALLOCATE(valores(numb,gg+1))
     valores(:,1:gg)=valores_aux(:,1:gg)
     valores(:,gg+1)=bonded%lev_supsets(:)
     gg=gg+1
  END IF
   
  DEALLOCATE(valores_aux)
   
  CALL SORT_INT_MATRIX (numb,gg,bonded%order,valores,num_crits)
   
  CALL REORDER_BONDED (bonded)
   
  DEALLOCATE(valores)

END SUBROUTINE order_bonded


SUBROUTINE order_bonded_w_next_shells(bonded,shell2nd,nnods,num_crits)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::nnods
  TYPE(p_bonded),INTENT(INOUT)::bonded
  TYPE(p_shell),DIMENSION(nnods),TARGET,INTENT(IN)::shell2nd
  INTEGER,INTENT(INOUT)::num_crits

  INTEGER::numb,dim_matrix
  INTEGER::ii,jj,gg,mm,kk,aa,qq,aaa,nn
  INTEGER,DIMENSION(:,:),ALLOCATABLE::valores
  TYPE(p_shell),POINTER::shell_aux

  dim_matrix=0
  numb=bonded%num

  gg=0
  DO ii=1,num_crits
     gg=gg+1
     DO jj=1,symm_crit(gg)
        gg=gg+1
        kk=bonded%order(symm_crit(gg))
        mm=bonded%aux2sh(kk)
        IF (dim_matrix<shell2nd(mm)%nats2) dim_matrix=shell2nd(mm)%nats2
     END DO
  END DO

  ALLOCATE(valores(numb,dim_matrix))
  valores(:,:)=0

  gg=0
  DO ii=1,num_crits
     gg=gg+1
     DO jj=1,symm_crit(gg)
        gg=gg+1
        kk=bonded%order(symm_crit(gg))
        mm=bonded%aux2sh(kk)
        shell_aux=>shell2nd(mm)
        aa=0
        DO nn=1,shell_aux%nats
           qq=shell_aux%ats_ord(nn)
           aa=aa+1
           valores(kk,aa)=shell_aux%ats(qq)%hbs%num
           aa=aa+1
           valores(kk,aa)=shell_aux%ats(qq)%bs%num
        END DO
     END DO
  END DO

  CALL SORT_INT_MATRIX (numb,dim_matrix,bonded%order,valores,num_crits)
  DEALLOCATE(valores)
  CALL REORDER_BONDED (bonded)

  IF (num_crits>0) THEN

     dim_matrix=0
     numb=bonded%num
     
     gg=0
     DO ii=1,num_crits
        gg=gg+1
        DO jj=1,symm_crit(gg)
           gg=gg+1
           kk=bonded%order(symm_crit(gg))
           mm=bonded%aux2sh(kk)
           IF (dim_matrix<shell2nd(mm)%nnods) dim_matrix=shell2nd(mm)%nnods
        END DO
     END DO
     
     IF (dim_matrix>0) THEN

        ALLOCATE(valores(numb,dim_matrix))
        valores(:,:)=0

        gg=0
        DO ii=1,num_crits
           gg=gg+1
           DO jj=1,symm_crit(gg)
              gg=gg+1
              kk=bonded%order(symm_crit(gg))
              mm=bonded%aux2sh(kk)
              shell_aux=>shell2nd(mm)
              aa=0
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

        CALL SORT_INT_MATRIX (numb,dim_matrix,bonded%order,valores,num_crits)
        DEALLOCATE(valores)
        CALL REORDER_BONDED (bonded)

        IF ((superfiltro_supsets.eqv..TRUE.).and.(num_crits>0)) THEN

           dim_matrix=dim_matrix*2
           ALLOCATE(valores(numb,dim_matrix))
           valores(:,:)=0

           gg=0
           DO ii=1,num_crits
              gg=gg+1
              DO jj=1,symm_crit(gg)
                 gg=gg+1
                 kk=bonded%order(symm_crit(gg))
                 mm=bonded%aux2sh(kk)
                 shell_aux=>shell2nd(mm)
                 aa=0
                 DO nn=1,shell_aux%nats
                    qq=shell_aux%ats_ord(nn)
                    aaa=aa+1
                    aa=aa+shell_aux%ats(qq)%hbs%num
                    valores(kk,aaa:aa)=shell_aux%ats(qq)%hbs%lev_cantsets(:)
                    aaa=aa+1
                    aa=aa+shell_aux%ats(qq)%bs%num
                    valores(kk,aaa:aa)=shell_aux%ats(qq)%bs%lev_cantsets(:)
                    aaa=aa+1
                    aa=aa+shell_aux%ats(qq)%hbs%num
                    valores(kk,aaa:aa)=shell_aux%ats(qq)%hbs%lev_nods(:)
                    aaa=aa+1
                    aa=aa+shell_aux%ats(qq)%bs%num
                    valores(kk,aaa:aa)=shell_aux%ats(qq)%bs%lev_nods(:)
                 END DO
              END DO
           END DO

           
           CALL SORT_INT_MATRIX (numb,dim_matrix,bonded%order,valores,num_crits)
           DEALLOCATE(valores)
           CALL REORDER_BONDED (bonded)

        END IF


     END IF
     
  END IF


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
  IF (bonded%filtro_supsets.eqv..TRUE.) THEN
     box(:)=bonded%lev_sets(:)
     bonded%lev_sets(:)=box(bonded%order(:))
     box(:)=bonded%lev_nods(:)
     bonded%lev_nods(:)=box(bonded%order(:))
     box(:)=bonded%lev_cantsets(:)
     bonded%lev_cantsets(:)=box(bonded%order(:))
  END IF
  DEALLOCATE(box)

  numb=bonded%num
  bonded%order(:)=(/(ii,ii=1,numb)/)

END SUBROUTINE REORDER_BONDED

SUBROUTINE build_mss(shell1st)

  IMPLICIT NONE
 
  TYPE(p_shell),INTENT(IN)::shell1st
 
  INTEGER::ntot,nats,nats2,nnods
  INTEGER::ii,jj,kk
  INTEGER,DIMENSION(:),ALLOCATABLE::oo
  
 
  IF (ALLOCATED(mss_ind_ats))   DEALLOCATE(mss_ind_ats)
  IF (ALLOCATED(mss_ind_nods))  DEALLOCATE(mss_ind_nods)

  ntot = shell1st%ntot
  nats = shell1st%nats
  nats2 = shell1st%nats2
  nnods = shell1st%nnods
 
  ALLOCATE(oo(nats))
  ALLOCATE(mss_ind_ats(ntot),mss_ind_nods(ntot))
 
  oo(:)=shell1st%ats_ord(:)
 
  jj=1
  mss_ind_ats(1)  = nats
  mss_ind_nods(1) = nats
  ii=jj+1
  jj=jj+nats
  mss_ind_ats(ii:jj)  = (/(shell1st%ats(oo(kk))%ind,kk=1,nats)/) !! ojo
  mss_ind_nods(ii:jj) = at2nod(mss_ind_ats(ii:jj))
  mss_ind_ats(ii:jj)  = trad2py_at(mss_ind_ats(ii:jj))
  mss_ind_nods(ii:jj) = trad2py_nod(mss_ind_nods(ii:jj))
  ii=jj+1
  jj=jj+nats2
  mss_ind_ats(ii:jj)  = (/((/shell1st%ats(oo(kk))%hbs%num,shell1st%ats(oo(kk))%bs%num/),kk=1,nats)/)
  mss_ind_nods(ii:jj) = mss_ind_ats(ii:jj)
  ii=jj+1
  jj=jj+nnods
  mss_ind_ats(ii:jj)  = (/((/shell1st%ats(oo(kk))%hbs%bonded_ats(:),shell1st%ats(oo(kk))%bs%bonded_ats(:)/),kk=1,nats)/)
  mss_ind_nods(ii:jj) = at2nod(mss_ind_ats(ii:jj))
  mss_ind_ats(ii:jj)  = trad2py_at(mss_ind_ats(ii:jj))
  mss_ind_nods(ii:jj) = trad2py_nod(mss_ind_nods(ii:jj))
 
  DEALLOCATE(oo)
 
END SUBROUTINE build_mss


SUBROUTINE build_mss_wout(shell2nd,aux_ind_ats,aux_ind_nods,aux_symm)
 
  IMPLICIT NONE

  TYPE(p_shell),INTENT(IN)::shell2nd
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

  aux_symm(:)=0
 
  jj=1
  aux_ind_ats(1)  = nats
  aux_ind_nods(1) = nats
  ii=jj+1
  jj=jj+nats
  aux_symm(ii:jj)     = (/(shell2nd%wsymm(kk),kk=1,nats)/)
  aux_ind_ats(ii:jj)  = (/(shell2nd%ats(oo(kk))%ind,kk=1,nats)/) !! ojo
  aux_ind_nods(ii:jj) = at2nod(aux_ind_ats(ii:jj))
  aux_ind_ats(ii:jj)  = trad2py_at(aux_ind_ats(ii:jj))
  aux_ind_nods(ii:jj) = trad2py_nod(aux_ind_nods(ii:jj))
  ii=jj+1
  jj=jj+nats2
  aux_ind_ats(ii:jj)  = (/((/shell2nd%ats(oo(kk))%hbs%num,shell2nd%ats(oo(kk))%bs%num/),kk=1,nats)/)
  aux_ind_nods(ii:jj) = aux_ind_ats(ii:jj)
  ii=jj+1
  jj=jj+nnods
  aux_symm(ii:jj)     = (/((/shell2nd%ats(oo(kk))%hbs%wsymm(:),shell2nd%ats(oo(kk))%bs%wsymm(:)/),kk=1,nats)/)
  aux_ind_ats(ii:jj)  = (/((/shell2nd%ats(oo(kk))%hbs%bonded_ats(:),shell2nd%ats(oo(kk))%bs%bonded_ats(:)/),kk=1,nats)/)
  aux_ind_nods(ii:jj) = at2nod(aux_ind_ats(ii:jj))
  aux_ind_ats(ii:jj)  = trad2py_at(aux_ind_ats(ii:jj))
  aux_ind_nods(ii:jj) = trad2py_nod(aux_ind_nods(ii:jj))
 
  DEALLOCATE(oo)

END SUBROUTINE build_mss_wout

SUBROUTINE build_mss_wout_word(shell2nd,aux_ind_ats,aux_ind_nods,aux_symm,aux_ord)
  
  IMPLICIT NONE

  TYPE(p_shell),INTENT(IN)::shell2nd
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
 
  aux_symm(:)=0

  jj=1
  aux_ind_ats(1)  = nats
  aux_ind_nods(1) = nats
  ii=jj+1
  jj=jj+nats
  aux_symm(ii:jj)     = (/(shell2nd%wsymm(kk),kk=1,nats)/)
  aux_ind_ats(ii:jj)  = (/(shell2nd%ats(oo(kk))%ind,kk=1,nats)/) !! ojo
  aux_ind_nods(ii:jj) = at2nod(aux_ind_ats(ii:jj))
  aux_ind_ats(ii:jj)  = trad2py_at(aux_ind_ats(ii:jj))
  aux_ind_nods(ii:jj) = trad2py_nod(aux_ind_nods(ii:jj))
  ii=jj+1
  jj=jj+nats2
  aux_ind_ats(ii:jj)  = (/((/shell2nd%ats(oo(kk))%hbs%num,shell2nd%ats(oo(kk))%bs%num/),kk=1,nats)/)
  aux_ind_nods(ii:jj) = aux_ind_ats(ii:jj)
  ii=jj+1
  jj=jj+nnods
  aux_symm(ii:jj)     = (/((/shell2nd%ats(oo(kk))%hbs%wsymm(:),shell2nd%ats(oo(kk))%bs%wsymm(:)/),kk=1,nats)/)
  aux_ord(:)          = (/((/shell2nd%ats(oo(kk))%hbs%bonded_ats(:),shell2nd%ats(oo(kk))%bs%bonded_ats(:)/),kk=1,nats)/)
  aux_ind_ats(ii:jj)  = trad2py_at(aux_ord(:))
  aux_ord(:)          = at2nod(aux_ord(:))
  aux_ind_nods(ii:jj) = aux_ord(:)
  aux_ind_nods(ii:jj) = trad2py_nod(aux_ord(:))
  aux_ord(:)          = (/((/shell2nd%ats(oo(kk))%hbs%aux2sh(:),shell2nd%ats(oo(kk))%bs%aux2sh(:)/),kk=1,nats)/)

  DEALLOCATE(oo)
 
END SUBROUTINE build_mss_wout_word


SUBROUTINE encode_mss()

  IMPLICIT NONE

  INTEGER::ntot

  IF (ALLOCATED(mss))  DEALLOCATE(mss)
  
  ntot=SIZE(mss_ind_nods)

  ALLOCATE(mss(ntot))

  mss(:)=0

END SUBROUTINE encode_mss


!!!!!!!!!!!!!!!!!!!!##########
!!!#################################
!!!####  GENERAL FUNCTIONS TO SORT  
!!!#################################


SUBROUTINE SORT_INT_MATRIX (num_ats,dim_matrix,order,valores,num_crits)

  IMPLICIT NONE
 
  INTEGER,INTENT(IN)::num_ats,dim_matrix
  INTEGER,DIMENSION(num_ats),INTENT(INOUT)::order
  INTEGER,DIMENSION(num_ats,dim_matrix),INTENT(IN)::valores
  INTEGER,INTENT(INOUT)::num_crits
 
  INTEGER::pp,ii,jj,kk,gg,ll,mm,idim,tope,new_num_crits
  INTEGER,DIMENSION(:),ALLOCATABLE::val_aux,ind_aux,symm_aux
  INTEGER,DIMENSION(:),ALLOCATABLE::vals,inds,new_symm_aux,cajon
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro
  LOGICAL::interruptor

  pp=1
  DO WHILE ((num_crits>0).AND.(pp<=dim_matrix))
     new_num_crits=0
     gg=0
     DO ii=1,num_crits
        gg=gg+1
        idim=symm_crit(gg)
        ALLOCATE(val_aux(idim),ind_aux(idim),symm_aux(idim))
        DO jj=1,idim
           gg=gg+1
           kk=symm_crit(gg)
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
     num_crits=new_num_crits
     IF (num_crits>0) THEN
        DEALLOCATE(symm_crit)
        ALLOCATE(symm_crit(tope))
        symm_crit(1:tope)=new_symm_aux(1:tope)
        DEALLOCATE(new_symm_aux)
     END IF
     pp=pp+1
  END DO

END SUBROUTINE SORT_INT_MATRIX



END MODULE GLOB


