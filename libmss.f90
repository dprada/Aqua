!!!#################################
!!!####  COMMON VARIABLES AND
!!!####  SUBROUTINES TO UPLOAD THEM
!!!#################################

MODULE GLOB

TYPE p_bonded
   INTEGER,DIMENSION(:),ALLOCATABLE::bonded_ats,bonded_nods
   INTEGER,DIMENSION(:),ALLOCATABLE::bonded_ord
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
   TYPE(p_at),DIMENSION(:),POINTER::ats
   INTEGER,DIMENSION(:),ALLOCATABLE::ats_ord
   INTEGER::nats,nats2,nnods,ntot
   INTEGER::ind
   INTEGER::symm_num_crits
   INTEGER,DIMENSION(:),ALLOCATABLE::symm
   LOGICAL::filtro_supsets
   INTEGER::lev_supsets,lev_sets,lev_nods

END TYPE p_shell

!f2py   intent(hide)::list_shells
TYPE(p_shell),DIMENSION(:),POINTER::list_shells

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

INTEGER,DIMENSION(:),ALLOCATABLE::mss_ind_ats,mss_ind_nods,mss


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
  IF (ASSOCIATED(list_shells))     DEALLOCATE(list_shells) !! cuidado aqui

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
     shell_aux%ats_ord(:)=(/(jj,jj=1,kk)/)
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
        ALLOCATE(at_aux%hbs%bonded_ats(numhbs),at_aux%hbs%bonded_nods(numhbs),at_aux%hbs%bonded_ord(numhbs))
        ALLOCATE(at_aux%hbs%lev_supsets(numhbs),at_aux%hbs%lev_sets(numhbs))
        ALLOCATE(at_aux%hbs%lev_cantsets(numhbs),at_aux%hbs%lev_nods(numhbs))
        ALLOCATE(at_aux%bs%bonded_ats(numbs),at_aux%bs%bonded_nods(numbs),at_aux%bs%bonded_ord(numbs))
        ALLOCATE(at_aux%bs%lev_supsets(numbs),at_aux%bs%lev_sets(numbs))
        ALLOCATE(at_aux%bs%lev_cantsets(numbs),at_aux%bs%lev_nods(numbs))
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
        at_aux%hbs%lev_cantsets(:) = 0
        at_aux%bs%lev_cantsets(:)  = 0
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
  TYPE(p_shell)::shell1st
  TYPE(p_at),POINTER::at_aux

  INTEGER,DIMENSION(:),ALLOCATABLE::privilegios
  INTEGER::dim_privil

  shell1st=list_shells(core)

  !! order bonds

  dim_privil=1
  ALLOCATE(privilegios(1))
  privilegios(1)=core

  DO ii=1,shell1st%nats
     at_aux=>shell1st%at(ii)
       IF (at_aux%hbs%num>1) THEN
          aux_num=at_aux%hbs%num
          ALLOCATE(symm_crit(aux_num))
          num_crits=1
          symm_crit(:)=(/aux_num,(ii,ii=1,aux_num)/)
          CALL order_bonded(dim_privil,privilegios,at_aux%hbs%bonded,num_crits)

       END IF
       IF (at_aux%bs%num>1) THEN
          CALL order_bonded(dim_privil,privilegios,at_aux%bs%bonded)
       END IF

  !CALL order_bonded_1st(shell1st)
  !CALL order_ats_1st(shell1st)
  CALL build_mss_shell1st(shell1st)
 
END SUBROUTINE build_shell1st

SUBROUTINE order_bonded(dim_privil,privilegios,bonded,num_crits)

  INTEGER,INTENT(IN)::dim_privil
  INTEGER,DIMENSION(dim_privil),INTENT(IN)::privilegios
  TYPE(p_bonded),INTENT(INOUT)::bonded
  INTEGER,INTENT(INOUT)::num_crits
  
  INTEGER::ii,jj,gg,priv,aa,bb,numb
  INTEGER,DIMENSION(:,:),ALLOCATABLE::valores_aux,valores
  LOGICAL::sihay

  numb=bonded%num

  IF (superfiltro_supersets.eqv..TRUE.) THEN
     ALLOCATE(valores_aux(dim_privil*2,numb))
  ELSE
     ALLOCATE(valores_aux(dim_privil,numb))
  END IF

  valores_aux(:,:)=0
  gg=0

  IF (superfiltro_supersets.eqv..TRUE.) THEN !! si hay supersets en el sistema
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
                    valores_aux(ggg,hh)=1
                    IF (priv==bonded%bonded_nods(hh)) THEN
                       sihay2=.TRUE.
                       valores_aux(gggg,hh)=1
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
                    valores_aux(ggg,hh)=1
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
                    valores_aux(ggg,hh)=1
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
     ggg=gg+1
     sihay1=.FALSE.
     DO hh=1,numb
        IF (priv==bonded%bonded_nods(hh)) THEN
           sihay1=.TRUE.
           valores_aux(ggg,hh)=1
        END IF
     END DO
     IF (sihay1.eqv..TRUE.) THEN
        gg=ggg
     END IF
  END IF

  IF (bonded%filtro_supsets) THEN
     ALLOCATE(valores(gg+3,numb))
     valores(1:gg,:)=valores_aux(1:gg,:)
     valores(gg+1,:)=bonded%lev_supsets(:)
     valores(gg+2,:)=bonded%lev_cantsets(:)
     valores(gg+3,:)=bonded%lev_nods(:)
     gg=gg+3
     DEALLOCATE(valores_aux)
  ELSE
     ALLOCATE(valores(gg+1,numb))
     valores(1:gg,:)=valores_aux(1:gg,:)
     valores(gg+1,:)=bonded%lev_supsets(:)
     gg=gg+1
  END IF

  CALL SORT_INT_MATRIX (numb,gg,bonded%order,valores,num_crits)


END SUBROUTINE order_bonded

SUBROUTINE build_mss_shell1st(shell1st)
 
  TYPE(p_shell),INTENT(IN)::shell1st
 
  INTEGER::ntot,nats,nats2,nnods
  INTEGER::ii,jj,kk
  INTEGER,DIMENSION(:),ALLOCATABLE::oo
  
 
  IF (ALLOCATED(mss_ind_ats))   DEALLOCATE(mss_ind_ats)
  IF (ALLOCATED(mss_ind_nods))  DEALLOCATE(mss_ind_nods)
  IF (ALLOCATED(mss))           DEALLOCATE(mss)
 
  ntot = shell1st%ntot
  nats = shell1st%nats
  nats2 = shell1st%nats2
  nnods = shell1st%nnods
 
  ALLOCATE(oo(nats))
  ALLOCATE(mss_ind_ats(ntot),mss_ind_nods(ntot),mss(ntot))
 
  mss(:)=0
 
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
 
END SUBROUTINE build_mss_shell1st
 
!!!!!!!!!!!!!!!!!!!!##########

SUBROUTINE SORT_INT_MATRIX (num_ats,dim_matrix,order,valores,num_crits)

  IMPLICIT NONE
 
  INTEGER,INTENT(IN)::num_ats,dim_matrix
  INTEGER,DIMENSION(num_ats),INTENT(INOUT)::order
  INTEGER,DIMENSION(num_ats,dim_matrix),INTENT(IN)::valores
  INTEGER,INTENT(INOUT)::num_crits
 
  pp=1
  DO WHILE ((num_crits>0).AND.(pp<=dim_matrix))

     gg=0
     DO ii=1,num_crits
        gg=gg+1
        idim=symm_crit(gg)
        ALLOCATE(val_aux(idim),ind_aux(idim),symm_aux(idim))
        DO jj=1,idim
           kk=symm_crit(gg+jj)
           ll=order(kk)
           val_aux(jj)=valores(ll,pp)
           symm_aux(jj)=kk
           ind_aux(jj)=ll
        END DO
        IF (ALL(val_aux(2:).eq.val_aux(1)).eqv..FALSE.) THEN
           ALLOCATE(filtro(idim),vals(idim),inds(idim))
           filtro=.TRUE.
           DO jj=1,idim
              kk=MAXLOC(val_aux,DIM=1,MASK=filtro(:))
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
                       cajon(:)=new_symm_aux(:)
                       DEALLOCATE(new_symm_aux)
                       ALLOCATE(new_symm_aux(tope+1+ll))
                       new_symm_aux(1:tope)=cajon(:)
                       tope=tope+1
                       new_symm_aux(tope)=ll
                       new_symm_aux((tope+1):(tope+ll))=order_aux(mm:(mm+ll-1))
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
                 cajon(:)=new_symm_aux(:)
                 DEALLOCATE(new_symm_aux)
                 ALLOCATE(new_symm_aux(tope+1+ll))
                 new_symm_aux(1:tope)=cajon(:)
                 tope=tope+1
                 new_symm_aux(tope)=ll
                 new_symm_aux((tope+1):(tope+ll))=order_aux(mm:(mm+ll-1))
                 tope=tope+ll
                 DEALLOCATE(cajon)
              END IF
           END IF
           DEALLOCATE(filtro,vals,inds)
        END IF
        DEALLOCATE(val_aux,ind_aux,symm_aux)
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

END SUBROUTINE SORT_INT_MATRIX

SUBROUTINE SORT_INT_MATRIX_ATS (num_ats,dim_matrix,order,valores,num_crits)
 
  IMPLICIT NONE
 
  INTEGER,INTENT(IN)::num_ats,dim_matrix
  INTEGER,DIMENSION(num_ats),INTENT(INOUT)::order
  INTEGER,DIMENSION(num_ats,dim_matrix),INTENT(IN)::valores
  INTEGER,INTENT(INOUT)::num_crits
 
  INTEGER::ii,jj,kk,ll,gg,nn,pp,idim,new_num_crits,tope
  INTEGER,DIMENSION(:),ALLOCATABLE::val_aux,ind_aux,order_aux
  INTEGER,DIMENSION(:),ALLOCATABLE::vals,inds
  INTEGER,DIMENSION(:),ALLOCATABLE::new_symm_aux,cajon
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro
  LOGICAL::interruptor
 
  pp=1
  DO WHILE ((num_crits>0).AND.(pp<=dim_matrix))
 
     new_num_crits=0
     
     gg=0
     DO ii=1,num_crits
        gg=gg+1
        idim=symm_crit(gg)
        ALLOCATE(val_aux(idim))
        DO jj=1,idim
           val_aux(jj)=valores(symm_crit(gg+jj),pp)
        END DO
        IF (ALL(val_aux(2:).eq.val_aux(1)).eqv..TRUE.) THEN !! si no entra no debería hacer nada
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





!!!#################################
!!!####  BUILD FIRST SHELL
!!!#################################





!!$  INTEGER::s1_nats,s1_nats2,s1_nnods,s1_ntot,s1_ats_symm_num_crits
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::s1_ats_symm,s1_order_ats,s1_num_bonds,s1_order_bonded_ats
!!$
!!$
!!$
!!$  ALLOCATE(s1_num_bonds(s1_nats2))
!!$  ALLOCATE(s1_order_bonded_ats(s1_nnods))
!!$  kk=1
!!$  ll=0
!!$  DO ii=1,s1_nats
!!$     jj=s1_order_ats(ii)
!!$     aa=num_Hbs_at(jj)
!!$     bb=num_Bs_at(jj)
!!$     cc=aa+bb
!!$     dd=in_Hbs(jj)
!!$     ee=in_Bs(jj)
!!$     s1_num_bonds(kk)=aa
!!$     s1_num_bonds(kk+1)=bb
!!$     kk=kk+2
!!$     s1_order_bonded_ats((ll+1):(ll+cc))=(/Hbs(dd+1:dd+aa),Bs(ee+1:ee+bb)/)
!!$     ll=ll+cc
!!$  END DO
!!$
!!$  IF (ALLOCATED(mss_ind_ats))   DEALLOCATE(mss_ind_ats)
!!$  IF (ALLOCATED(mss_ind_nods))  DEALLOCATE(mss_ind_nods)
!!$  IF (ALLOCATED(mss))           DEALLOCATE(mss)
!!$
!!$  ALLOCATE(mss_ind_ats(s1_ntot),mss_ind_nods(s1_ntot),mss(s1_ntot))
!!$  mss_ind_ats(:)=0
!!$  mss_ind_nods(:)=0
!!$  mss(:)=0
!!$
!!$  jj=1
!!$  mss_ind_ats(jj)  = s1_nats
!!$  mss_ind_nods(jj) = s1_nats
!!$  ii=jj+1
!!$  jj=jj+s1_nats
!!$  mss_ind_ats(ii:jj)  = s1_order_ats(:)
!!$  mss_ind_nods(ii:jj) = at2nod(s1_order_ats(:))
!!$  mss_ind_ats(ii:jj)  = trad2py_at(s1_order_ats(:))
!!$  mss_ind_nods(ii:jj) = trad2py_nod(mss_ind_nods(ii:jj))
!!$  ii=jj+1
!!$  jj=jj+s1_nats2
!!$  mss_ind_ats(ii:jj)  = s1_num_bonds(:)
!!$  mss_ind_nods(ii:jj) = s1_num_bonds(:)
!!$  ii=jj+1
!!$  jj=jj+s1_nnods
!!$  mss_ind_ats(ii:jj)  = s1_order_bonded_ats(:)
!!$  mss_ind_nods(ii:jj) = at2nod(s1_order_bonded_ats(:))
!!$  mss_ind_ats(ii:jj)  = trad2py_at(s1_order_bonded_ats(:))
!!$  mss_ind_nods(ii:jj) = trad2py_nod(mss_ind_nods(ii:jj))
!!$
!!$  DEALLOCATE(s1_ats_symm,s1_order_ats,s1_num_bonds,s1_order_bonded_ats)
!!$
!!$END SUBROUTINE build_shell1st


!!$SUBROUTINE build_shell2nd (core)
!!$
!!$  IMPLICIT NONE
!!$
!!$  TYPE iarray_pointer
!!$     INTEGER,DIMENSION(:),POINTER::p1
!!$  END TYPE iarray_pointer
!!$
!!$  INTEGER,INTENT(IN)::core
!!$
!!$  INTEGER::aa,bb,cc,dd,ee,ii,jj,kk,ll,mm,nn,gg,iii,jjj,aaa,bbb,ccc,ggg
!!$  INTEGER::core2
!!$
!!$  INTEGER::s1_nats,s1_nats2,s1_nnods,s1_ntot,s1_ats_symm_num_crits
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::s1_ats_symm,s1_order_ats,s1_num_bonds,s1_order_bonded_ats
!!$
!!$  INTEGER::si_nats,si_nats2,si_nnods
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::s2_nats,s2_nats2,s2_nnods,s2_ntot,s2_ats_symm_num_crits
!!$  TYPE(iarray_pointer),DIMENSION(:),POINTER::s2_ats_symm,s2_order_ats,s2_num_bonds,s2_order_bonded_ats
!!$
!!$  INTEGER::num_crits,dim_symm_crits,dim_mat
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::aux_order,aux_nod,aux_at
!!$  INTEGER,DIMENSION(:,:),ALLOCATABLE::valores,box
!!$  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro,filtro2
!!$
!!$  INTEGER::ntot
!!$
!!$  s1_nats=num_ats_nod(core)
!!$  s1_nats2=num_ats_nod2(core)
!!$  s1_nnods=num_Hbs_Bs_nod(core)
!!$  s1_ntot=s1_nats+s1_nats2+s1_nnods+1
!!$  ll=ats_symm_length(core)
!!$  ii=ats_symm_in(core)
!!$  jj=ii+ll
!!$  ALLOCATE(s1_ats_symm(ll))
!!$  s1_ats_symm_num_crits=ats_symm_num_crits(core)
!!$  s1_ats_symm(:)=ats_symm(ii:jj)
!!$  ALLOCATE(s1_order_ats(s1_nats))
!!$  ii=in_at_nod(core)
!!$  s1_order_ats(:)=(/(jj,jj=1,s1_nats)/)+ii
!!$  ALLOCATE(s1_num_bonds(s1_nats2))
!!$  ALLOCATE(s1_order_bonded_ats(s1_nnods))
!!$  kk=1
!!$  ll=0
!!$  DO ii=1,s1_nats
!!$     jj=s1_order_ats(ii)
!!$     aa=num_Hbs_at(jj)
!!$     bb=num_Bs_at(jj)
!!$     cc=aa+bb
!!$     dd=in_Hbs(jj)
!!$     ee=in_Bs(jj)
!!$     s1_num_bonds(kk)=aa
!!$     s1_num_bonds(kk+1)=bb
!!$     kk=kk+2
!!$     s1_order_bonded_ats((ll+1):(ll+cc))=(/Hbs(dd+1:dd+aa),Bs(ee+1:ee+bb)/)
!!$     ll=ll+cc
!!$  END DO
!!$
!!$  ALLOCATE(s2_nats(s1_nnods),s2_nats2(s1_nnods),s2_nnods(s1_nnods),s2_ntot(s1_nnods),s2_ats_symm_num_crits(s1_nnods))
!!$  ALLOCATE(s2_ats_symm(s1_nnods),s2_order_ats(s1_nnods),s2_num_bonds(s1_nnods),s2_order_bonded_ats(s1_nnods))
!!$
!!$  DO iii=1,s1_nnods
!!$     core2=s1_order_bonded_ats(iii)
!!$     core2=at2nod(core2)
!!$     si_nats=num_ats_nod(core2)
!!$     si_nats2=num_ats_nod2(core2)
!!$     si_nnods=num_Hbs_Bs_nod(core2)
!!$     s2_nats(iii)=si_nats
!!$     s2_nats2(iii)=si_nats2
!!$     s2_nnods(iii)=si_nnods
!!$     s2_ntot(iii)=si_nats+si_nats2+si_nnods+1
!!$     ll=ats_symm_length(core)
!!$     ii=ats_symm_in(core)
!!$     jj=ii+ll
!!$     ALLOCATE(s2_ats_symm(iii)%p1(ll))
!!$     s2_ats_symm_num_crits(iii)=ats_symm_num_crits(core2)
!!$     s2_ats_symm(iii)%p1(:)=ats_symm(ii:jj)
!!$     ALLOCATE(s2_order_ats(iii)%p1(si_nats))
!!$     ii=in_at_nod(core2)
!!$     s2_order_ats(iii)%p1(:)=(/(jj,jj=1,si_nats)/)+ii
!!$     ALLOCATE(s2_num_bonds(iii)%p1(si_nats2))
!!$     ALLOCATE(s2_order_bonded_ats(iii)%p1(si_nnods))
!!$     kk=1
!!$     ll=0
!!$     DO ii=1,si_nats
!!$        jj=s2_order_ats(iii)%p1(ii)
!!$        aa=num_Hbs_at(jj)
!!$        bb=num_Bs_at(jj)
!!$        cc=aa+bb
!!$        dd=in_Hbs(jj)
!!$        ee=in_Bs(jj)
!!$        s2_num_bonds(iii)%p1(kk)=aa
!!$        s2_num_bonds(iii)%p1(kk+1)=bb
!!$        kk=kk+2
!!$        s2_order_bonded_ats(iii)%p1((ll+1):(ll+cc))=(/Hbs(dd+1:dd+aa),Bs(ee+1:ee+bb)/)
!!$        ll=ll+cc
!!$     END DO
!!$  END DO
!!$
!!$
!!$  !!! Quito simetrias en 2nd shell order bonded
!!$
!!$  DO iii=1,s1_nnods
!!$
!!$     jj=0
!!$     DO jjj=1,s2_nats2(iii)
!!$        si_nats=s2_num_bonds(iii)%p1(jjj)
!!$        ii=jj+1
!!$        jj=jj+si_nats
!!$        IF (si_nats>1) THEN
!!$           ALLOCATE(aux_at(si_nats))
!!$           aux_at(:)=s2_order_bonded_ats(iii)%p1(ii:jj)
!!$           CALL bueno(core,aux_at,si_nats)
!!$           s2_order_bonded_ats(iii)%p1(ii:jj)=aux_at(:)
!!$
!!$           mm=ii
!!$           DO nn=1,si_nats              
!!$              s2_order_bonded_ats(iii)%p1(mm)=aux_at(aux_order(nn))
!!$              mm=mm+1
!!$           END DO
!!$           DEALLOCATE(aux_order,aux_nod,aux_at,filtro)
!!$           IF (ALLOCATED(symm_crit)) DEALLOCATE(symm_crit)
!!$        END IF
!!$     END DO
!!$
!!$  END DO
!!$
!!$  !!! Quito simetrias en 2nd shell order ats (Y si lo engancho con el hecho de no tener antes armado el microestado, puedo ahorrar calculo)
!!$  !!! Deberia armar el microestado justo despues de este paso.
!!$
!!$  !DO iii=1,s1_nnods
!!$
!!$
!!$  !END DO
!!$
!!$  ! Quito simetrias order bonded 1st
!!$
!!$  ! Quito simetrias order ats 1st
!!$
!!$
!!$  ntot=s1_ntot+SUM(s2_ntot(:),DIM=1)
!!$
!!$  IF (ALLOCATED(mss_ind_ats))   DEALLOCATE(mss_ind_ats)
!!$  IF (ALLOCATED(mss_ind_nods))  DEALLOCATE(mss_ind_nods)
!!$  IF (ALLOCATED(mss))           DEALLOCATE(mss)
!!$
!!$  ALLOCATE(mss_ind_ats(ntot),mss_ind_nods(ntot),mss(ntot))
!!$  mss_ind_ats(:)=0
!!$  mss_ind_nods(:)=0
!!$  mss(:)=0
!!$
!!$  jj=1
!!$  mss_ind_ats(jj)  = s1_nats
!!$  mss_ind_nods(jj) = s1_nats
!!$  ii=jj+1
!!$  jj=jj+s1_nats
!!$  mss_ind_ats(ii:jj)  = s1_order_ats(:)
!!$  mss_ind_nods(ii:jj) = at2nod(s1_order_ats(:))
!!$  mss_ind_ats(ii:jj)  = trad2py_at(s1_order_ats(:))
!!$  mss_ind_nods(ii:jj) = trad2py_nod(mss_ind_nods(ii:jj))
!!$  ii=jj+1
!!$  jj=jj+s1_nats2
!!$  mss_ind_ats(ii:jj)  = s1_num_bonds(:)
!!$  mss_ind_nods(ii:jj) = s1_num_bonds(:)
!!$  ii=jj+1
!!$  jj=jj+s1_nnods
!!$  mss_ind_ats(ii:jj)  = s1_order_bonded_ats(:)
!!$  mss_ind_nods(ii:jj) = at2nod(s1_order_bonded_ats(:))
!!$  mss_ind_ats(ii:jj)  = trad2py_at(s1_order_bonded_ats(:))
!!$  mss_ind_nods(ii:jj) = trad2py_nod(mss_ind_nods(ii:jj))
!!$
!!$  DO iii=1,s1_nnods
!!$     jj=jj+1
!!$     si_nats=s2_nats(iii) 
!!$     mss_ind_ats(jj)  = si_nats
!!$     mss_ind_nods(jj) = si_nats
!!$     ii=jj+1
!!$     jj=jj+si_nats
!!$     mss_ind_ats(ii:jj)  = s2_order_ats(iii)%p1(:)
!!$     mss_ind_nods(ii:jj) = at2nod(mss_ind_ats(ii:jj))
!!$     mss_ind_ats(ii:jj)  = trad2py_at(mss_ind_ats(ii:jj))
!!$     mss_ind_nods(ii:jj) = trad2py_nod(mss_ind_nods(ii:jj))
!!$     ii=jj+1
!!$     jj=jj+s2_nats2(iii)
!!$     mss_ind_ats(ii:jj)  = s2_num_bonds(iii)%p1(:)
!!$     mss_ind_nods(ii:jj) = mss_ind_ats(ii:jj)
!!$     ii=jj+1
!!$     jj=jj+s2_nnods(iii)
!!$     mss_ind_ats(ii:jj)  = s2_order_bonded_ats(iii)%p1(:)
!!$     mss_ind_nods(ii:jj) = at2nod(mss_ind_ats(ii:jj))
!!$     mss_ind_ats(ii:jj)  = trad2py_at(mss_ind_ats(ii:jj))
!!$     mss_ind_nods(ii:jj) = trad2py_nod(mss_ind_nods(ii:jj))
!!$     DEALLOCATE(s2_ats_symm(iii)%p1,s2_order_ats(iii)%p1,s2_num_bonds(iii)%p1,s2_order_bonded_ats(iii)%p1)
!!$  END DO
!!$
!!$
!!$  DEALLOCATE(s1_ats_symm,s1_order_ats,s1_num_bonds,s1_order_bonded_ats)
!!$  DEALLOCATE(s2_nats,s2_nats2,s2_nnods,s2_ntot,s2_ats_symm_num_crits)
!!$  DEALLOCATE(s2_ats_symm,s2_order_ats,s2_num_bonds,s2_order_bonded_ats)
!!$
!!$END SUBROUTINE build_shell2nd
!!$
!!$SUBROUTINE bueno()
!!$
!!$
!!$  num_crits=1
!!$  dim_symm_crits=si_nats+1
!!$  ALLOCATE(aux_order(si_nats),aux_nod(si_nats),aux_at(si_nats),filtro(si_nats),symm_crit(dim_symm_crits))
!!$  aux_order(:)=(/(aa,aa=1,si_nats)/)
!!$  symm_crit(:)=(/si_nats,aux_order(:)/)
!!$  aux_nod(:)=at2nod(aux_at(:))
!!$  IF (ANY(aux_nod(:).eq.core).eqv..TRUE.) THEN ! SI EL CORE ESTA EN ORDER
!!$     IF (superfiltro_supsets.eqv..TRUE.) THEN !SI EXISTEN SUPSETS
!!$        IF (filtro_supsets(core).eqv..TRUE.) THEN !SI EL CORE ESTA EN SUPSETS
!!$           dim_mat=3
!!$           ALLOCATE(valores(dim_symm_crits,dim_mat))
!!$           valores(:,:)=0
!!$           aa=lev_supsets(core)
!!$           bb=lev_sets(core)
!!$           cc=core
!!$           gg=0
!!$           DO nn=1,num_crits
!!$              gg=gg+1
!!$              DO mm=1,symm_crit(gg)
!!$                 gg=gg+1
!!$                 kk=aux_nod(symm_crit(gg))
!!$                 aaa=lev_supsets(kk)
!!$                 bbb=lev_sets(kk)
!!$                 ccc=kk
!!$                 IF (aa==aaa) THEN
!!$                    valores(gg,1)=1
!!$                    IF (bb==bbb) THEN
!!$                       valores(gg,2)=1
!!$                       IF (cc==ccc) THEN
!!$                          valores(gg,3)=1
!!$                       END IF
!!$                    END IF
!!$                 END IF
!!$              END DO
!!$           END DO
!!$           CALL SORT_INT_MATRIX_ATS (si_nats,dim_symm_crits,dim_mat,aux_order,valores,num_crits)
!!$           DEALLOCATE(valores)
!!$        ELSE !SI EL CORE NO ESTA EN SUPSETS
!!$           dim_mat=2
!!$           ALLOCATE(valores(dim_symm_crits,dim_mat))
!!$           valores(:,:)=0
!!$           aa=lev_supsets(core)
!!$           bb=core
!!$           gg=0
!!$           DO nn=1,num_crits
!!$              gg=gg+1
!!$              DO mm=1,symm_crit(gg)
!!$                 gg=gg+1
!!$                 kk=aux_nod(symm_crit(gg))
!!$                 aaa=lev_supsets(kk)
!!$                 bbb=kk
!!$                 IF (aa==aaa) THEN
!!$                    valores(gg,1)=1
!!$                    IF (bb==bbb) THEN
!!$                       valores(gg,2)=1
!!$                    END IF
!!$                 END IF
!!$              END DO
!!$           END DO
!!$           CALL SORT_INT_MATRIX_ATS (si_nats,dim_symm_crits,dim_mat,aux_order,valores,num_crits)
!!$           DEALLOCATE(valores)
!!$        END IF
!!$     ELSE ! SI NO EXISTEN SUPSETS
!!$        dim_mat=2
!!$        ALLOCATE(valores(dim_symm_crits,dim_mat))
!!$        valores(:,:)=0
!!$        aa=lev_supsets(core)
!!$        bb=core
!!$        gg=0
!!$        DO nn=1,num_crits
!!$           gg=gg+1
!!$           DO mm=1,symm_crit(gg)
!!$              gg=gg+1
!!$              kk=aux_nod(symm_crit(gg))
!!$              aaa=lev_supsets(kk)
!!$              bbb=kk
!!$              IF (aa==aaa) THEN
!!$                 valores(gg,1)=1
!!$                 IF (bb==bbb) THEN
!!$                    valores(gg,2)=1
!!$                 END IF
!!$              END IF
!!$           END DO
!!$        END DO
!!$        CALL SORT_INT_MATRIX_ATS (si_nats,dim_symm_crits,dim_mat,aux_order,valores,num_crits)
!!$        DEALLOCATE(valores)
!!$     END IF
!!$  END IF
!!$  !! aqui podria poner mas condiciones como si hay loops en 2nd shell
!!$  IF (num_crits>0) THEN
!!$     IF (superfiltro_supsets.eqv..TRUE.) THEN
!!$        filtro(:)=.FALSE.
!!$        gg=0
!!$        DO nn=1,num_crits
!!$           gg=gg+1
!!$           DO mm=1,symm_crit(gg)
!!$              gg=gg+1
!!$              filtro(symm_crit(gg))=.TRUE.
!!$           END DO
!!$        END DO
!!$        IF (ANY(filtro_supsets(PACK(aux_nod(:),MASK=filtro))).eqv..TRUE.) THEN
!!$           dim_mat=3
!!$           ALLOCATE(valores(dim_symm_crits,dim_mat),filtro2(dim_symm_crits))
!!$           filtro2=.TRUE.
!!$           gg=0
!!$           DO nn=1,num_crits
!!$              gg=gg+1
!!$              filtro2(gg)=.FALSE.
!!$              DO mm=1,symm_crit(gg)
!!$                 gg=gg+1
!!$                 kk=aux_nod(symm_crit(gg))
!!$                 valores(gg,1)=lev_supsets(kk)
!!$                 valores(gg,2)=lev_sets(kk)
!!$                 valores(gg,3)=lev_nods(kk)
!!$              END DO
!!$           END DO
!!$           aa=MAXVAL(valores(:,1),DIM=1,MASK=filtro2)
!!$           bb=MAXVAL(valores(:,2),DIM=1,MASK=filtro2)
!!$           DEALLOCATE(filtro2)
!!$           ALLOCATE(box(aa,bb))
!!$           gg=0
!!$           DO nn=1,num_crits
!!$              gg=gg+1
!!$              box(:,:)=0
!!$              ggg=gg
!!$              DO mm=1,symm_crit(gg)
!!$                 gg=gg+1
!!$                 aa=valores(gg,1)
!!$                 bb=valores(gg,2)
!!$                 box(aa,bb)=box(aa,bb)+1
!!$              END DO
!!$              gg=ggg
!!$              DO mm=1,symm_crit(gg)
!!$                 gg=gg+1
!!$                 aa=valores(gg,1)
!!$                 bb=valores(gg,2)
!!$                 valores(gg,2)=box(aa,bb)
!!$              END DO
!!$           END DO
!!$           CALL SORT_INT_MATRIX_ATS (si_nats,dim_symm_crits,dim_mat,aux_order,valores,num_crits)
!!$           DEALLOCATE(valores,box)
!!$        ELSE
!!$           dim_mat=1
!!$           ALLOCATE(valores(dim_symm_crits,dim_mat))
!!$           gg=0
!!$           DO nn=1,num_crits
!!$              gg=gg+1
!!$              DO mm=1,symm_crit(gg)
!!$                 gg=gg+1
!!$                 kk=aux_nod(symm_crit(gg))
!!$                 valores(gg,1)=lev_supsets(kk)
!!$              END DO
!!$           END DO
!!$           CALL SORT_INT_MATRIX_ATS (si_nats,dim_symm_crits,dim_mat,aux_order,valores,num_crits)
!!$           DEALLOCATE(valores)
!!$        END IF
!!$     ELSE
!!$        dim_mat=1
!!$        ALLOCATE(valores(dim_symm_crits,dim_mat))
!!$        gg=0
!!$        DO nn=1,num_crits
!!$           gg=gg+1
!!$           DO mm=1,symm_crit(gg)
!!$              gg=gg+1
!!$              kk=aux_nod(symm_crit(gg))
!!$              valores(gg,1)=lev_supsets(kk)
!!$           END DO
!!$        END DO
!!$        CALL SORT_INT_MATRIX_ATS (si_nats,dim_symm_crits,dim_mat,aux_order,valores,num_crits)
!!$        DEALLOCATE(valores)
!!$     END IF
!!$  END IF
!!$  !aqui reordeno y recoloco
!!$
!!$  aux_at(:)=aux_at(aux_order(nn))
!!$
!!$  DEALLOCATE(aux_order,aux_nod,filtro)
!!$  IF (ALLOCATED(symm_crit)) DEALLOCATE(symm_crit)
!!$
!!$END SUBROUTINE bueno

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
!!$  IF (superfiltro_supsets.eqv..TRUE.) THEN
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


