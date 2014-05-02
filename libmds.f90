MODULE GLOB

!f2py   intent(hide)::dists
!f2py   intent(hide)::Nnods,Ktot
!f2py   intent(hide)::T_start,T_ind
!f2py   intent(hide)::T_tau
!f2py   intent(hide)::dists_up,net_up
!f2py   intent(hide)::iseed
DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::dists
INTEGER::Nnods,Ktot
INTEGER,DIMENSION(:),ALLOCATABLE::T_start,T_ind
DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::T_tau,Pe
DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::pre_dijkstra
LOGICAL::dists_up=.FALSE.
LOGICAL::net_up=.FALSE.
LOGICAL::pre_dijkstra_up=.FALSE.
INTEGER,DIMENSION(4)::iseed

!f2py   intent(hide)::num_pivots,ind_pivots
INTEGER::num_pivots
INTEGER,DIMENSION(:),ALLOCATABLE::ind_pivots

CONTAINS

SUBROUTINE LOAD_NET (xT_start,xT_ind,xT_tau,xNnods,xKtot)

  INTEGER,INTENT(IN)::xNnods,xKtot
  INTEGER,DIMENSION(xKtot),INTENT(IN)::xT_ind
  DOUBLE PRECISION,DIMENSION(xKtot),INTENT(IN)::xT_tau
  INTEGER,DIMENSION(xNnods+1),INTENT(IN)::xT_start

  INTEGER::ii,jj
  DOUBLE PRECISION::ww,wl

  IF (net_up.eqv..TRUE.) THEN
     DEALLOCATE(T_start,T_ind,T_tau,Pe)
  END IF
  
  Nnods=xNnods
  Ktot=xKtot
  ALLOCATE(T_start(Nnods+1),T_ind(Ktot),T_tau(Ktot),Pe(Nnods))
  T_start(:)=xT_start(:)
  T_ind(:)  =xT_ind(:)
  T_tau(:)  =xT_tau(:)
  Pe(:)     =0.0d0
  net_up    =.TRUE.

  ww=0.0d0
  DO ii=1,Nnods
     DO jj=T_start(ii)+1,T_start(ii+1)
        wl=T_tau(jj)
        ww=ww+wl
        Pe(ii)=Pe(ii)+wl
     END DO
     wl=Pe(ii)
     DO jj=T_start(ii)+1,T_start(ii+1)
        T_tau(jj)=T_tau(jj)/wl
     END DO
  END DO

  Pe(:)=Pe(:)/ww

  iseed(1)=123
  iseed(2)=456
  iseed(3)=2651
  iseed(4)=3455 !impar (limite 4094)

END SUBROUTINE LOAD_NET

SUBROUTINE PRE_RELAX_1ST_ORDER()

  INTEGER::ii,jj,mm
  DOUBLE PRECISION::wl,time_rel

  IF (net_up.eqv..FALSE.) THEN
     PRINT*,'Load net to MDS'
     stop
  END IF

  IF (pre_dijkstra_up.eqv..TRUE.) THEN
     DEALLOCATE(pre_dijkstra)
  END IF

  ALLOCATE(pre_dijkstra(Ktot))
  pre_dijkstra_up=.TRUE.
  pre_dijkstra=0.0d0

  DO ii=1,Nnods
     DO jj=T_start(ii)+1,T_start(ii+1)
        wl=T_tau(jj)
        mm=T_ind(jj)
        IF (mm/=ii) THEN
           time_rel=1.0d0-wl-(Pe(ii)/Pe(mm))*wl
           time_rel=-1.0d0/log(time_rel)
        ELSE
           time_rel=0.0d0
        END IF
        pre_dijkstra(jj)=time_rel
     END DO
  END DO

END SUBROUTINE PRE_RELAX_1ST_ORDER

SUBROUTINE PRE_RELAX_1ST_ORDER2()

  INTEGER::ii,jj,mm
  DOUBLE PRECISION::wl,time_rel

  IF (net_up.eqv..FALSE.) THEN
     PRINT*,'Load net to MDS'
     stop
  END IF

  IF (pre_dijkstra_up.eqv..TRUE.) THEN
     DEALLOCATE(pre_dijkstra)
  END IF

  ALLOCATE(pre_dijkstra(Ktot))
  pre_dijkstra_up=.TRUE.
  pre_dijkstra=0.0d0

  DO ii=1,Nnods
     DO jj=T_start(ii)+1,T_start(ii+1)
        wl=T_tau(jj)
        mm=T_ind(jj)
        IF (mm/=ii) THEN
           time_rel=wl+(Pe(ii)/Pe(mm))*wl
           time_rel=1.0d0/time_rel
        ELSE
           time_rel=0.0d0
        END IF
        pre_dijkstra(jj)=time_rel
     END DO
  END DO

END SUBROUTINE PRE_RELAX_1ST_ORDER2

SUBROUTINE PRE_RELAX_HAMM()

  INTEGER::ii,jj,mm
  DOUBLE PRECISION::wl,time_rel

  IF (net_up.eqv..FALSE.) THEN
     PRINT*,'Load net to MDS'
     stop
  END IF

  IF (pre_dijkstra_up.eqv..TRUE.) THEN
     DEALLOCATE(pre_dijkstra)
  END IF

  ALLOCATE(pre_dijkstra(Ktot))
  pre_dijkstra_up=.TRUE.
  pre_dijkstra=0.0d0

  DO ii=1,Nnods
     DO jj=T_start(ii)+1,T_start(ii+1)
        wl=T_tau(jj)
        mm=T_ind(jj)
        IF (mm/=ii) THEN
           time_rel=(Pe(mm)/wl)*(1.0d0/(Pe(ii)+Pe(mm)))
        ELSE
           time_rel=0.0d0
        END IF
        pre_dijkstra(jj)=time_rel
     END DO
  END DO

END SUBROUTINE PRE_RELAX_HAMM

SUBROUTINE PRE_INV_FLUX()

  INTEGER::ii,jj,mm
  DOUBLE PRECISION::inv_flux

  IF (net_up.eqv..FALSE.) THEN
     PRINT*,'Load net to MDS'
     stop
  END IF

  IF (pre_dijkstra_up.eqv..TRUE.) THEN
     DEALLOCATE(pre_dijkstra)
  END IF

  ALLOCATE(pre_dijkstra(Ktot))
  pre_dijkstra_up=.TRUE.
  pre_dijkstra=0.0d0

  DO ii=1,Nnods
     DO jj=T_start(ii)+1,T_start(ii+1)
        mm=T_ind(jj)
        IF (mm/=ii) THEN
           inv_flux=log(1.0d0/(Pe(ii)*T_tau(jj)))
        ELSE
           inv_flux=0.0d0
        END IF
        pre_dijkstra(jj)=inv_flux
     END DO
  END DO

END SUBROUTINE PRE_INV_FLUX

SUBROUTINE FPT_POINT (inicio,fin,tiempo)

  INTEGER,INTENT(IN)::inicio,fin
  DOUBLE PRECISION,INTENT(OUT)::tiempo

  INTEGER::ii,jj
  DOUBLE PRECISION::dice,bandera
  INTEGER,DIMENSION(4)::hseed

  tiempo=0.0d0

  hseed=iseed

  ii=inicio
  DO WHILE (ii/=fin)
     CALL dlarnv(1,hseed,1,dice)
     bandera=0.0d0
     DO jj=T_start(ii)+1,T_start(ii+1)
        bandera=bandera+T_tau(jj)
        IF (dice<=bandera) THEN
           ii=T_ind(jj)
           EXIT
        END IF
     END DO
     tiempo=tiempo+1.0d0
  END DO

  iseed=hseed

END SUBROUTINE FPT_POINT

SUBROUTINE PRE_MIN_FPT()

  INTEGER::ii,jj,mm,ll,nn,claro
  DOUBLE PRECISION::time_fpt,old_fpt,eps,new_fpt,out_fpt

  IF (net_up.eqv..FALSE.) THEN
     PRINT*,'Load net to MDS'
     stop
  END IF

  IF (pre_dijkstra_up.eqv..TRUE.) THEN
     DEALLOCATE(pre_dijkstra)
  END IF

  ALLOCATE(pre_dijkstra(Ktot))
  pre_dijkstra_up=.TRUE.
  pre_dijkstra=0.0d0

  DO ii=1,Nnods
     DO jj=T_start(ii)+1,T_start(ii+1)
        mm=T_ind(jj)
        IF (mm/=ii) THEN
           time_fpt=0.0
           old_fpt=0.0
           ll=0
           eps=1.0
           claro=0
           DO WHILE (claro/=4)
              CALL fpt_point (ii,mm,out_fpt)
              ll=ll+1
              time_fpt=time_fpt+out_fpt
              new_fpt=(1.0d0*time_fpt)/(1.0d0*ll)
              eps=ABS(new_fpt-old_fpt)
              old_fpt=new_fpt
              IF (eps<0.01) THEN
                 claro=claro+1
              ELSE
                 claro=0
              END IF
           END DO
           time_fpt=new_fpt
        ELSE
           time_fpt=0.0d0
        END IF
        pre_dijkstra(jj)=time_fpt
     END DO
  END DO

  DO ii=1,Nnods
     DO jj=T_start(ii)+1,T_start(ii+1)
        mm=T_ind(jj)
        IF (mm/=ii) THEN
           DO ll=T_start(mm)+1,T_start(mm+1)
              nn=T_ind(ll)
              IF (nn==ii) THEN
                 time_fpt=MIN(pre_dijkstra(jj),pre_dijkstra(ll))
                 pre_dijkstra(jj)=time_fpt
                 pre_dijkstra(ll)=time_fpt
                 EXIT
              END IF
           END DO
        END IF
     END DO
  END DO

END SUBROUTINE PRE_MIN_FPT

SUBROUTINE PRE_AVE_FPT()

  INTEGER::ii,jj,mm,ll,nn,claro
  DOUBLE PRECISION::time_fpt,old_fpt,eps,new_fpt,out_fpt

  IF (net_up.eqv..FALSE.) THEN
     PRINT*,'Load net to MDS'
     stop
  END IF

  IF (pre_dijkstra_up.eqv..TRUE.) THEN
     DEALLOCATE(pre_dijkstra)
  END IF

  ALLOCATE(pre_dijkstra(Ktot))
  pre_dijkstra_up=.TRUE.
  pre_dijkstra=0.0d0

  DO ii=1,Nnods
     DO jj=T_start(ii)+1,T_start(ii+1)
        mm=T_ind(jj)
        IF (mm/=ii) THEN
           time_fpt=0.0
           old_fpt=0.0
           ll=0
           eps=1.0d0
           claro=0
           DO WHILE (claro/=4)
              CALL fpt_point (ii,mm,out_fpt)
              ll=ll+1
              time_fpt=time_fpt+out_fpt
              new_fpt=(1.0d0*time_fpt)/(1.0d0*ll)
              eps=ABS(new_fpt-old_fpt)
              old_fpt=new_fpt
              IF (eps<0.01) THEN
                 claro=claro+1
              ELSE
                 claro=0
              END IF
           END DO
           time_fpt=new_fpt
        ELSE
           time_fpt=0.0d0
        END IF
        pre_dijkstra(jj)=time_fpt
     END DO
  END DO

  DO ii=1,Nnods
     DO jj=T_start(ii)+1,T_start(ii+1)
        mm=T_ind(jj)
        IF (mm/=ii) THEN
           DO ll=T_start(mm)+1,T_start(mm+1)
              nn=T_ind(ll)
              IF (nn==ii) THEN
                 time_fpt=(pre_dijkstra(jj)+pre_dijkstra(ll))/2.0
                 pre_dijkstra(jj)=time_fpt
                 pre_dijkstra(ll)=time_fpt
                 EXIT
              END IF
           END DO
        END IF
     END DO
  END DO

END SUBROUTINE PRE_AVE_FPT



SUBROUTINE DIJKSTRA ()
 
  IMPLICIT NONE

  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::vect_aux
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro,filtro2
  
  DOUBLE PRECISION::azero
  INTEGER::ii,jj,gg,kk

  IF (dists_up.eqv..TRUE.) THEN
     DEALLOCATE(dists)
  END IF

  ALLOCATE(filtro(Nnods),vect_aux(Nnods))
  azero=0.0d0

  ALLOCATE(dists(Nnods,Nnods))
  dists=0.0d0
  dists_up=.TRUE.

  ALLOCATE(filtro2(Nnods))
  DO ii=1,Nnods
     vect_aux=1.0d0/azero
     filtro=.true.
     filtro2=.true.
     vect_aux(ii)=0.0d0
     DO jj=1,ii-1
        vect_aux(jj)=dists(ii,jj)
        filtro2(jj)=.false.
     END DO
     DO WHILE (COUNT(filtro)>0)
        jj=MINLOC(vect_aux(:),DIM=1,MASK=filtro)
        DO kk=T_start(jj)+1,T_start(jj+1)
           gg=T_ind(kk)
           IF (filtro2(gg).eqv..true.) THEN
              IF (vect_aux(gg)>(vect_aux(jj)+pre_dijkstra(kk))) THEN
                 vect_aux(gg)=vect_aux(jj)+pre_dijkstra(kk)
              END IF
           END IF
        END DO
        
        filtro(jj)=.false.
        filtro2(jj)=.false.
     END DO
     dists(:,ii)=vect_aux(:)
  END DO
  DEALLOCATE(filtro2,filtro,vect_aux)
 
END SUBROUTINE DIJKSTRA


SUBROUTINE DIJKSTRA_PIVOTS ()
 
  IMPLICIT NONE

  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::vect_aux
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro,filtro2
  
  DOUBLE PRECISION::azero
  INTEGER::ii,jj,gg,kk,hh

  IF (dists_up.eqv..TRUE.) THEN
     DEALLOCATE(dists)
  END IF

  ALLOCATE(filtro(Nnods),vect_aux(Nnods))
  azero=0.0d0

  ALLOCATE(dists(Nnods,num_pivots))
  dists=0.0d0
  dists_up=.TRUE.

  ALLOCATE(filtro2(Nnods))
  DO hh=1,num_pivots
     ii=ind_pivots(hh)
     vect_aux=1.0d0/azero
     filtro=.true.
     filtro2=.true.
     vect_aux(ii)=0.0d0
     DO WHILE (COUNT(filtro)>0)
        jj=MINLOC(vect_aux(:),DIM=1,MASK=filtro)
        DO kk=T_start(jj)+1,T_start(jj+1)
           gg=T_ind(kk)
           IF (filtro2(gg).eqv..true.) THEN
              IF (vect_aux(gg)>(vect_aux(jj)+pre_dijkstra(kk))) THEN
                 vect_aux(gg)=vect_aux(jj)+pre_dijkstra(kk)
              END IF
           END IF
        END DO
        
        filtro(jj)=.false.
        filtro2(jj)=.false.
     END DO
     dists(:,hh)=vect_aux(:)
  END DO
  DEALLOCATE(filtro2,filtro,vect_aux)
 
END SUBROUTINE DIJKSTRA_PIVOTS

SUBROUTINE DIJKSTRA_A_TARGET (objetivo,xNnods,aux_dists)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::objetivo,xNnods
  DOUBLE PRECISION,DIMENSION(xNnods),INTENT(OUT)::aux_dists
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro,filtro2
  
  DOUBLE PRECISION::azero
  INTEGER::ii,jj,gg,kk,hh


  ALLOCATE(filtro(Nnods))
  azero=0.0d0


  ALLOCATE(filtro2(Nnods))

  ii=objetivo
  aux_dists=1.0d0/azero
  filtro=.true.
  filtro2=.true.
  aux_dists(ii)=0.0d0
  DO WHILE (COUNT(filtro)>0)
     jj=MINLOC(aux_dists(:),DIM=1,MASK=filtro)
     DO kk=T_start(jj)+1,T_start(jj+1)
        gg=T_ind(kk)
        IF (filtro2(gg).eqv..true.) THEN
           IF (aux_dists(gg)>(aux_dists(jj)+pre_dijkstra(kk))) THEN
              aux_dists(gg)=aux_dists(jj)+pre_dijkstra(kk)
           END IF
        END IF
     END DO
     filtro(jj)=.false.
     filtro2(jj)=.false.
  END DO

  DEALLOCATE(filtro2,filtro)
 
END SUBROUTINE DIJKSTRA_A_TARGET


SUBROUTINE MDS (coordinates,eigenvals,eigenvects,stress,opt_stress,dim,lout,xNnods)
 
  IMPLICIT NONE
  
  INTEGER,INTENT(IN)::lout,dim,opt_stress,xNnods
 
  DOUBLE PRECISION,DIMENSION(lout),INTENT(OUT)::eigenvals
  DOUBLE PRECISION,DIMENSION(xNnods,lout),INTENT(OUT)::eigenvects
  DOUBLE PRECISION,DIMENSION(xNnods,dim),INTENT(OUT)::coordinates
  DOUBLE PRECISION,DIMENSION(lout),INTENT(OUT)::stress
 
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::vect_aux,coors_aux
  DOUBLE PRECISION::dd,norm,salida
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::di
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::CC
  INTEGER::num_val,info,Lwork
  INTEGER, DIMENSION(:), ALLOCATABLE::iwork,ifail
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::work
 
  INTEGER::i,j,abajo,ii,jj,gg
  DOUBLE PRECISION::maxL

  eigenvals=0.0d0
  eigenvects=0.0d0
  coordinates=0.0d0
  stress=0.0d0
 
  ALLOCATE(CC(Nnods,Nnods),di(Nnods))
  di=0.0d0
  dd=0.0d0
  CC=0.0d0

  CC(:,:)=dists(:,:)
 
  !! Deberia poner aqui una opcion para hacer un guardado opcional distancias=CC
  !! y usarlas luego como opcion al calcular el stress.
  print*,'calcula matriz'

  DO i=1,Nnods
     DO j=1,Nnods
        CC(j,i)=CC(j,i)**2    ! Necesario en el metodo
        di(i)=di(i)+CC(j,i)
        dd=dd+CC(j,i)
     END DO
     di(i)=di(i)/(Nnods*1.0d0)
  END DO
  dd=dd/(Nnods*Nnods*1.0d0)
 
  DO i=1,Nnods
     DO j=1,Nnods
        CC(j,i)=-0.50d0*(CC(j,i)-di(i)-di(j)+dd)
     END DO
  END DO
 
  Lwork=8*Nnods
 
  ALLOCATE (work(Lwork),iwork(5*Nnods),ifail(Nnods))
 
  eigenvals=0.0d0
  eigenvects=0.0d0
  work=0.0d0
  iwork=0
  ifail=0
  abajo=Nnods-lout+1
 
 
  print*,'va a diagonalizar'
 
  CALL dsyevx ('V','I','U',Nnods,CC,Nnods,0,0,abajo,Nnods,0.0d0,num_val,&
       eigenvals,eigenvects,Nnods,work,Lwork,iwork,ifail,info)
 
  IF (info/=0) THEN
     print*,"#Error with the diagonalization -MDS.f90:dsyevx-"
     print*,"#the array 'work' should has the dimension:",work(1)
  END IF
 
  DEALLOCATE (work,iwork,ifail,CC)
  DEALLOCATE(di)
 
  coordinates=0.0d0
 
  DO j=1,dim
     coordinates(:,j)=sqrt(eigenvals(lout-j+1))*eigenvects(:,lout-j+1)
  END DO
 
  IF (opt_stress==1) THEN
     ALLOCATE(vect_aux(Ktot),coors_aux(Nnods))
     vect_aux=0.0d0
     coors_aux=0.0d0
     norm=0.0d0
     DO ii=1,Nnods
        DO jj=T_start(ii)+1,T_start(ii+1)
           norm=norm+T_tau(jj)**2
        END DO
     END DO
     DO j=1,lout
        if (eigenvals(lout-j+1)>0.0d0) THEN
           salida=0.0d0
           dd=sqrt(eigenvals(lout-j+1))
           coors_aux(:)=dd*eigenvects(:,lout-j+1)
           DO ii=1,Nnods
              DO jj=T_start(ii)+1,T_start(ii+1)
                 gg=T_ind(jj)
                 vect_aux(jj)=vect_aux(jj)+(coors_aux(ii)-coors_aux(gg))**2
                 salida=salida+(sqrt(vect_aux(jj))-T_tau(jj))**2
              END DO
           END DO
           stress(j)=sqrt(salida/norm)
        END IF
     END DO
     DEALLOCATE(vect_aux,coors_aux)
  END IF
  
  maxL=MAXVAL(coordinates)
  maxL=maxL/50.0d0
  coordinates(:,:)=coordinates(:,:)/maxL

 
END SUBROUTINE MDS


SUBROUTINE MDS_PIVOTS (coordinates,dim,xNnods)
!SUBROUTINE MDS_PIVOTS (coordinates,eigenvals,eigenvects,stress,opt_stress,dim,lout,xNnods)
 
  IMPLICIT NONE
  
  INTEGER,INTENT(IN)::dim,xNnods

  DOUBLE PRECISION,DIMENSION(xNnods,dim),INTENT(OUT)::coordinates
 
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::eigenvects_p
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::eigenvals_p
  !DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::vect_aux,coors_aux
  DOUBLE PRECISION::dd !,norm,salida
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::di,dp
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::CC,CCC
  INTEGER::num_val,info,Lwork
  INTEGER, DIMENSION(:), ALLOCATABLE::iwork,ifail
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::work
 
  INTEGER::abajo,ii,jj,gg,hh
  DOUBLE PRECISION::maxL,auxval

  coordinates=0.0d0
 
  ALLOCATE(CC(Nnods,num_pivots),di(Nnods),dp(num_pivots))
  di=0.0d0
  dd=0.0d0
  dp=0.0d0
  CC=0.0d0

  CC(:,:)=dists(:,:)
 
  !! Deberia poner aqui una opcion para hacer un guardado opcional distancias=CC
  !! y usarlas luego como opcion al calcular el stress.
  print*,'calcula matriz'

  DO ii=1,num_pivots
     DO jj=1,Nnods
        CC(jj,ii)=CC(jj,ii)**2    ! Necesario en el metodo
        di(jj)=di(jj)+CC(jj,ii)
        dp(ii)=dp(ii)+CC(jj,ii)
        dd=dd+CC(jj,ii)
     END DO
  END DO
  dd=dd/(Nnods*num_pivots*1.0d0)
  di=di/(num_pivots*1.0d0)
  dp=dp/(Nnods*1.0d0)
 
  DO ii=1,num_pivots
     DO jj=1,Nnods
        CC(jj,ii)=-0.50d0*(CC(jj,ii)-di(jj)-dp(ii)+dd)
     END DO
  END DO

  DEALLOCATE(di,dp)

  print*,'construye CCC'

  ALLOCATE(CCC(num_pivots,num_pivots))

  DO ii=1,num_pivots
     DO jj=1,num_pivots
        auxval=0.0d0
        DO gg=1,Nnods
              auxval=auxval+CC(gg,jj)*CC(gg,ii)
        END DO
        CCC(jj,ii)=auxval
     END DO
  END DO

  print*,'diagonaliza'

  Lwork=8*num_pivots

  ALLOCATE (work(Lwork),iwork(5*num_pivots),ifail(num_pivots))
  ALLOCATE(eigenvects_p(num_pivots,num_pivots),eigenvals_p(num_pivots))

  eigenvals_p=0.0d0
  eigenvects_p=0.0d0
  work=0.0d0
  iwork=0
  ifail=0
  abajo=num_pivots-num_pivots+1
  
!  print*,'va a diagonalizar'

  CALL dsyevx ('V','I','U',num_pivots,CCC,num_pivots,0,0,abajo,num_pivots,0.0d0,num_val,&
       eigenvals_p,eigenvects_p,num_pivots,work,Lwork,iwork,ifail,info)
 
  IF (info/=0) THEN
     print*,"#Error with the diagonalization -MDS.f90:dsyevx-"
     print*,"#the array 'work' should has the dimension:",work(1)
  END IF
 
  DEALLOCATE (work,iwork,ifail,CCC)

  coordinates=0.0d0
 
  DO ii=1,dim
     DO jj=1,num_pivots
        DO gg=1,Nnods
           coordinates(gg,ii)=coordinates(gg,ii)+CC(gg,jj)*eigenvects_p(jj,num_pivots-ii+1)
        END DO
     END DO
  END DO

  maxL=MAXVAL(coordinates)
  maxL=maxL/50.0d0
  coordinates(:,:)=coordinates(:,:)/maxL


!  coordinates=0.0d0
! 
!  DO j=1,dim
!     coordinates(:,j)=sqrt(eigenvals(lout-j+1))*eigenvects(:,lout-j+1)
!  END DO
! 
!  IF (opt_stress==1) THEN
!     !ALLOCATE(vect_aux(Ktot),coors_aux(Nnods))
!     !vect_aux=0.0d0
!     !coors_aux=0.0d0
!     !norm=0.0d0
!     !DO ii=1,Nnods
!     !   DO jj=T_start(ii)+1,T_start(ii+1)
!     !      norm=norm+T_tau(jj)**2
!     !   END DO
!     !END DO
!     !DO j=1,lout
!     !   if (eigenvals(lout-j+1)>0.0d0) THEN
!     !      salida=0.0d0
!     !      dd=sqrt(eigenvals(lout-j+1))
!     !      coors_aux(:)=dd*eigenvects(:,lout-j+1)
!     !      DO ii=1,Nnods
!     !         DO jj=T_start(ii)+1,T_start(ii+1)
!     !            gg=T_ind(jj)
!     !            vect_aux(jj)=vect_aux(jj)+(coors_aux(ii)-coors_aux(gg))**2
!     !            salida=salida+(sqrt(vect_aux(jj))-T_tau(jj))**2
!     !         END DO
!     !      END DO
!     !      stress(j)=sqrt(salida/norm)
!     !   END IF
!     !END DO
!     !DEALLOCATE(vect_aux,coors_aux)
!  END IF
  

 
END SUBROUTINE MDS_PIVOTS


SUBROUTINE CHOOSE_RANDOM_PIVOTS_1(xnum_pivots,list_pivots)

  INTEGER,INTENT(IN)::xnum_pivots
  INTEGER,DIMENSION(xnum_pivots),INTENT(OUT)::list_pivots

  INTEGER::ii,jj
  INTEGER,DIMENSION(4)::hseed
  DOUBLE PRECISION::dice
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::acumulado
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro
  LOGICAL::salida

  num_pivots=xnum_pivots
  IF (ALLOCATED(ind_pivots)) DEALLOCATE(ind_pivots)
  ALLOCATE(ind_pivots(num_pivots))

  ALLOCATE(acumulado(Nnods),filtro(Nnods))
  acumulado=0.0d0
  acumulado(1)=Pe(1)
  DO ii=2,Nnods
     acumulado(ii)=acumulado(ii-1)+Pe(ii)
  END DO
  filtro=.FALSE.

  hseed=iseed

  ii=1
  DO WHILE (ii<=num_pivots)
     CALL dlarnv(1,hseed,1,dice)
     salida=.FALSE.
     DO jj=1,Nnods
        IF (filtro(jj).eqv..FALSE.) THEN
           IF (dice<acumulado(jj)) THEN
              salida=.TRUE.
              ind_pivots(ii)=jj
              filtro(jj)=.TRUE.
              ii=ii+1
              EXIT
           END IF
        END IF
     END DO
     IF (salida.eqv..FALSE.) THEN
        IF (filtro(Nnods).eqv..FALSE.) THEN
           ind_pivots(ii)=Nnods
           filtro(Nnods)=.TRUE.
           ii=ii+1
        END IF
     END IF
  END DO

  iseed=hseed

  print*,'listo',COUNT(filtro)

  DEALLOCATE(filtro,acumulado)

  list_pivots=ind_pivots-1

END SUBROUTINE CHOOSE_RANDOM_PIVOTS_1


SUBROUTINE CHOOSE_RANDOM_PIVOTS_2_W_DIJKSTRA(xnum_pivots,list_pivots)

  INTEGER,INTENT(IN)::xnum_pivots
  INTEGER,DIMENSION(xnum_pivots),INTENT(OUT)::list_pivots

  INTEGER::ii,jj,gg,num_candidatos
  INTEGER,DIMENSION(4)::hseed
  DOUBLE PRECISION::dice
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::acumulado
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::aux_dists
  INTEGER,DIMENSION(:),ALLOCATABLE::ind_candidatos
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro
  LOGICAL::salida

  num_pivots=xnum_pivots
  IF (ALLOCATED(ind_pivots)) DEALLOCATE(ind_pivots)
  ALLOCATE(ind_pivots(num_pivots))
  ALLOCATE(aux_dists(Nnods))

  IF (dists_up.eqv..TRUE.) THEN
     DEALLOCATE(dists)
  END IF

  ALLOCATE(dists(Nnods,num_pivots))
  dists=0.0d0
  dists_up=.TRUE.

  ALLOCATE(filtro(Nnods))

  num_candidatos=Nnods
  filtro=.TRUE.
  ALLOCATE(acumulado(num_candidatos),ind_candidatos(num_candidatos))
  aux_cut=0.0d0
  jj=0
  DO ii=1,Nnods
     IF (filtro(ii).eqv..TRUE.) THEN
        aux_cut=aux_cut+Pe(ii)
        jj=jj+1
        acumulado(jj)=aux_cut
        ind_candidatos(jj)=ii
     END IF
  END DO
  acumulado=acumulado/aux_cut

  hseed=iseed

  DO ii=1,num_pivots
     
     CALL dlarnv(1,hseed,1,dice)
     salida=.FALSE.
     DO jj=1,num_candidatos
        IF (dice<acumulado(jj)) THEN
           salida=.TRUE.
           ind_pivots(ii)=ind_candidatos(jj)
           EXIT
        END IF
     END DO
     IF (salida.eqv..FALSE.) THEN
        ind_pivots(ii)=ind_candidatos(num_candidatos)
     END IF

     CALL DIJKSTRA_A_TARGET(ind_pivots(ii),Nnods,aux_dists)
     dists(:,ii)=aux_dists(:)

     aux_dists(:)=0.0d0
     DO jj=1,Nnods
        aux_dists(jj)=MINVAL(dists(jj,1:ii))
     END DO
     aux_cut=0.30d0*MAXVAL(aux_dists(:))
     filtro=.FALSE.
     num_candidatos=0
     DO jj=1,Nnods
        IF (aux_dists(jj)>aux_cut) THEN
           filtro(jj)=.TRUE.
           num_candidatos=num_candidatos+1
        END IF
     END DO
     DEALLOCATE(acumulado,ind_candidatos)
     ALLOCATE(acumulado(num_candidatos),ind_candidatos(num_candidatos))
     aux_cut=0.0d0
     jj=0
     DO gg=1,Nnods
        IF (filtro(gg).eqv..TRUE.) THEN
           aux_cut=aux_cut+Pe(gg)
           jj=jj+1
           acumulado(jj)=aux_cut
           ind_candidatos(jj)=gg
        END IF
     END DO
     acumulado=acumulado/aux_cut

  END DO

  iseed=hseed

  DEALLOCATE(filtro,acumulado,ind_candidatos)

  list_pivots=ind_pivots-1


END SUBROUTINE CHOOSE_RANDOM_PIVOTS_2_W_DIJKSTRA


END MODULE GLOB
