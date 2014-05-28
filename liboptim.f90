MODULE GLOB

DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::dists
DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::coors

!f2py   intent(hide)::Nnods,Ktot
!f2py   intent(hide)::T_start,T_ind
!f2py   intent(hide)::T_tau
!f2py   intent(hide)::Pe
!f2py   intent(hide)::dists_up,net_up
INTEGER::Nnods,Ktot
INTEGER,DIMENSION(:),ALLOCATABLE::T_start,T_ind
DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::T_tau,Pe
LOGICAL::dists_up=.FALSE.
LOGICAL::net_up=.FALSE.
INTEGER,DIMENSION(4)::iseed

CONTAINS

SUBROUTINE INIT_ISEED()

  iseed(1)=123
  iseed(2)=456
  iseed(3)=2651
  iseed(4)=3455 !impar (limite 4094)

END SUBROUTINE INIT_ISEED

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

END SUBROUTINE LOAD_NET


SUBROUTINE LOAD_DISTS(xNnods,xdists)

  INTEGER,INTENT(IN)::xNnods
  DOUBLE PRECISION,DIMENSION(xNnods,xNnods),INTENT(IN)::xdists

  Nnods=xNnods

  IF (dists_up.eqv..TRUE.) THEN
     DEALLOCATE(dists)
  END IF

  ALLOCATE(dists(Nnods,Nnods))
  dists=xdists

  dists_up=.TRUE.

END SUBROUTINE LOAD_DISTS

SUBROUTINE RANDOM_DISTRIBUTION(tipo,dim)

  INTEGER,INTENT(IN)::tipo,dim
  DOUBLE PRECISION::lmax
  DOUBLE PRECISION,DIMENSION(dim)::dice3

  ALLOCATE(coors(dim,Nnods))

  IF (tipo==1) THEN

     lmax=MAXVAL(dists)
     lmax=(lmax)**(1.0d0/dim)
     DO ii=1,Nnods
        CALL dlarnv(1,iseed,dim,dice3)
        coors(:,ii)=lmax*dice3(:)
     END DO

  END IF

END SUBROUTINE RANDOM_DISTRIBUTION

SUBROUTINE DIFFUSION_DISTANCE (nn,xnods,distancias)

  IMPLICIT NONE
  
  INTEGER,INTENT(IN)::nn,xnods
  DOUBLE PRECISION,DIMENSION(xnods,xnods),INTENT(OUT)::distancias
  
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::transits,transits2

  INTEGER::ii,jj,kk,gg,ll,mm
  DOUBLE PRECISION::aux

  ALLOCATE(transits(xnods,xnods),transits2(xnods,xnods))

  distancias=0.0d0
  transits=0.0d0
  transits2=0.0d0

  DO ii=1,Nnods
     DO jj=T_start(ii)+1,T_start(ii+1)
        transits(T_ind(jj),ii)=T_tau(jj)
     END DO
  END DO

  distancias=transits
  DO kk=1,nn-1
     print*,kk+1
     transits2=0.0d0
     DO ii=1,Nnods
        DO jj=1,Nnods
           DO mm=T_start(jj)+1,T_start(jj+1)
              ll=T_ind(mm)
              transits2(ll,ii)=transits2(ll,ii)+T_tau(mm)*distancias(jj,ii)
           END DO
        END DO
     END DO
     distancias=transits2
     print*,SUM(distancias(:,346))
  END DO
  transits=distancias
  distancias=0.0d0

!  ALLOCATE(pesotes(Nnods,Nnods))
!  pesotes=transits


  !distancias=transits
  !DO kk=1,nn-1
  !   print*,kk+1
  !   transits2=0.0d0
  !   DO ii=1,Nnods
  !      DO jj=1,Nnods
  !         aux=0.0d0
  !         DO gg=1,Nnods
  !            aux=aux+transits(ii,gg)*distancias(gg,jj)
  !         END DO
  !         transits2(ii,jj)=aux
  !      END DO
  !   END DO
  !   distancias=transits2
  !   print*,SUM(distancias(:,346))
  !END DO
  !transits=distancias
  !distancias=0.0d0


  DO ii=1,Nnods
     DO jj=1,Nnods
        aux=0.0d0
        DO kk=1,Nnods
           aux=aux+((transits(kk,ii)-transits(kk,jj))**2)/Pe(kk)
        END DO
        distancias(jj,ii)=aux
     END DO
  END DO


END SUBROUTINE DIFFUSION_DISTANCE


SUBROUTINE MAJORIZATION (tipo_pesos,coorsin,xnods,fcoors)

  INTEGER,INTENT(IN)::tipo_pesos,xnods
  DOUBLE PRECISION,DIMENSION(3,xnods),INTENT(IN)::coorsin
  DOUBLE PRECISION,DIMENSION(3,xnods),INTENT(OUT)::fcoors

  INTEGER::ii,jj,kk
  DOUBLE PRECISION::nq,estres0,estresf,wij,denom,aa
  DOUBLE PRECISION,DIMENSION(3)::q,newq
  DOUBLE PRECISION,DIMENSION(3,xnods)::coors0
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::pesos

  ALLOCATE(pesos(xnods,xnods))

  pesos=0.0d0
  
  IF (tipo_pesos==1) THEN
     DO ii=1,xnods
        DO jj=1,xnods
           pesos(ii,jj)=1.0d0
           !pesos(ii,jj)=Pe(ii)*Pe(jj)
           !pesos(ii,jj)=dists(ii,jj)**(-2)
           !pesos(ii,jj)=PijPi
        END DO
     END DO
  ELSE IF (tipo_pesos==2) THEN
     DO ii=1,xnods
        DO jj=1,xnods
           pesos(ii,jj)=Pe(ii)*Pe(jj)
        END DO
     END DO
  ELSE IF (tipo_pesos==3) THEN
     DO ii=1,xnods
        DO jj=1,xnods
           pesos(ii,jj)=dists(ii,jj)**(-2)
        END DO
     END DO
  ELSE IF (tipo_pesos==4) THEN
     DO ii=1,xnods
        DO jj=1,xnods
           pesos(ii,jj)=Pe(ii)*Pe(jj)*dists(ii,jj)**(-2)
        END DO
     END DO
  ELSE IF (tipo_pesos==5) THEN
     DO ii=1,xnods
        DO jj=T_start(ii)+1,T_start(ii+1)
           pesos(T_ind(jj),ii)=1.0
        END DO
     END DO
  ELSE IF (tipo_pesos==6) THEN
     DO ii=1,xnods
        DO jj=T_start(ii)+1,T_start(ii+1)
           pesos(T_ind(jj),ii)=Pe(ii)*T_tau(jj)
        END DO
     END DO
  END IF


  !DO ii=1,Nnods
  !   DO jj=T_start(ii)+1,T_start(ii+1)
  !      pesos(T_ind(jj),ii)=T_tau(jj)
  !   END DO
  !END DO
  !DO ii=1,xnods
  !   DO jj=1,xnods
  !      pesos(ii,jj)=pesotes(jj,ii)
  !   END DO
  !END DO

  aa=1.0d0
  DO ii=1,xnods
     DO jj=1,xnods
        IF (ii/=jj) THEN
           aa=aa+pesos(ii,jj)*dists(ii,jj)**2
        END IF
     END DO
  END DO
  pesos(:,:)=pesos(:,:)/aa


  coors0=coorsin

  estres0=0.0d0
  DO ii=1,xnods
     DO jj=ii+1,xnods
        q(:)=coors0(:,ii)-coors0(:,jj)
        nq=sqrt(q(1)*q(1)+q(2)*q(2)+q(3)*q(3))
        nq=(nq-dists(ii,jj))**2
        estres0=estres0+pesos(ii,jj)*nq
     END DO
  END DO

  DO kk=1,2000

     DO ii=1,xnods
        denom=0.0d0
        newq=0.0d0
        DO jj=1,xnods
           IF (ii/=jj) THEN
              denom=denom+pesos(ii,jj)
              q(:)=coors0(:,ii)-coors0(:,jj)
              nq=sqrt(q(1)*q(1)+q(2)*q(2)+q(3)*q(3))
              if (nq>0.0d0) THEN
                 sij=dists(ii,jj)/nq
              else
                 sij=0.0d0
              end if
              q(:)=sij*q(:)+coors0(:,jj)
              newq(:)=newq(:)+pesos(ii,jj)*q(:)
           END IF
        END DO
        newq(:)=newq(:)/denom
        fcoors(:,ii)=newq(:)
     END DO
 
     estresf=0.0d0
     DO ii=1,xnods
        DO jj=ii+1,xnods
           q(:)=fcoors(:,ii)-fcoors(:,jj)
           nq=sqrt(q(1)*q(1)+q(2)*q(2)+q(3)*q(3))
           nq=(nq-dists(ii,jj))**2
           estresf=estresf+pesos(ii,jj)*nq
        END DO
     END DO

     print*,kk,(estres0-estresf)/estres0
     estres0=estresf
     coors0(:,:)=fcoors(:,:)

  END DO

  DEALLOCATE(pesos)

END SUBROUTINE MAJORIZATION


!SUBROUTINE CLASSICAL_SGD ()
! 
! 
!END SUBROUTINE CLASSICAL_SGD


END MODULE GLOB
