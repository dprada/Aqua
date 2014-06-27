MODULE GLOB

DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::dists
DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::coors

!f2py   intent(hide)::Nnods,Ktot,dim
!f2py   intent(hide)::T_start,T_ind
!f2py   intent(hide)::T_tau
!f2py   intent(hide)::Pe
!f2py   intent(hide)::dists_up,net_up,coors_up
INTEGER::Nnods,Ktot,dim
INTEGER,DIMENSION(:),ALLOCATABLE::T_start,T_ind
DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::T_tau,Pe
DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::Pmat
LOGICAL::dists_up=.FALSE.
LOGICAL::net_up=.FALSE.
LOGICAL::pmat_up=.FALSE.
LOGICAL::coors_up=.FALSE.
INTEGER,DIMENSION(4)::iseed

CONTAINS

SUBROUTINE INIT_ISEED(ss1,ss2,ss3,ss4)

  INTEGER,INTENT(IN)::ss1,ss2,ss3,ss4
  iseed(1)=ss1
  iseed(2)=ss2
  iseed(3)=ss3
  iseed(4)=ss4 !impar (limite 4094)

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

SUBROUTINE RESET_COORS(xdim)

  INTEGER,INTENT(IN)::xdim

  IF (coors_up.eqv..TRUE.) THEN
     DEALLOCATE(coors)
  END IF

  dim=xdim
  ALLOCATE(coors(dim,Nnods))
  coors=0.0d0

END SUBROUTINE RESET_COORS

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

SUBROUTINE POST_LOAD_DIFF_DIST()

  INTEGER::ii,jj,kk
  DOUBLE PRECISION::aux

  IF (dists_up.eqv..TRUE.) THEN
     DEALLOCATE(dists)
  END IF

  ALLOCATE(dists(Nnods,Nnods))
  dists_up=.TRUE.


  dists=0.0d0
  DO ii=1,Nnods
     DO jj=1,Nnods
        aux=0.0d0
        DO kk=1,Nnods
           aux=aux+((pmat(kk,ii)-pmat(kk,jj))**2)/Pe(kk)
        END DO
        dists(jj,ii)=sqrt(aux)
     END DO
  END DO

END SUBROUTINE POST_LOAD_DIFF_DIST

SUBROUTINE RANDOM_DISTRIBUTION(tipo)

  INTEGER,INTENT(IN)::tipo
  DOUBLE PRECISION::lmax
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::dice3
  INTEGER::ii

  ALLOCATE(dice3(dim))

  IF (tipo==1) THEN

     lmax=MAXVAL(dists)
     lmax=lmax**(1.0d0/dim)
     DO ii=1,Nnods
        CALL dlarnv(1,iseed,dim,dice3)
        coors(:,ii)=lmax*dice3(:)
     END DO

  ELSE IF (tipo==2) THEN
     lmax=1.0d0
     DO ii=1,Nnods
        CALL dlarnv(1,iseed,dim,dice3)
        coors(:,ii)=lmax*dice3(:)
     END DO
  END IF

  DEALLOCATE(dice3)

END SUBROUTINE RANDOM_DISTRIBUTION

SUBROUTINE DIFFUSION_DISTANCE (nn)

  IMPLICIT NONE
  
  INTEGER,INTENT(IN)::nn

  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::transits

  INTEGER::ii,jj,kk,gg,ll,mm
  DOUBLE PRECISION::aux,aux2,auxpe

  ALLOCATE(transits(Nnods,Nnods))
  IF (dists_up.eqv..TRUE.) THEN
     DEALLOCATE(dists)
  END IF

  ALLOCATE(dists(Nnods,Nnods))
  dists_up=.TRUE.

  transits=0.0d0

  DO ii=1,Nnods
     DO jj=T_start(ii)+1,T_start(ii+1)
        transits(T_ind(jj),ii)=T_tau(jj)
     END DO
  END DO

  DO kk=1,nn-1
     print*,kk+1
     dists=transits
     transits=0.0d0
     DO ii=1,Nnods
        DO jj=1,Nnods
           DO mm=T_start(jj)+1,T_start(jj+1)
              ll=T_ind(mm)
              transits(ll,ii)=transits(ll,ii)+T_tau(mm)*dists(jj,ii)
           END DO
        END DO
     END DO
  END DO

  !DO ii=1,Nnods
  !   DO jj=1,Nnods
  !      transits(jj,ii)=transits(jj,ii)**2
  !   END DO
  !END DO
  ! 
  !DO kk=1,Nnods
  !   auxpe=Pe(kk)
  !   DO ii=1,Nnods
  !      aux=transits(kk,ii)/auxpe
  !      aux2=sqrt(aux)/auxpe
  !      IF (aux>0.0d0) THEN
  !         DO jj=1,Nnods
  !            IF (transits(kk,jj)>0.0d0) THEN
  !               dists(jj,ii)=dists(jj,ii)+aux+transits(kk,jj)/auxpe-aux2*sqrt(transits(kk,jj))
  !            ELSE
  !               dists(jj,ii)=dists(jj,ii)+aux
  !            END IF
  !         END DO
  !      ELSE
  !         DO jj=1,Nnods
  !            dists(jj,ii)=dists(jj,ii)+transits(kk,jj)/auxpe
  !         END DO
  !      END IF
  !   END DO
  !END DO

  dists=0.0d0
  DO ii=1,Nnods
     DO jj=1,Nnods
        aux=0.0d0
        DO kk=1,Nnods
           aux=aux+((transits(kk,ii)-transits(kk,jj))**2)/Pe(kk)
        END DO
        dists(jj,ii)=sqrt(aux)
     END DO
  END DO

  DEALLOCATE(transits)

END SUBROUTINE DIFFUSION_DISTANCE


SUBROUTINE STOCH_GRAD_DESC_MOD(tipo_pesos,iters)

  INTEGER,INTENT(IN)::tipo_pesos,iters

  INTEGER::ii,jj,kk,pin
  DOUBLE PRECISION::nq,estres0,estresf,alpha,raux
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::q,newq
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::coors0
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::pesos

  ALLOCATE(pesos(Nnods,Nnods))
  ALLOCATE(q(dim),newq(dim),coors0(dim,Nnods))

  pesos=0.0d0
  
  IF (tipo_pesos==1) THEN
     DO ii=1,Nnods
        DO jj=1,Nnods
           pesos(ii,jj)=0.10d0
        END DO
     END DO
  ELSE IF (tipo_pesos==2) THEN
     DO ii=1,Nnods
        DO jj=1,Nnods
           pesos(ii,jj)=Pe(ii)*Pe(jj)
        END DO
     END DO
  ELSE IF (tipo_pesos==3) THEN
     DO ii=1,Nnods
        DO jj=1,Nnods
           pesos(ii,jj)=dists(ii,jj)**(-2)
        END DO
     END DO
  ELSE IF (tipo_pesos==4) THEN
     DO ii=1,Nnods
        DO jj=1,Nnods
          pesos(ii,jj)=Pe(ii)*Pe(jj)*dists(ii,jj)**(-2)
        END DO
     END DO
  ELSE IF (tipo_pesos==5) THEN
     DO ii=1,Nnods
        DO jj=T_start(ii)+1,T_start(ii+1)
           pesos(T_ind(jj),ii)=1.0
        END DO
     END DO
  ELSE IF (tipo_pesos==6) THEN
     DO ii=1,Nnods
        DO jj=T_start(ii)+1,T_start(ii+1)
           pesos(T_ind(jj),ii)=Pe(ii)*T_tau(jj)
        END DO
     END DO
  END IF


  coors0=coors

  estres0=0.0d0
  DO ii=1,Nnods
     DO jj=ii+1,Nnods
        q(:)=coors0(:,ii)-coors0(:,jj)
        nq=sqrt(q(1)*q(1)+q(2)*q(2)+q(3)*q(3))
        nq=(nq-dists(ii,jj))**2
        estres0=estres0+pesos(ii,jj)*nq
     END DO
  END DO

  print*, estres0

  DO kk=1,iters

     ! #Elijo nodo para hacer pin at random
     CALL dlarnv(1,iseed,1,raux)
     pin=INT(Nnods*raux)+1
     !alpha=(2.0d0/(1.0d0+kk))
     alpha=0.5
     IF (kk>1000) THEN
        alpha=0.50d0/(1.0d0+kk-1000.0d0)
     END IF

     DO ii=1,Nnods
        IF (ii/=pin) THEN
           q(:)=coors0(:,ii)-coors0(:,pin)
           nq=sqrt(q(1)*q(1)+q(2)*q(2)+q(3)*q(3))
           q(:)=alpha*(dists(ii,pin)-nq)*(q(:)/nq)
           coors0(:,ii)=coors0(:,ii)+q(:)
        END IF
     END DO

     estresf=0.0d0
     DO ii=1,Nnods
        DO jj=ii+1,Nnods
           q(:)=coors0(:,ii)-coors0(:,jj)
           nq=sqrt(q(1)*q(1)+q(2)*q(2)+q(3)*q(3))
           nq=(nq-dists(ii,jj))**2
           estresf=estresf+pesos(ii,jj)*nq
        END DO
     END DO

     print*,estresf
     print*,kk,(estres0-estresf)/estres0
     estres0=estresf

  END DO
        
  coors=coors0
  DEALLOCATE(pesos)
  DEALLOCATE(q,newq,coors0)


END SUBROUTINE STOCH_GRAD_DESC_MOD


SUBROUTINE STOCH_GRAD_DESC_TSNE(iters1, iters2)
 
  INTEGER,INTENT(IN)::iters1,iters2
 
  INTEGER::kk,ii,jj,max_iter,pin
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::PP,QQ
  DOUBLE PRECISION::val_aux_tot,val_aux_i,raux,nq,alpha,estresf,estres0
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::q

  print*,'entra'
  ALLOCATE(PP(Nnods,Nnods),QQ(Nnods,Nnods),q(dim))
  PP=0.0d0
  QQ=0.0d0
  q=0.0d0


  val_aux_tot=0.0d0
  DO ii=1,Nnods
     DO jj=1,Nnods
        val_aux_i=Pmat(jj,ii)+Pmat(ii,jj)
        PP(jj,ii)=val_aux_i
        val_aux_tot=val_aux_tot+val_aux_i
     END DO
  END DO

  PP=PP/val_aux_tot
  PP=PP*4.0d0
  !P = Math.maximum(P, 1e-12);
  print*,'sale2'


  !compute pairwise affinities
  val_aux_tot=0.0d0
  DO ii=1,Nnods
     DO jj=ii+1,Nnods
        q(:)=coors(:,ii)-coors(:,jj)
        nq=q(1)*q(1)+q(2)*q(2)+q(3)*q(3)
        val_aux_i=1.0d0/(1.0d0+nq)
        QQ(ii,jj)=val_aux_i
        QQ(jj,ii)=val_aux_i
        val_aux_tot=val_aux_tot+2.0d0*val_aux_i
     END DO
     QQ(ii,ii)=0.0d0
  END DO
  QQ=QQ/val_aux_tot

  estres0=0.0d0
  DO ii=1,Nnods
     DO jj=1,Nnods
        IF ((PP(ii,jj)>0.0d0).and.(ii/=jj)) THEN
           val_aux_i=PP(ii,jj)*log(PP(ii,jj)/QQ(ii,jj))
           IF (ISNAN(val_aux_i).eqv..TRUE.) val_aux_i=0.0d0
           estres0=estres0+val_aux_i
        END IF
     END DO
  END DO

  print*,'0', estres0

  max_iter=iters2
 
  DO kk=1,max_iter
     CALL dlarnv(1,iseed,1,raux)
     pin=INT(Nnods*raux)+1
     alpha=500.0
     IF (kk>500) THEN
        alpha=500.0d0/(1.0d0+kk-500.0d0)
     END IF
 
     !DO ii=1,Nnods
     ! 
     !   q(:)=0.0d0
     !   DO jj=1,Nnods
     !      val_aux_i=(PP(ii,jj)-QQ(ii,jj))*QQ(ii,jj)*val_aux_tot
     !      q(:)=q(:)+val_aux_i*(coors(:,ii)-coors(:,jj))
     !   END DO
     ! 
     !   delta_q(:,ii)=q
     ! 
     !END DO


     DO ii=1,Nnods
        IF (ii/=pin) THEN
           val_aux_i=(PP(ii,pin)-QQ(ii,pin))*QQ(ii,pin)*val_aux_tot
           q(:)=alpha*val_aux_i*q(:)
           coors(:,ii)=coors(:,ii)-alpha*val_aux_i*(coors(:,ii)-coors(:,pin))
        END IF
     END DO

     !compute pairwise affinities
     val_aux_tot=0.0d0
     DO ii=1,Nnods
        DO jj=ii+1,Nnods
           q(:)=coors(:,ii)-coors(:,jj)
           nq=q(1)*q(1)+q(2)*q(2)+q(3)*q(3)
           val_aux_i=1.0d0/(1.0d0+nq)
           QQ(ii,jj)=val_aux_i
           QQ(jj,ii)=val_aux_i
           val_aux_tot=val_aux_tot+2.0d0*val_aux_i
        END DO
        QQ(ii,ii)=0.0d0
     END DO
     QQ=QQ/val_aux_tot
     ! Q = Math.maximum(Q, 1e-12)

     estresf=0.0d0
     DO ii=1,Nnods
        DO jj=1,Nnods
           IF ((PP(ii,jj)>0.0d0).and.(ii/=jj)) THEN
              val_aux_i=PP(ii,jj)*log(PP(ii,jj)/QQ(ii,jj))
              IF (ISNAN(val_aux_i).eqv..TRUE.) val_aux_i=0.0d0
              estresf=estresf+val_aux_i
           END IF
        END DO
     END DO

     print*,kk,estresf

     estres0=estresf
 
  END DO
 
 
END SUBROUTINE STOCH_GRAD_DESC_TSNE



SUBROUTINE STOCH_GRAD_DESC_TSNE_CLAS(iters1, iters2)
 
  INTEGER,INTENT(IN)::iters1,iters2

  INTEGER::kk,ii,jj,max_iter,aa,bb
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::PP,QQ,delta_q,i_q,gains
  DOUBLE PRECISION::val_aux_tot,val_aux_i,raux,min_gain,eta,estres0,estresf
  DOUBLE PRECISION::momentum,initial_momentum,final_momentum,nq
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::q

  min_gain = 0.010d0
  eta = 500.0d0
  initial_momentum=0.50d0
  final_momentum=0.80d0

  !# X2P
  print*,'entra'
  ALLOCATE(PP(Nnods,Nnods),QQ(Nnods,Nnods),q(dim))
  ALLOCATE(delta_q(dim,Nnods),i_q(dim,Nnods),gains(dim,Nnods))
  PP=0.0d0
  QQ=0.0d0
  q=0.0d0
  delta_q=0.0d0
  i_q=0.0d0
  gains=1.0d0

  print*,'sale'

  val_aux_tot=0.0d0
  DO ii=1,Nnods
     DO jj=1,Nnods
        val_aux_i=Pmat(jj,ii)+Pmat(ii,jj)
        PP(jj,ii)=val_aux_i
        val_aux_tot=val_aux_tot+val_aux_i
     END DO
  END DO

  PP=PP/val_aux_tot
  PP=PP*4.0d0
  !P = Math.maximum(P, 1e-12);
  print*,'sale2'


  !compute pairwise affinities
  val_aux_tot=0.0d0
  DO ii=1,Nnods
     DO jj=ii+1,Nnods
        q(:)=coors(:,ii)-coors(:,jj)
        nq=q(1)*q(1)+q(2)*q(2)+q(3)*q(3)
        val_aux_i=1.0d0/(1.0d0+nq)
        QQ(ii,jj)=val_aux_i
        QQ(jj,ii)=val_aux_i
        val_aux_tot=val_aux_tot+2.0d0*val_aux_i
     END DO
     QQ(ii,ii)=0.0d0
  END DO
  QQ=QQ/val_aux_tot

  estres0=0.0d0
  DO ii=1,Nnods
     DO jj=1,Nnods
        IF ((PP(ii,jj)>0.0d0).and.(ii/=jj)) THEN
           val_aux_i=PP(ii,jj)*log(PP(ii,jj)/QQ(ii,jj))
           IF (ISNAN(val_aux_i).eqv..TRUE.) val_aux_i=0.0d0
           estres0=estres0+val_aux_i
        END IF
     END DO
  END DO

  print*,'0', estres0

  max_iter=iters2

  DO kk=1,max_iter

     DO ii=1,Nnods

        q(:)=0.0d0
        DO jj=1,Nnods
           val_aux_i=(PP(ii,jj)-QQ(ii,jj))*QQ(ii,jj)*val_aux_tot
           q(:)=q(:)+val_aux_i*(coors(:,ii)-coors(:,jj))
        END DO

        delta_q(:,ii)=q

     END DO

     IF (kk<20) THEN

        momentum=initial_momentum

     ELSE

        momentum=final_momentum

     END IF

     DO ii=1,Nnods
        DO jj=1,3
           aa=0
           bb=0
           IF (delta_q(jj,ii)>0) aa=1
           IF (i_q(jj,ii)>0) bb=1
           IF (aa/=bb) THEN
              gains(jj,ii)=gains(jj,ii)+0.20d0
           ELSE
              gains(jj,ii)=gains(jj,ii)*0.80d0
           END IF
           IF (gains(jj,ii)<min_gain) gains(jj,ii)=min_gain
           i_q(jj,ii)=momentum*i_q(jj,ii)-eta*gains(jj,ii)*delta_q(jj,ii)
        END DO
     END DO

     coors=coors+i_q

     q=0.0d0
     DO ii=1,Nnods
        q=q+coors(:,ii)
     END DO
     q=q/(Nnods*1.0d0)
     DO ii=1,Nnods
        coors(:,ii)=coors(:,ii)-q
     END DO

     !compute pairwise affinities
     val_aux_tot=0.0d0
     DO ii=1,Nnods
        DO jj=ii+1,Nnods
           q(:)=coors(:,ii)-coors(:,jj)
           nq=q(1)*q(1)+q(2)*q(2)+q(3)*q(3)
           val_aux_i=1.0d0/(1.0d0+nq)
           QQ(ii,jj)=val_aux_i
           QQ(jj,ii)=val_aux_i
           val_aux_tot=val_aux_tot+2.0d0*val_aux_i
        END DO
        QQ(ii,ii)=0.0d0
     END DO
     QQ=QQ/val_aux_tot
     ! Q = Math.maximum(Q, 1e-12)

     estresf=0.0d0
     DO ii=1,Nnods
        DO jj=1,Nnods
           IF ((PP(ii,jj)>0.0d0).and.(ii/=jj)) THEN
              val_aux_i=PP(ii,jj)*log(PP(ii,jj)/QQ(ii,jj))
              IF (ISNAN(val_aux_i).eqv..TRUE.) val_aux_i=0.0d0
              estresf=estresf+val_aux_i
           END IF
        END DO
     END DO

     print*,kk,estresf
     !print*,kk,(estres0-estresf)/estres0
     !estres0=estresf

     IF (kk==iters1) THEN
        PP=PP*0.25
     END IF

  END DO

 
END SUBROUTINE STOCH_GRAD_DESC_TSNE_CLAS

SUBROUTINE STOCH_GRAD_DESC_TSNE_CLAS2(iters1, iters2,alpha)
 
  INTEGER,INTENT(IN)::iters1,iters2
  DOUBLE PRECISION,INTENT(IN)::alpha

  INTEGER::kk,ii,jj,max_iter,aa,bb,pin,bandera
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::PP,QQ,delta_q,i_q,gains
  DOUBLE PRECISION::val_aux_tot,val_aux_i,raux,min_gain,eta
  DOUBLE PRECISION::momentum,initial_momentum,final_momentum,nq
  DOUBLE PRECISION::estresf_tsne,estresf_mds,estresf,estres0
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::q,q1,q2,qt

  bandera=50
  
  min_gain = 0.010d0
  eta = 500.0d0
  initial_momentum=0.50d0
  final_momentum=0.80d0

  !# X2P
  print*,'entra'
  ALLOCATE(PP(Nnods,Nnods),QQ(Nnods,Nnods),q(dim),q1(dim),q2(dim),qt(dim))
  ALLOCATE(delta_q(dim,Nnods),i_q(dim,Nnods),gains(dim,Nnods))
  PP=0.0d0
  QQ=0.0d0
  q=0.0d0
  q2=0.0d0
  delta_q=0.0d0
  i_q=0.0d0
  gains=1.0d0

  print*,'sale'

  val_aux_tot=0.0d0
  DO ii=1,Nnods
     DO jj=1,Nnods
        val_aux_i=Pmat(jj,ii)+Pmat(ii,jj)
        PP(jj,ii)=val_aux_i
        val_aux_tot=val_aux_tot+val_aux_i
     END DO
  END DO

  PP=PP/val_aux_tot
  PP=PP*4.0d0
  !P = Math.maximum(P, 1e-12);
  print*,'sale2'


  !compute pairwise affinities
  val_aux_tot=0.0d0
  DO ii=1,Nnods
     DO jj=ii+1,Nnods
        q(:)=coors(:,ii)-coors(:,jj)
        nq=q(1)*q(1)+q(2)*q(2)+q(3)*q(3)
        val_aux_i=1.0d0/(1.0d0+nq)
        QQ(ii,jj)=val_aux_i
        QQ(jj,ii)=val_aux_i
        val_aux_tot=val_aux_tot+2.0d0*val_aux_i
     END DO
     QQ(ii,ii)=0.0d0
  END DO
  QQ=QQ/val_aux_tot

  estres0=0.0d0
  estresf_tsne=0.0d0
  estresf_mds=0.0d0
  print*,Nnods
  print*,alpha

  DO ii=1,Nnods
     DO jj=1,Nnods
        IF ((PP(ii,jj)>0.0d0).and.(ii/=jj)) THEN
           val_aux_i=PP(ii,jj)*log(PP(ii,jj)/QQ(ii,jj))
           IF (ISNAN(val_aux_i).eqv..TRUE.) val_aux_i=0.0d0
           estresf_tsne=estresf_tsne+val_aux_i
        END IF
     END DO
  END DO

  DO ii=1,Nnods
     DO jj=ii+1,Nnods
        q2(:)=coors(:,ii)-coors(:,jj)
        nq=sqrt(q2(1)*q2(1)+q2(2)*q2(2)+q2(3)*q2(3))
        nq=(nq-dists(ii,jj))**2
        estresf_mds=estresf_mds+nq
     END DO
  END DO
  estresf_mds=estresf_mds*0.001d0

  estres0=(1.0d0-alpha)*estresf_tsne+alpha*estresf_mds

  print*,'0', estres0

  max_iter=iters2

  DO kk=1,max_iter

     DO ii=1,Nnods

        q=0.0d0
        q1=0.0d0
        q2=0.0d0

        !tsne
        DO jj=1,Nnods

           qt=(coors(:,ii)-coors(:,jj))

           val_aux_i=(PP(ii,jj)-QQ(ii,jj))*QQ(ii,jj)*val_aux_tot
           q1(:)=q1(:)+val_aux_i*qt(:)

           IF (ii/=jj) THEN
              nq=sqrt(qt(1)*qt(1)+qt(2)*qt(2)+qt(3)*qt(3))
              q2(:)=q2(:)+(dists(ii,jj)-nq)*(qt(:)/nq)
           END IF

        END DO

        nq=sqrt(q2(1)*q2(1)+q2(2)*q2(2)+q2(3)*q2(3))
        q2=0.00010*(q2/nq)

        q=(1.0d0-alpha)*q1-alpha*q2
        delta_q(:,ii)=q

     END DO

     IF (kk<20) THEN

        momentum=initial_momentum

     ELSE

        momentum=final_momentum

     END IF

     DO ii=1,Nnods
        DO jj=1,3
           aa=0
           bb=0
           IF (delta_q(jj,ii)>0) aa=1
           IF (i_q(jj,ii)>0) bb=1
           IF (aa/=bb) THEN
              gains(jj,ii)=gains(jj,ii)+0.20d0
           ELSE
              gains(jj,ii)=gains(jj,ii)*0.80d0
           END IF
           IF (gains(jj,ii)<min_gain) gains(jj,ii)=min_gain
           i_q(jj,ii)=momentum*i_q(jj,ii)-eta*gains(jj,ii)*delta_q(jj,ii)
        END DO
     END DO

     coors=coors+i_q

     q=0.0d0
     DO ii=1,Nnods
        q=q+coors(:,ii)
     END DO
     q=q/(Nnods*1.0d0)
     DO ii=1,Nnods
        coors(:,ii)=coors(:,ii)-q
     END DO

     !compute pairwise affinities
     val_aux_tot=0.0d0
     DO ii=1,Nnods
        DO jj=ii+1,Nnods
           q(:)=coors(:,ii)-coors(:,jj)
           nq=q(1)*q(1)+q(2)*q(2)+q(3)*q(3)
           val_aux_i=1.0d0/(1.0d0+nq)
           QQ(ii,jj)=val_aux_i
           QQ(jj,ii)=val_aux_i
           val_aux_tot=val_aux_tot+2.0d0*val_aux_i
        END DO
        QQ(ii,ii)=0.0d0
     END DO
     QQ=QQ/val_aux_tot
     ! Q = Math.maximum(Q, 1e-12)



     IF (kk==bandera) THEN

        estresf=0.0d0
        estresf_tsne=0.0d0
        estresf_mds=0.0d0

        DO ii=1,Nnods
           DO jj=ii+1,Nnods
              IF ((PP(ii,jj)>0.0d0).and.(ii/=jj)) THEN
                 val_aux_i=PP(ii,jj)*log(PP(ii,jj)/QQ(ii,jj))
                 IF (ISNAN(val_aux_i).eqv..TRUE.) val_aux_i=0.0d0
                 estresf_tsne=estresf_tsne+val_aux_i
              END IF
              q2(:)=coors(:,ii)-coors(:,jj)
              nq=sqrt(q2(1)*q2(1)+q2(2)*q2(2)+q2(3)*q2(3))
              nq=(nq-dists(ii,jj))**2
              estresf_mds=estresf_mds+nq
           END DO
        END DO

        estresf_mds=estresf_mds*0.001d0
        estresf=(1.0d0-alpha)*estresf_tsne+alpha*estresf_mds

        print*,kk,estresf
        !print*,kk,(estres0-estresf)/estres0
        !estres0=estresf
        bandera=bandera+50

     END IF

     IF (kk==iters1) THEN
        PP=PP*0.25
     END IF

  END DO

 
END SUBROUTINE STOCH_GRAD_DESC_TSNE_CLAS2



SUBROUTINE STOCH_GRAD_DESC_CCA(tipo_pesos,iters)

  INTEGER,INTENT(IN)::tipo_pesos,iters

  INTEGER::ii,jj,kk,pin
  DOUBLE PRECISION::nq,estres0,estresf,alpha,raux,EFE,sigma,sigma_init,sigma_fin,delta_sigma
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::q,newq
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::coors0
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::pesos

  ALLOCATE(pesos(Nnods,Nnods))
  ALLOCATE(q(dim),newq(dim),coors0(dim,Nnods))

  pesos=0.0d0
  EFE=0.0d0

  IF (tipo_pesos==1) THEN
     DO ii=1,Nnods
        DO jj=1,Nnods
           pesos(ii,jj)=0.10d0
        END DO
     END DO
  ELSE IF (tipo_pesos==2) THEN
     DO ii=1,Nnods
        DO jj=1,Nnods
           pesos(ii,jj)=Pe(ii)*Pe(jj)
        END DO
     END DO
  ELSE IF (tipo_pesos==3) THEN
     DO ii=1,Nnods
        DO jj=1,Nnods
           pesos(ii,jj)=dists(ii,jj)**(-2)
        END DO
     END DO
  ELSE IF (tipo_pesos==4) THEN
     DO ii=1,Nnods
        DO jj=1,Nnods
          pesos(ii,jj)=Pe(ii)*Pe(jj)*dists(ii,jj)**(-2)
        END DO
     END DO
  ELSE IF (tipo_pesos==5) THEN
     DO ii=1,Nnods
        DO jj=T_start(ii)+1,T_start(ii+1)
           pesos(T_ind(jj),ii)=1.0
        END DO
     END DO
  ELSE IF (tipo_pesos==6) THEN
     DO ii=1,Nnods
        DO jj=T_start(ii)+1,T_start(ii+1)
           pesos(T_ind(jj),ii)=Pe(ii)*T_tau(jj)
        END DO
     END DO
  END IF


  coors0=coors

  estres0=0.0d0
  DO ii=1,Nnods
     DO jj=ii+1,Nnods
        q(:)=coors0(:,ii)-coors0(:,jj)
        nq=sqrt(q(1)*q(1)+q(2)*q(2)+q(3)*q(3))
        nq=(nq-dists(ii,jj))**2
        estres0=estres0+pesos(ii,jj)*nq
     END DO
  END DO

  print*, estres0

  sigma_init=MAXVAL(dists)
  sigma_fin=0.0d0
  delta_sigma=(sigma_init-sigma_fin)/(iters-500.0d0)
  sigma=sigma_init


  DO kk=1,iters

     ! #Elijo nodo para hacer pin at random
     CALL dlarnv(1,iseed,1,raux)
     pin=INT(Nnods*raux)+1
     !alpha=(2.0d0/(1.0d0+kk))
     alpha=1.0
     IF (kk>500) THEN
        alpha=1.0d0/(1.0d0+kk-500.0d0)
        !sigma=sigma-delta_sigma
     END IF



     DO ii=1,Nnods
        IF (ii/=pin) THEN
           q(:)=coors0(:,ii)-coors0(:,pin)
           nq=sqrt(q(1)*q(1)+q(2)*q(2)+q(3)*q(3))
           q(:)=alpha*(dists(ii,pin)-nq)*(q(:)/nq)
           !IF (nq<=sigma) THEN
           !   EFE=1.0d0
           !ELSE
           !   EFE=0.0d0
           !END IF
           coors0(:,ii)=coors0(:,ii)+q(:)
        END IF
     END DO

     estresf=0.0d0
     DO ii=1,Nnods
        DO jj=ii+1,Nnods
           q(:)=coors0(:,ii)-coors0(:,jj)
           nq=sqrt(q(1)*q(1)+q(2)*q(2)+q(3)*q(3))
           nq=(nq-dists(ii,jj))**2
           estresf=estresf+pesos(ii,jj)*nq
        END DO
     END DO

     print*,estresf
     print*,kk,(estres0-estresf)/estres0
     estres0=estresf

  END DO
        
  coors=coors0
  DEALLOCATE(pesos)
  DEALLOCATE(q,newq,coors0)


END SUBROUTINE STOCH_GRAD_DESC_CCA



SUBROUTINE MAJORIZATION (tipo_pesos,iters)

  INTEGER,INTENT(IN)::tipo_pesos,iters

  INTEGER::ii,jj,kk
  DOUBLE PRECISION::nq,estres0,estresf,wij,denom,aa

  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::q,newq
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::coors0
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::pesos

  ALLOCATE(pesos(Nnods,Nnods))
  ALLOCATE(q(dim),newq(dim),coors0(dim,Nnods))

  pesos=0.0d0
  
  IF (tipo_pesos==1) THEN
     DO ii=1,Nnods
        DO jj=1,Nnods
           pesos(ii,jj)=1.0d0
        END DO
     END DO
  ELSE IF (tipo_pesos==2) THEN
     DO ii=1,Nnods
        DO jj=1,Nnods
           pesos(ii,jj)=Pe(ii)*Pe(jj)
        END DO
     END DO
  ELSE IF (tipo_pesos==3) THEN
     DO ii=1,Nnods
        DO jj=1,Nnods
           pesos(ii,jj)=dists(ii,jj)**(-2)
        END DO
     END DO
  ELSE IF (tipo_pesos==4) THEN
     DO ii=1,Nnods
        DO jj=1,Nnods
           pesos(ii,jj)=Pe(ii)*Pe(jj)*dists(ii,jj)**(-2)
        END DO
     END DO
  ELSE IF (tipo_pesos==5) THEN
     DO ii=1,Nnods
        DO jj=T_start(ii)+1,T_start(ii+1)
           pesos(T_ind(jj),ii)=1.0
        END DO
     END DO
  ELSE IF (tipo_pesos==6) THEN
     DO ii=1,Nnods
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
  !DO ii=1,Nnods
  !   DO jj=1,Nnods
  !      pesos(ii,jj)=pesotes(jj,ii)
  !   END DO
  !END DO

  aa=1.0d0
  DO ii=1,Nnods
     DO jj=1,Nnods
        IF (ii/=jj) THEN
           aa=aa+pesos(ii,jj)*dists(ii,jj)**2
        END IF
     END DO
  END DO
  pesos(:,:)=pesos(:,:)/aa


  coors0=coors

  estres0=0.0d0
  DO ii=1,Nnods
     DO jj=ii+1,Nnods
        q(:)=coors0(:,ii)-coors0(:,jj)
        nq=sqrt(q(1)*q(1)+q(2)*q(2)+q(3)*q(3))
        nq=(nq-dists(ii,jj))**2
        estres0=estres0+pesos(ii,jj)*nq
     END DO
  END DO

  print*, estres0

  DO kk=1,iters

     DO ii=1,Nnods
        denom=0.0d0
        newq=0.0d0
        DO jj=1,Nnods
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
        coors(:,ii)=newq(:)
     END DO
 
     estresf=0.0d0
     DO ii=1,Nnods
        DO jj=ii+1,Nnods
           q(:)=coors(:,ii)-coors(:,jj)
           nq=sqrt(q(1)*q(1)+q(2)*q(2)+q(3)*q(3))
           nq=(nq-dists(ii,jj))**2
           estresf=estresf+pesos(ii,jj)*nq
        END DO
     END DO

     print*,kk,estresf,(estres0-estresf)/estres0
     estres0=estresf
     coors0(:,:)=coors(:,:)

  END DO

  DEALLOCATE(pesos)
  DEALLOCATE(q,newq,coors0)

END SUBROUTINE MAJORIZATION




!SUBROUTINE CLASSICAL_SGD ()
! 
! 
!END SUBROUTINE CLASSICAL_SGD

SUBROUTINE LOAD_PES_250()

  OPEN(UNIT=250,NAME='fort.250',STATUS='OLD',ACTION='READ')

  IF (pmat_up.eqv..FALSE.) THEN
     pmat_up=.TRUE.
     ALLOCATE(pmat(Nnods,Nnods))
  ELSE
     DEALLOCATE(pmat)
     ALLOCATE(pmat(Nnods,Nnods))
  END IF

  DO ii=1,Nnods
     DO jj=1,Nnods
        READ(250,*) pmat(jj,ii)
     END DO
  END DO

  CLOSE(250)

END SUBROUTINE LOAD_PES_250

SUBROUTINE EVOLVE_PES(iters)

  INTEGER,INTENT(IN)::iters

  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::mataux

  IF (pmat_up.eqv..FALSE.) THEN
     pmat_up=.TRUE.
     ALLOCATE(pmat(Nnods,Nnods))
  ELSE
     DEALLOCATE(pmat)
     ALLOCATE(pmat(Nnods,Nnods))
  END IF

  pmat=0.0d0

  DO ii=1,Nnods
     DO jj=T_start(ii)+1,T_start(ii+1)
        pmat(T_ind(jj),ii)=T_tau(jj)
     END DO
  END DO

  IF (iters>1) THEN
     ALLOCATE(mataux(Nnods,Nnods))
  END IF

  DO kk=1,iters-1
     print*,kk+1
     mataux=pmat
     pmat=0.0d0
     DO ii=1,Nnods
        DO jj=1,Nnods
           DO mm=T_start(jj)+1,T_start(jj+1)
              ll=T_ind(mm)
              pmat(ll,ii)=pmat(ll,ii)+T_tau(mm)*mataux(jj,ii)
           END DO
        END DO
     END DO
  END DO

  IF (iters>1) THEN
     DEALLOCATE(mataux)
  END IF


  !DO ii=1,Nnods
  !   DO jj=1,Nnods
  !      WRITE(250,*) pmat(jj,ii)
  !   END DO
  !END DO


END SUBROUTINE EVOLVE_PES

!SUBROUTINE TSNE (iters)
! 
!  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::qji
!  DOUBLE PRECISION::Q
! 
!  initial_momentum=0.5
!  final_momentum=0.8
! 
!  DO ii=1,Nnods
!     DO jj=ii+1,Nnods
!        valaux=(pmat(jj,ii)+pmat(ii,jj))/(2*Nnods)
!        pmat(jj,ii)=valaux
!        pmat(ii,jj)=valaux
!     END DO
!  END DO
! 
!  DO itt=1,iters
! 
!     !Computo las afinidades por parejas
! 
!     Q=0.0d0
!     DO ii=1,Nnods
!        DO jj=1,Nnods
!           IF (ii==jj) THEN
!              qji(jj,ii)=0.0d0
!           ELSE
!              q(:)=coors(:,ii)-coors(:,jj)
!              nq=sqrt(q(1)*q(1)+q(2)*q(2)+q(3)*q(3))
!              valaux=1.0d0/(1.0d0+nq)
!              qji(jj,ii)=valaux
!              Q=Q+valaux
!           END IF
!        END DO
!     END DO
!     qji=qji/Q
!     
!     !Computo los gradientes
! 
!     grady=0.0d0
!     DO ii=1,Nnods
!        DO jj=1,Nnods
!           IF (ii/=jj) THEN
!              valaux=(pmat(jj,ii)-qmat(jj,ii))*qmat(jj,ii)*Q
!              grady(:,ii)=grady(:,ii)+valaux*(coors(:,ii)-coors(:,jj))
!           END IF
!        END DO
!     END DO
! 
!     IF (itt<250) THEN
!        momentum=initial_momentum
!     ELSE
!        momentum=final_momentum
!     END IF
! 
!     
! 
! 
!     
! 
!END SUBROUTINE TSNE


END MODULE GLOB
