SUBROUTINE contact_map (cutoff,coors,num_atoms,map)
  
  IMPLICIT NONE

  INTEGER,INTENT(IN)::num_atoms
  REAL,INTENT(IN)::cutoff
  REAL,DIMENSION(num_atoms,3),INTENT(IN)::coors


  INTEGER::i,j
  REAL,DIMENSION(num_atoms,num_atoms),INTENT(OUT)::map
  REAL,DIMENSION(3)::qaux1,qaux2
  REAL::dd,co

  map=0.0d0
  co=cutoff**2

  DO i=1,num_atoms
     qaux1=coors(i,:)
     DO j=i+1,num_atoms
        qaux2=qaux1-coors(j,:)
        dd=dot_product(qaux2,qaux2)
        IF (dd<=co) THEN
           map(i,j)=1.0d0
           map(j,i)=1.0d0
        END IF
     END DO
  END DO

END SUBROUTINE contact_map

SUBROUTINE gnm (contact,num_nodes,values,vectors,freqs,bfact,invers,correl)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::num_nodes
  REAL,DIMENSION(num_nodes,num_nodes),INTENT(IN)::contact

  DOUBLE PRECISION,DIMENSION(num_nodes),INTENT(OUT)::values,freqs,bfact
  DOUBLE PRECISION,DIMENSION(num_nodes,num_nodes),INTENT(OUT)::vectors
  DOUBLE PRECISION,DIMENSION(num_nodes,num_nodes),INTENT(OUT)::invers,correl

  INTEGER::i,j,k
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::matrix,matrix_bck

  !Para diagonalizar:
  INTEGER::num_val,info
  INTEGER, DIMENSION(:), ALLOCATABLE::iwork,ifail
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::work,invvalues

  DOUBLE PRECISION,DIMENSION(num_nodes,num_nodes)::aux

  ALLOCATE (matrix(num_nodes,num_nodes),matrix_bck(num_nodes,num_nodes))
  ALLOCATE (invvalues(num_nodes),work(8*num_nodes),iwork(5*num_nodes),ifail(num_nodes))


  matrix=0.0d0
  matrix_bck=0.0d0
  invvalues=0.0d0
  work=0.0d0
  iwork=0
  ifail=0
  values=0.0d0
  vectors=0.0d0
  freqs=0.0d0
  bfact=0.0d0

  matrix=-contact

  !contact comes from a boolean matrix in python. F2py always translates booleans as integers (false=0, true=1)

  DO i=1,num_nodes
     DO j=i+1,num_nodes
        IF (contact(i,j)>0.0d0) THEN
           matrix(i,i)=matrix(i,i)+contact(i,j)
           matrix(j,j)=matrix(j,j)+contact(i,j)
        END IF
     END DO
  END DO

  matrix_bck=matrix

  CALL dsyevx ('V','I','U',num_nodes,matrix,num_nodes,0,0,1,num_nodes,0.0d0,num_val&
       &,values,vectors,num_nodes,work,8*num_nodes,iwork,ifail,info) 
  
  IF (info/=0) THEN
     print*,"Error with the diagonalization."
     print*,"The array 'work' should has the dimension:",work(1)
  END IF
  

  freqs=sqrt(abs(values))
  
  !Inverse matrix:

  DO i=1,num_nodes
     invvalues(i)=1.0d0/values(i)
  END DO

  invers=0.0d0

  DO k=2,num_nodes
     
     DO i=1,num_nodes
        DO j=1,num_nodes
           
           invers(i,j)=invers(i,j)+invvalues(k)*vectors(i,k)*vectors(j,k)
           
        END DO
     END DO
  END DO

  !Compruebo la inversa:
  aux=0.0d0
  DO i=1,num_nodes
     DO j=1,num_nodes
        aux(i,j)=dot_product(invers(i,:),matrix_bck(:,j))
     END DO
  END DO


  !Total correlation:
  DO i=1,num_nodes
     DO j=1,num_nodes
        correl(i,j)=invers(i,j)/(sqrt(invers(i,i)*invers(j,j)))
     END DO
  END DO

  !B Factors:
  DO i=1,num_nodes
     bfact(i)=invers(i,i)
  END DO

  !Fixing the output:
  aux=0.0d0
  aux=vectors
  vectors=0.0d0
  vectors=TRANSPOSE(aux)


END SUBROUTINE gnm

SUBROUTINE correlation (vectors,values,lista,num_nodes,dim_lista,correl)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::num_nodes,dim_lista
  DOUBLE PRECISION,DIMENSION(num_nodes*3,num_nodes*3),INTENT(IN)::vectors
  DOUBLE PRECISION,DIMENSION(num_nodes*3),INTENT(IN)::values
  INTEGER,DIMENSION(dim_lista),INTENT(IN)::lista
  INTEGER,DIMENSION(dim_lista)::flista

  DOUBLE PRECISION,DIMENSION(num_nodes,num_nodes),INTENT(OUT)::correl
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::auto_correl

  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::invers

  INTEGER::ii,jj,k,i,j,dim_aux

  dim_aux=num_nodes*3

  ALLOCATE(invers(dim_aux,dim_aux),auto_correl(num_nodes))
  correl=0.0d0
  invers=0.0d0
  auto_correl=0.0d0

  flista(:)=lista(:)+1

  DO ii=1,dim_lista
     k=flista(ii)
     DO i=1,dim_aux
        DO j=1,dim_aux
           
           invers(i,j)=invers(i,j)+(1.0d0/values(k))*vectors(k,i)*vectors(k,j)
           
        END DO
     END DO
  END DO
  
  auto_correl=0.0d0
  !autocorrelation:
  DO i=1,num_nodes
     ii=3*(i-1)
     DO j=1,3
        auto_correl(i)=auto_correl(i)+invers(ii+j,ii+j)
     END DO
  END DO

  correl=0.0d0
  !Total correlation:
  DO i=1,num_nodes
     ii=3*(i-1)
     DO j=1,num_nodes
        jj=3*(j-1)

        DO k=1,3
           correl(i,j)=correl(i,j)+invers(ii+k,jj+k)
        END DO

        correl(i,j)=correl(i,j) !/(sqrt(auto_correl(i)*auto_correl(j)))

     END DO
  END DO

END SUBROUTINE correlation

SUBROUTINE anm (contact,coors,num_nodes,values,vectors,freqs,bfact,invers,correl)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::num_nodes
  REAL,DIMENSION(num_nodes,num_nodes),INTENT(IN)::contact
  REAL,DIMENSION(num_nodes,3),INTENT(IN)::coors

  DOUBLE PRECISION,DIMENSION(num_nodes*3),INTENT(OUT)::values,freqs
  DOUBLE PRECISION,DIMENSION(num_nodes),INTENT(OUT)::bfact
  DOUBLE PRECISION,DIMENSION(num_nodes*3,num_nodes*3),INTENT(OUT)::vectors
  DOUBLE PRECISION,DIMENSION(num_nodes*3,num_nodes*3),INTENT(OUT)::invers
  DOUBLE PRECISION,DIMENSION(num_nodes,num_nodes),INTENT(OUT)::correl


  INTEGER::i,j,k,dim_aux,dim,ii,jj,g
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::matrix,matrix_bck
  DOUBLE PRECISION::val_aux,val_aux2
  DOUBLE PRECISION,DIMENSION(3)::vect_aux

  !Para diagonalizar:
  INTEGER::num_val,info
  INTEGER, DIMENSION(:), ALLOCATABLE::iwork,ifail
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::work,invvalues
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::aux

  values=0.0d0
  bfact=0.0d0
  vectors=0.0d0
  invers=0.0d0
  correl=0.0d0
  freqs=0.0d0

  dim=3
  dim_aux=dim*num_nodes

  ALLOCATE (matrix(dim_aux,dim_aux),matrix_bck(dim_aux,dim_aux))
  ALLOCATE (invvalues(dim_aux),work(8*dim_aux),iwork(5*dim_aux),ifail(dim_aux))
  ALLOCATE (aux(dim_aux,dim_aux))

  matrix=0.0d0
  matrix_bck=0.0d0
  invvalues=0.0d0
  work=0.0d0
  iwork=0
  ifail=0
  values=0.0d0
  vectors=0.0d0
  freqs=0.0d0
  bfact=0.0d0
  val_aux=0.0d0
  val_aux2=0.0d0
  vect_aux=0.0d0
  invers=0.0d0
  correl=0.0d0
  aux=0.0d0


  ! The contact matrix is a real object. (false=0.0, true=1.0 or kij)

  !! Sub matrix Hij:

  DO i=1,num_nodes
     DO j=i+1,num_nodes
        IF (contact(i,j)>0.0d0) THEN
           ii=dim*(i-1)
           jj=dim*(j-1)
           vect_aux=coors(i,:)-coors(j,:)
           val_aux=dot_product(vect_aux,vect_aux)
           DO k=1,dim
              DO g=1,dim
                 val_aux2=0.0d0
                 val_aux2=-(coors(i,k)-coors(j,k))*(coors(i,g)-coors(j,g))
                 val_aux2=(val_aux2/val_aux)*contact(i,j)   !! Here it is the for constant kij
                 matrix(ii+k,jj+g)=val_aux2
                 matrix(jj+g,ii+k)=val_aux2
              END DO
           END DO
        END IF
     END DO
  END DO

  !! Sub matrix Hii

  DO i=1,num_nodes
     ii=dim*(i-1)
     DO j=1,i-1
        jj=dim*(j-1)
        DO k=1,dim
           DO g=1,dim
              matrix(ii+k,ii+g)=matrix(ii+k,ii+g)-matrix(ii+k,jj+g)
           END DO
        END DO
     END DO

     DO j=i+1,num_nodes
        jj=dim*(j-1)
        DO k=1,dim
           DO g=1,dim
              matrix(ii+k,ii+g)=matrix(ii+k,ii+g)-matrix(ii+k,jj+g)
           END DO
        END DO
     END DO
  END DO

  matrix_bck=matrix


  CALL dsyevx ('V','I','U',dim_aux,matrix,dim_aux,0,0,1,dim_aux,0.0d0,num_val&
       &,values,vectors,dim_aux,work,8*dim_aux,iwork,ifail,info) 

  IF (info/=0) THEN
     print*,"Error with the diagonalization."
     print*,"The array 'work' should has the dimension:",work(1)
  END IF
  
  freqs=sqrt(abs(values))

  
  !Inverse matrix:

  DO i=1,dim_aux
     invvalues(i)=1.0d0/values(i)
  END DO

  invers=0.0d0

  DO k=7,dim_aux
     DO i=1,dim_aux
        DO j=1,dim_aux
           
           invers(i,j)=invers(i,j)+invvalues(k)*vectors(i,k)*vectors(j,k)
           
        END DO
     END DO
  END DO

  !Checking the inverse matrix:
  !  aux=0.0d0
  !  DO i=1,dim_aux
  !     DO j=1,dim_aux!        aux(i,j)=dot_product(invers(i,:),matrix_bck(:,j))
  !        print*,i,j,aux(i,j)
  !     END DO
  !  END DO
  !  aux=0.0d0
  !  aux=MATMUL(matrix_bck,invers)
  !  print*,aux

  bfact=0.0d0
  !B Factors:
  DO i=1,num_nodes
     ii=3*(i-1)
     DO j=1,3
        bfact(i)=bfact(i)+invers(ii+j,ii+j)
     END DO
  END DO

  correl=0.0d0
  !Total correlation:
  DO i=1,num_nodes
     ii=3*(i-1)
     DO j=1,num_nodes
        jj=3*(j-1)

        DO k=1,3
           correl(i,j)=correl(i,j)+invers(ii+k,jj+k)
        END DO

        correl(i,j)=correl(i,j)/(sqrt(bfact(i)*bfact(j)))

     END DO
  END DO


!!$  !Fixing the output:
  aux=0.0d0
  aux(:,:)=vectors(:,:)
  vectors=0.0d0
  vectors=TRANSPOSE(aux)


END SUBROUTINE anm



!!!! f2py -c -m pyn_fort_gnm pyn_fort_gnm.f90 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_def -lpthread
!!!! f2py2 -c -m pyn_fort_gnm pyn_fort_gnm.f90 -llapack

