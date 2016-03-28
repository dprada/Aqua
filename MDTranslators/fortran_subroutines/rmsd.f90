MODULE para_rmsd

IMPLICIT NONE

CONTAINS

subroutine min_rmsd(N,struct_ref,struct_2,U,center_ref,center_2,rmsd,g)

INTEGER,INTENT(IN)::N
DOUBLE PRECISION,DIMENSION(N,3),INTENT(IN)::struct_ref,struct_2
DOUBLE PRECISION,DIMENSION(3,3),INTENT(OUT)::U
DOUBLE PRECISION,DIMENSION(3),INTENT(OUT)::center_ref,center_2
DOUBLE PRECISION,INTENT(OUT)::rmsd
DOUBLE PRECISION,DIMENSION(N,3),INTENT(OUT)::g

INTEGER::i,j,k
DOUBLE PRECISION,DIMENSION(N,3)::x,y
DOUBLE PRECISION,DIMENSION(N)::w
DOUBLE PRECISION::sw,msd,x_norm,y_norm
DOUBLE PRECISION,DIMENSION(3,3)::R
DOUBLE PRECISION,DIMENSION(4,4)::F
DOUBLE PRECISION,DIMENSION(3)::tmp

!para diagonalizar:
DOUBLE PRECISION,DIMENSION(4,4)::CC
INTEGER::num_val,info
INTEGER, DIMENSION(:),ALLOCATABLE::iwork,ifail
DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE::values,work
DOUBLE PRECISION, DIMENSION(:,:),ALLOCATABLE::vectors


ALLOCATE(values(4),vectors(4,4),work(8*4),iwork(5*4),ifail(4))

w=0.0d0
x=0.0d0
y=0.0d0
CC=0.0d0
rmsd=0.0d0
msd=0.0d0
sw=0.0d0
values=0.0d0
vectors=0.0d0
center_ref=0.0d0
center_2=0.0d0
U=0.0d0
g=0.0d0
x_norm=0.0d0
y_norm=0.0d0
R=0.0d0
F=0.0d0
tmp=0.0d0


!!! copio y peso las coordenadas:
w=1.0d0
DO i=1,N
   sw=w(i)
   x(i,:)=sw*struct_ref(i,:)
   y(i,:)=sw*struct_2(i,:)
END DO

!!! calculo baricentros, centroides y normas:

DO i=1,3
   center_ref(i)=sum(x(:,i))/dble(N)
   center_2(i)=sum(y(:,i))/dble(N)
   x(:,i)=x(:,i)-center_ref(i)
   y(:,i)=y(:,i)-center_2(i)
   x_norm=x_norm+dot_product(x(:,i),x(:,i))
   y_norm=y_norm+dot_product(y(:,i),y(:,i))
END DO

!!! calculo la matriz R
DO i=1,3
   DO j=1,3
      R(i,j)=dot_product(x(:,i),y(:,j))
   END DO
END DO

!!! construimos la matriz F:

F(1,1)=R(1,1)+R(2,2)+R(3,3)
F(2,1)=R(2,3)-R(3,2)
F(3,1)=R(3,1)-R(1,3)
F(4,1)=R(1,2)-R(2,1)
F(1,2)=F(2,1)
F(2,2)=R(1,1)-R(2,2)-R(3,3)
F(3,2)=R(1,2)+R(2,1)
F(4,2)=R(1,3)+R(3,1)
F(1,3)=F(3,1)
F(2,3)=F(3,2)
F(3,3)=-R(1,1)+R(2,2)-R(3,3)
F(4,3)=R(2,3)+R(3,2)
F(1,4)=F(4,1)
F(2,4)=F(4,2)
F(3,4)=F(4,3)
F(4,4)=-R(1,1)-R(2,2)+R(3,3) 

!!! calculos los autovalores y autovectores:
CC=F
call dsyevx ('V','I','U',4,CC,4,0,0,1,4,0.0d0,num_val&
     &,values,vectors,4,work,8*4,iwork,ifail,info)

!!! computo el rmsd, la matriz de rotacion y g

msd=max(0.0d0,((x_norm+y_norm)-2.0d0*values(4)))/dble(N)
rmsd=sqrt(msd)


call rotation_matrix(vectors(:,4),U)

DO i=1,N
   DO j=1,3
      tmp(:)=matmul(transpose(U(:,:)),y(i,:))
      g(i,j)=(x(i,j)-tmp(j))/(rmsd*dble(N))
   END DO
END DO

!!! calculo las nuevas posiciones con la traslacion y rotacion

!DO i=1,N
!   pos_new(i,:)=matmul(transpose(U(:,:)),struct_2(i,:)-center_2)+center_ref
!END DO


END subroutine min_rmsd


subroutine fit_min_rmsd(N,struct_ref,struct_2,rmsd)

INTEGER,INTENT(IN)::N
DOUBLE PRECISION,DIMENSION(N,3),INTENT(IN)::struct_ref
DOUBLE PRECISION,DIMENSION(N,3),INTENT(INOUT)::struct_2
DOUBLE PRECISION,DIMENSION(3,3)::U
DOUBLE PRECISION,DIMENSION(3)::center_ref,center_2
DOUBLE PRECISION,INTENT(OUT)::rmsd
!DOUBLE PRECISION,DIMENSION(N,2),INTENT(OUT)::g

INTEGER::i,j,k
DOUBLE PRECISION,DIMENSION(N,3)::x,y,pos_new
DOUBLE PRECISION,DIMENSION(N)::w
DOUBLE PRECISION::sw,msd,x_norm,y_norm
DOUBLE PRECISION,DIMENSION(3,3)::R
DOUBLE PRECISION,DIMENSION(4,4)::F
DOUBLE PRECISION,DIMENSION(3)::tmp

!para diagonalizar:
DOUBLE PRECISION,DIMENSION(4,4)::CC
INTEGER::num_val,info
INTEGER, DIMENSION(:),ALLOCATABLE::iwork,ifail
DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE::values,work
DOUBLE PRECISION, DIMENSION(:,:),ALLOCATABLE::vectors


ALLOCATE(values(4),vectors(4,4),work(8*4),iwork(5*4),ifail(4))

w=0.0d0
x=0.0d0
y=0.0d0
CC=0.0d0
rmsd=0.0d0
msd=0.0d0
sw=0.0d0
values=0.0d0
vectors=0.0d0
center_ref=0.0d0
center_2=0.0d0
U=0.0d0
!g=0.0d0
x_norm=0.0d0
y_norm=0.0d0
R=0.0d0
F=0.0d0
tmp=0.0d0
pos_new=0.0d0

!!! copio y peso las coordenadas:
w=1.0d0
DO i=1,N
   sw=w(i)
   x(i,:)=sw*struct_ref(i,:)
   y(i,:)=sw*struct_2(i,:)
END DO

!!! calculo baricentros, centroides y normas:

DO i=1,3
   center_ref(i)=sum(x(:,i))/dble(N)
   center_2(i)=sum(y(:,i))/dble(N)
   x(:,i)=x(:,i)-center_ref(i)
   y(:,i)=y(:,i)-center_2(i)
   x_norm=x_norm+dot_product(x(:,i),x(:,i))
   y_norm=y_norm+dot_product(y(:,i),y(:,i))
END DO

!!! calculo la matriz R
DO i=1,3
   DO j=1,3
      R(i,j)=dot_product(x(:,i),y(:,j))
   END DO
END DO

!!! construimos la matriz F:

F(1,1)=R(1,1)+R(2,2)+R(3,3)
F(2,1)=R(2,3)-R(3,2)
F(3,1)=R(3,1)-R(1,3)
F(4,1)=R(1,2)-R(2,1)
F(1,2)=F(2,1)
F(2,2)=R(1,1)-R(2,2)-R(3,3)
F(3,2)=R(1,2)+R(2,1)
F(4,2)=R(1,3)+R(3,1)
F(1,3)=F(3,1)
F(2,3)=F(3,2)
F(3,3)=-R(1,1)+R(2,2)-R(3,3)
F(4,3)=R(2,3)+R(3,2)
F(1,4)=F(4,1)
F(2,4)=F(4,2)
F(3,4)=F(4,3)
F(4,4)=-R(1,1)-R(2,2)+R(3,3) 

!!! calculos los autovalores y autovectores:
CC=F
call dsyevx ('V','I','U',4,CC,4,0,0,1,4,0.0d0,num_val&
     &,values,vectors,4,work,8*4,iwork,ifail,info)

!!! computo el rmsd, la matriz de rotacion y g

msd=max(0.0d0,((x_norm+y_norm)-2.0d0*values(4)))/dble(N)
rmsd=sqrt(msd)


call rotation_matrix(vectors(:,4),U)

!DO i=1,N
!   DO j=1,3
!      tmp(:)=matmul(transpose(U(:,:)),y(i,:))
!      g(i,j)=(x(i,j)-tmp(j))/(rmsd*dble(N))
!   END DO
!END DO

!!! calculo las nuevas posiciones con la traslacion y rotacion

DO i=1,N
   pos_new(i,:)=matmul(transpose(U(:,:)),struct_2(i,:)-center_2)+center_ref
END DO

struct_2=0.0d0
struct_2=pos_new

END subroutine fit_min_rmsd




subroutine rotation_matrix(q, U)

DOUBLE PRECISION,DIMENSION(:),INTENT(in)::q
DOUBLE PRECISION,DIMENSION(:,:),INTENT(out)::U
DOUBLE PRECISION::q0,q1,q2,q3,b0,b1,b2,b3,q00,q01,q02,q03,q11,q12,q13,q22,q23,q33

q0=q(1)
q1=q(2)
q2=q(3)
q3=q(4)

b0=2.0d0*q0
b1=2.0d0*q1
b2=2.0d0*q2
b3=2.0d0*q3

q00=b0*q0-1.0d0
q01=b0*q1
q02=b0*q2
q03=b0*q3

q11=b1*q1
q12=b1*q2
q13=b1*q3  

q22=b2*q2
q23=b2*q3

q33=b3*q3 

U(1,1)=q00+q11
U(1,2)=q12-q03
U(1,3)=q13+q02

U(2,1)=q12+q03
U(2,2)=q00+q22
U(2,3)=q23-q01

U(3,1)=q13-q02
U(3,2)=q23+q01
U(3,3)=q00+q33

end subroutine rotation_matrix

END MODULE para_rmsd
