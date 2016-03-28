PROGRAM INTEG_RDF

IMPLICIT NONE

INTEGER::i,j,k,IOS
INTEGER::num_args,num_wats
CHARACTER*40::filename,buffer,option
REAL::pi
DOUBLE PRECISION::density
DOUBLE PRECISION::lim_inf,lim_sup
DOUBLE PRECISION::x,y,z,delta_x,integral
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::rdf

pi=acos(-1.0d0)

ALLOCATE(rdf(3,3000))
rdf=0.0d0

num_args=IARG()

lim_inf=0.0d0
lim_sup=0.0d0

IF ((num_args/=5)) THEN
   PRINT*, ' '
   PRINT*, 'integ_rdf -option filename lim_min lim_max num_wats'
   PRINT*, ' '
   PRINT*, ' -option  :  -whole, -shell, -inclusions'
   PRINT*, ' filename :  name of the file'
   PRINT*, ' lim_min  :  inferior limit for the integral (real)'
   PRINT*, ' lim_sup  :  superior limit for the integral (real)'
   PRINT*, ' num_wats :  number of water molecules in the box (integer)'
   PRINT*, ' '
   stop
END IF

CALL getarg(1,buffer)
option=TRIM(buffer)

CALL getarg(2,buffer)
filename=TRIM(buffer)

CALL getarg(3,buffer)
READ(buffer,*) lim_inf

CALL getarg(4,buffer)
READ(buffer,*) lim_sup

CALL getarg(5,buffer)
READ(buffer,*) num_wats


IF (option=='-whole') THEN

   OPEN(12,FILE=filename,STATUS='OLD',ACTION='READ')
   
   i=0
   DO
      
      READ(12,*,IOSTAT=IOS) x, y ,z
      IF (IOS/=0) EXIT
      i=i+1
      rdf(1,i)=x
      rdf(2,i)=y
      
   END DO
   
   delta_x=0.0d0
   DO j=1,100
      delta_x=delta_x+rdf(1,j+1)-rdf(1,j)
   END DO

   density=0.0d0

   DO j=1,i
      integral=integral+rdf(2,j)*delta_x*4*pi*rdf(1,j)**2
   END DO
   density=num_wats/integral

   integral=0.0d0
   
   DO j=1,i
      IF ((rdf(1,j)>lim_inf).and.(rdf(1,j)<lim_sup)) THEN
         integral=integral+density*rdf(2,j)*delta_x*4*pi*rdf(1,j)**2
      END IF
   END DO
   
   PRINT*,integral
   

END IF

IF (option=='-shell') THEN



   OPEN(12,FILE=filename,STATUS='OLD',ACTION='READ')

   i=0
   DO
      
      READ(12,*,IOSTAT=IOS) x, y, z
      IF (IOS/=0) EXIT
      i=i+1
      rdf(1,i)=x
      rdf(2,i)=y
      rdf(3,i)=z
   END DO
   
   delta_x=0.0d0
   DO j=1,100
      delta_x=delta_x+rdf(1,j+1)-rdf(1,j)
   END DO
   
   density=0.0d0
   
   DO j=1,i
      integral=integral+rdf(2,j)*delta_x*4*pi*rdf(1,j)**2
   END DO
   density=num_wats/integral

   integral=0.0d0
   
   DO j=1,i
      IF ((rdf(1,j)>lim_inf).and.(rdf(1,j)<lim_sup)) THEN
         integral=integral+density*rdf(3,j)*delta_x*4*pi*rdf(1,j)**2
      END IF
   END DO
   
   PRINT*,integral

END IF

IF (option=='-inclusions') THEN



   OPEN(12,FILE=filename,STATUS='OLD',ACTION='READ')

   i=0
   DO
      
      READ(12,*,IOSTAT=IOS) x, y, z
      IF (IOS/=0) EXIT
      i=i+1
      rdf(1,i)=x
      rdf(2,i)=y
      rdf(3,i)=z
   END DO
   
   delta_x=0.0d0
   DO j=1,100
      delta_x=delta_x+rdf(1,j+1)-rdf(1,j)
   END DO
   
   density=0.0d0
   
   DO j=1,i
      integral=integral+rdf(2,j)*delta_x*4*pi*rdf(1,j)**2
   END DO
   density=num_wats/integral

   integral=0.0d0
   
   DO j=1,i
      IF ((rdf(1,j)>lim_inf).and.(rdf(1,j)<lim_sup)) THEN
         integral=integral+density*(rdf(2,j)-rdf(3,j))*delta_x*4*pi*rdf(1,j)**2
      END IF
   END DO
   
   PRINT*,integral

END IF



END PROGRAM INTEG_RDF
