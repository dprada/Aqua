MODULE pdb_in_out

IMPLICIT NONE

integer,parameter::long_max=80
integer,parameter::long_corta=30

CONTAINS

SUBROUTINE saco_num_atom(opcion,pdb,N,ios)

INTEGER,INTENT(OUT)::N,ios
CHARACTER(LEN=long_max),INTENT(IN)::pdb,opcion
CHARACTER(LEN=long_corta)::comienzo
INTEGER::ioerror,unit
CHARACTER(LEN=3)::atomname

unit=12

OPEN (unit,FILE=TRIM(pdb),iostat=ios,STATUS="OLD",ACTION="READ")

IF (ios/=0) THEN
   RETURN
END IF

DO !bucle sin fin, termina en ATOM
   READ(unit, 50) comienzo
   IF (comienzo=='ATOM') EXIT
END DO

BACKSPACE(unit)

N=0

DO 

   READ (unit,50,iostat=ioerror) comienzo

   IF (ioerror < 0) exit !se acabo el archivo

   IF (comienzo/='ATOM' .and. comienzo/='ANISOU' .and. comienzo/='HETATM'& 
          .and. comienzo/='SIGATM' .and. comienzo/='SIGUIJ') exit !el archivo no se acabo pero los atomos si

   IF (comienzo=='ATOM' .or. comienzo=='HETATM') THEN

      BACKSPACE(unit)
      read(unit,87) atomname
      
      IF (comienzo/='ATOM') EXIT !no lo cuento si no es atomo
      
      IF ((opcion=='ca'.and.atomname=='CA') .or. &
           (opcion=='backbone'.and. &
           (atomname=='CA'.or.atomname=='N'.or.atomname=='C')).or.&
           (opcion=='backbone2'.and. &
           (atomname=='CA'.or.atomname=='N'.or.atomname=='C'.or.atomname=='CH3')).or.&
           (opcion=='all')) THEN
         N=N+1
      END IF
   END IF

END DO

CLOSE (unit)



50 format (a6)
65 format (a26)
87 format (13x,a3)

END SUBROUTINE saco_num_atom

SUBROUTINE leo_pdb(opcion,pdb,coord,N)

INTEGER,INTENT(OUT)::N
CHARACTER(LEN=long_max),INTENT(IN)::pdb,opcion
DOUBLE PRECISION,DIMENSION(:,:),INTENT(OUT)::coord
CHARACTER(LEN=long_corta)::comienzo
CHARACTER(len=10)::atomname, res_name, cur_res
DOUBLE PRECISION::x,y,z
INTEGER::unit,cuento,ioerror

unit=12
N=0

OPEN(unit,FILE=TRIM(pdb),STATUS="OLD",ACTION="READ")

DO !bucle sin fin, termina en ATOM
   READ(unit, 50) comienzo
   IF (comienzo=='ATOM') EXIT
END DO

BACKSPACE(unit)

cuento=0

DO     !bucle que termina si comienzo no es ATOM

   READ (unit,50,iostat=ioerror) comienzo

   IF (ioerror < 0) exit !se acabo el archivo

   IF (comienzo/='ATOM' .and. comienzo/='ANISOU' .and. comienzo/='HETATM'& 
        .and. comienzo/='SIGATM' .and. comienzo/='SIGUIJ') exit !el archivo no se acabo pero los atomos si

   IF (comienzo=='ATOM' .or. comienzo=='HETATM') THEN
      
      BACKSPACE(unit)
      READ(unit,87) atomname, res_name, cur_res, x, y , z
      
      IF ((opcion=='ca'.and.atomname=='CA') .or. &
           (opcion=='backbone'.and. &
           (atomname=='CA'.or.atomname=='N'.or.atomname=='C')).or.&
           (opcion=='backbone2'.and. &
           (atomname=='CA'.or.atomname=='N'.or.atomname=='C'.or.atomname=='CH3')).or.&
           (opcion=='all')) THEN
         
         cuento=cuento+1
         coord(cuento,1:3) = (/ x, y, z /)
      END IF
   END IF
END DO

N=cuento

close(unit)

50 format (a6)
65 format (a26)
87 format (13x,a3,1x,a3,2x,a4,4x,3f8.3)

END SUBROUTINE leo_pdb

SUBROUTINE indices_pdb(opcion,pdb,indices,labels,N)

INTEGER,INTENT(IN)::N
CHARACTER(LEN=long_max),INTENT(IN)::pdb,opcion
INTEGER,DIMENSION(N),INTENT(OUT)::indices
CHARACTER(LEN=long_corta)::comienzo
CHARACTER(len=10)::atomname, res_name, cur_res
CHARACTER(len=10),DIMENSION(N),INTENT(OUT)::labels
DOUBLE PRECISION::x,y,z
INTEGER::unit,cuento,ioerror,cuento2

unit=12

OPEN(unit,FILE=TRIM(pdb),STATUS="OLD",ACTION="READ")

DO !bucle sin fin, termina en ATOM
   READ(unit, 50) comienzo
   IF (comienzo=='ATOM') EXIT
END DO

BACKSPACE(unit)

cuento=0
cuento2=0

DO     !bucle que termina si comienzo no es ATOM

   READ (unit,50,iostat=ioerror) comienzo

   IF (ioerror < 0) exit !se acabo el archivo

   IF (comienzo/='ATOM' .and. comienzo/='ANISOU' .and. comienzo/='HETATM'& 
        .and. comienzo/='SIGATM' .and. comienzo/='SIGUIJ') exit !el archivo no se acabo pero los atomos si

   IF (comienzo=='ATOM' .or. comienzo=='HETATM') THEN
      
      BACKSPACE(unit)
      READ(unit,87) atomname, res_name, cur_res, x, y , z
      
      cuento2=cuento2+1

      IF ((opcion=='ca'.and.atomname=='CA') .or. &
           (opcion=='backbone'.and. &
           (atomname=='CA'.or.atomname=='N'.or.atomname=='C')).or.&
           (opcion=='backbone2'.and. &
           (atomname=='CA'.or.atomname=='N'.or.atomname=='C'.or.atomname=='CH3')).or.&
           (opcion=='all')) THEN
         
         cuento=cuento+1
         indices(cuento) = cuento2
         labels(cuento) = atomname

      END IF
   END IF
END DO

close(unit)

50 format (a6)
65 format (a26)
87 format (13x,a3,1x,a3,2x,a4,4x,3f8.3)

END SUBROUTINE indices_pdb

SUBROUTINE residuos_pdb(opcion,pdb,indices,N)

INTEGER,INTENT(IN)::N
CHARACTER(LEN=long_max),INTENT(IN)::pdb,opcion
INTEGER,DIMENSION(N),INTENT(OUT)::indices
CHARACTER(LEN=long_corta)::comienzo
CHARACTER(len=10)::atomname, res_name, cur_res
DOUBLE PRECISION::x,y,z
INTEGER::unit,cuento,ioerror,cuento2

unit=12

OPEN(unit,FILE=TRIM(pdb),STATUS="OLD",ACTION="READ")

DO !bucle sin fin, termina en ATOM
   READ(unit, 50) comienzo
   IF (comienzo=='ATOM') EXIT
END DO

BACKSPACE(unit)

cuento=0
cuento2=0

DO     !bucle que termina si comienzo no es ATOM

   READ (unit,50,iostat=ioerror) comienzo

   IF (ioerror < 0) exit !se acabo el archivo

   IF (comienzo/='ATOM' .and. comienzo/='ANISOU' .and. comienzo/='HETATM'& 
        .and. comienzo/='SIGATM' .and. comienzo/='SIGUIJ') exit !el archivo no se acabo pero los atomos si

   IF (comienzo=='ATOM' .or. comienzo=='HETATM') THEN
      
      BACKSPACE(unit)
      READ(unit,87) atomname, res_name, cur_res, x, y , z
      
      IF ((opcion=='ca'.and.atomname=='CA') .or. &
           (opcion=='backbone'.and. &
           (atomname=='CA'.or.atomname=='N'.or.atomname=='C')).or.&
           (opcion=='backbone2'.and. &
           (atomname=='CA'.or.atomname=='N'.or.atomname=='C'.or.atomname=='CH3')).or.&
           (opcion=='all')) THEN
         
         cuento=cuento+1
         READ(cur_res,*) cuento2
         indices(cuento) = cuento2

      END IF
   END IF
END DO

close(unit)

50 format (a6)
65 format (a26)
87 format (13x,a3,1x,a3,2x,a4,4x,3f8.3)

END SUBROUTINE residuos_pdb


SUBROUTINE print_pdb (pdb,coord)

CHARACTER(LEN=long_max),INTENT(IN)::pdb
DOUBLE PRECISION,DIMENSION(:,:),INTENT(IN)::coord

CHARACTER(LEN=long_corta)::comienzo
CHARACTER(len=10)::atomname, res_name
DOUBLE PRECISION::x,y,z
INTEGER::unit,cuento,ioerror,res_no,indice

unit=12


OPEN(unit,FILE=TRIM(pdb),STATUS="OLD",ACTION="READ")

DO !bucle sin fin, termina en ATOM
   READ(unit, 50) comienzo
   IF (comienzo=='ATOM') EXIT
END DO

BACKSPACE(unit)
cuento=0

DO     !bucle que termina si comienzo no es ATOM
   
   READ (unit,50,iostat=ioerror) comienzo
   
   IF (ioerror < 0) exit !se acabo el archivo
   
   IF (comienzo/='ATOM' .and. comienzo/='ANISOU' .and. comienzo/='HETATM'& 
        .and. comienzo/='SIGATM' .and. comienzo/='SIGUIJ') exit !el archivo no se acabo pero los atomos si
   
   IF (comienzo=='ATOM' .or. comienzo=='HETATM') THEN
      
      BACKSPACE(unit)
      READ(unit,87) indice, atomname, res_name, res_no, x, y , z
      
      cuento=cuento+1
      
      PRINT 88,comienzo,indice,atomname,res_name,res_no, coord(cuento,1:3)
      
   END IF

END DO

print 89, comienzo,indice+1,res_name,res_no
print'(a3)','END'

close(unit)

50 format (a6)
65 format (a26)
87 format (6x,i5,1x,a4,1x,a3,2x,i4,4x,3f8.3)
60 format (a26, 4x, 3f8.3)
88 format (a6,i5,1x,a4,1x,a3,2x,i4,4x,3f8.3)
89 format (a6,i5,6x,a3,2x,i4)

END SUBROUTINE print_pdb


END MODULE pdb_in_out
