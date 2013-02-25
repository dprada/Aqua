SUBROUTINE open_read(len_ch,file_name,funit,o_vars,o_natom,o_delta_t,pos_o)

  IMPLICIT NONE
  INTEGER,INTENT(IN)::len_ch
  CHARACTER(80),INTENT(IN)::file_name
  INTEGER,INTENT(OUT)::funit,o_natom
  INTEGER(KIND=8),INTENT(OUT)::pos_o
  INTEGER,DIMENSION(20),INTENT(OUT)::o_vars
  DOUBLE PRECISION,INTENT(OUT)::o_delta_t

  LOGICAL:: UNITOP

  INTEGER::uno
  CHARACTER(1),DIMENSION(20)::word_uno
  equivalence (word_uno, uno)
  CHARACTER(6)::endianness

  INTEGER(KIND=4)::NTITLE, DCD_VARS(20),NATOM,HD(2)
  CHARACTER(4)::DCD_TYPE
  CHARACTER(80),ALLOCATABLE,DIMENSION(:)::TITLE
  REAL(KIND=4)::delta_t
  equivalence(delta_t,DCD_VARS(10))
  integer(4), parameter :: dcd_magic_little = X'44524f43', dcd_magic_big = X'434f5244'


  UNITOP=.False.

  funit=500

  do while (UNITOP)
     inquire (unit=funit,opened=UNITOP)
     if (UNITOP) funit=funit+1
  end do

  if (len_ch>80) then
     PRINT*, '# Error: Name of file too long.'
  end if

  OPEN (unit=funit,file=TRIM(file_name),status='old',action='read',form='unformatted',access='stream')

  !!!!!!!!! Checking the native Endianness

  uno=1
  IF (ichar(word_uno(1)) == 0) THEN
     !"Big Endian"
     endianness='BIG   '
  ELSE
     !"Little Endian"
     endianness='LITTLE'
  END IF
  !!!!!!!!!

  !!!!!!!!! Checking the file Endianness

  READ(funit) HD(:)
  IF (HD(1)==84) THEN
     IF (HD(2) == dcd_magic_little) THEN
        ! Detected standard 32-bit DCD file with little-endian data.
        IF (endianness=='BIG   ') THEN
           CLOSE(funit)
           OPEN (unit=funit,file='run.dcd',status='old',action='read',form='unformatted',access='stream',convert='LITTLE_ENDIAN')
           READ(funit) HD(:)
        END IF
     END IF
     IF (HD(2) == dcd_magic_big) THEN
        ! Detected standard 32-bit DCD file with big-endian data.
        IF (endianness=='LITTLE') THEN
           CLOSE(funit)
           OPEN (unit=funit,file='run.dcd',status='old',action='read',form='unformatted',access='stream',convert='BIG_ENDIAN')
           READ(funit) HD(:)
        END IF
     END IF
  ELSE
     IF ((HD(1)+HD(2))==84) THEN
        ! Detected CHARMM -i8 64-bit DCD file not supported.
        print*,'# Detected CHARMM -i8 64-bit DCD file. File not supported.'
        stop
     ELSE
        print*,'# Detected Unknown DCD format. File not supported.'
        stop
     END IF
  END IF

  !!!!!!!!!! Reading header

  READ(funit,pos=5) DCD_TYPE,DCD_VARS
  READ(funit,pos=97) NTITLE
  ALLOCATE(TITLE(NTITLE))
  READ(funit) TITLE(:)
  READ(funit) HD,NATOM
  
  IF (DCD_TYPE/='CORD') THEN
     print*,'# DCD type ',DCD_TYPE,' not supported.'
     stop
  END IF

  ! DCD_TYPE: CORD or VELD
  ! DCD_VARS(1):  Number of frames in this file
  ! DCD_VARS(2):  Number of previous integration steps
  ! DCD_VARS(3):  Frequency (integration steps) to save this file
  ! DCD_VARS(4):  Number of integration steps in the run to create this file
  ! DCD_VARS(5):  Frequency of coordinate saving
  ! DCD_VARS(6):  
  ! DCD_VARS(7):
  ! DCD_VARS(8):  Number of degrees of freedom during the run
  ! DCD_VARS(9):  Number of fixed atoms
  ! DCD_VARS(10): Timestep in AKMA-units. Bit-copy from the 32-bit real number
  ! DCD_VARS(11): 1 if crystal lattice information is present in the frames
  ! DCD_VARS(12): 1 if this is a 4D trajectory
  ! DCD_VARS(13): 1 if fluctuating charges are present
  ! DCD_VARS(14): 1 if trajectory is the result of merge without consistency checks
  ! DCD_VARS(15):
  ! DCD_VARS(16):
  ! DCD_VARS(17):
  ! DCD_VARS(18):
  ! DCD_VARS(19):
  ! DCD_VARS(20): CHARMM version number

  DEALLOCATE(TITLE)

  ! Output:
  INQUIRE(funit,pos=pos_o)
  o_natom=NATOM
  o_vars=DCD_VARS
  o_delta_t=dble(delta_t)

END SUBROUTINE open_read


SUBROUTINE read (funit,natom,with_cell,pos_i,pos_o,cell,coors,io_err,io_end)

  IMPLICIT NONE
  INTEGER,INTENT(IN)::funit,natom
  INTEGER(KIND=8),INTENT(IN)::pos_i
  INTEGER,INTENT(IN):: with_cell
  INTEGER,INTENT(OUT)::io_err,io_end
  INTEGER(KIND=8),INTENT(OUT)::pos_o

  DOUBLE PRECISION,DIMENSION(3,3),INTENT(OUT)::cell
  DOUBLE PRECISION,DIMENSION(natom,3),INTENT(OUT)::coors

  REAL*8::buffer_cell(3,3)
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:)::buffer
  INTEGER(KIND=4)::HD(2)


  cell=0.0d0
  pos_o=pos_i
  io_err=0
  io_end=1

  ALLOCATE(buffer(natom))

  IF (with_cell==1) THEN
     buffer_cell=0.0d0
     READ(funit,pos=pos_o,end=600) HD,buffer_cell(1,1), buffer_cell(1,2), buffer_cell(2,2), &
          buffer_cell(1,3), buffer_cell(2,3), buffer_cell(3,3)
     cell=dble(buffer_cell)
     INQUIRE(funit,pos=pos_o)
  END IF

  coors=0.0d0
  READ(funit,pos=pos_o,end=600) HD,buffer(:)
  coors(:,1)=dble(buffer(:))
  READ(funit) HD,buffer(:)
  coors(:,2)=dble(buffer(:))
  READ(funit) HD,buffer(:)
  coors(:,3)=dble(buffer(:))

  io_end=0

  INQUIRE(funit,pos=pos_o)
  
600  DEALLOCATE(buffer)

END SUBROUTINE read

SUBROUTINE open_write(len_ch,file_name,i_vars,i_natom,i_delta_t,origin_name,funit)

  IMPLICIT NONE
  INTEGER,INTENT(IN)::len_ch,i_natom
  CHARACTER(80),INTENT(IN)::file_name,origin_name
  INTEGER,DIMENSION(20),INTENT(IN)::i_vars
  DOUBLE PRECISION,INTENT(IN)::i_delta_t
  INTEGER,INTENT(OUT)::funit

  LOGICAL:: UNITOP

  INTEGER(KIND=4)::NTITLE, DCD_VARS(20),NATOM,HD(2)
  CHARACTER(4)::DCD_TYPE
  CHARACTER(80),ALLOCATABLE,DIMENSION(:)::TITLE
  REAL(KIND=4)::delta_t
  equivalence(DCD_VARS(10),delta_t)

  UNITOP=.False.
  funit=550
  do while (UNITOP)
     inquire (unit=funit,opened=UNITOP)
     if (UNITOP) funit=funit+1
  end do

  if (len_ch>80) then
     PRINT*, '# Error: Name of file too long.'
     return
  end if

  OPEN (unit=funit,file=TRIM(file_name),status='NEW',action='write',form='unformatted')

  DCD_TYPE='CORD'
  DCD_VARS=i_vars
  delta_t=i_delta_t
  NATOM=i_natom
  NTITLE=2
  ALLOCATE(TITLE(2))
  TITLE(1)='REMARK TRAJECTORY CREATED BY PYNORAMIX 0.1'
  TITLE(2)='REMARK FROM THE ORIGINAL TRAJECTORY NAMED '//TRIM(origin_name)

  WRITE(funit) DCD_TYPE, DCD_VARS
  WRITE(funit) NTITLE,TITLE(:)
  WRITE(funit) NATOM

END SUBROUTINE open_write

SUBROUTINE write(funit,cell,coors,i_natom)

  IMPLICIT NONE
  INTEGER,INTENT(IN)::funit,i_natom
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::cell
  DOUBLE PRECISION,DIMENSION(i_natom,3),INTENT(IN)::coors

  REAL*8,ALLOCATABLE,DIMENSION(:,:)::cell_buffer
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:)::buffer

  ALLOCATE(buffer(i_natom),cell_buffer(3,3))
  cell_buffer=cell                      !!! MMM.... There is something to be fixed here

  WRITE(funit) cell_buffer(1,1), cell_buffer(1,2), cell_buffer(2,2), cell_buffer(1,3), cell_buffer(2,3), cell_buffer(3,3)

  buffer=real(coors(:,1))
  WRITE(funit) buffer(:)
  buffer=real(coors(:,2))
  WRITE(funit) buffer(:)
  buffer=real(coors(:,3))
  WRITE(funit) buffer(:)

  DEALLOCATE(buffer)

END SUBROUTINE write

SUBROUTINE close(funit,io_err)

  IMPLICIT NONE
  INTEGER,INTENT(IN)::funit
  INTEGER,INTENT(OUT)::io_err

  CLOSE(funit)
  io_err=0 !good

END SUBROUTINE close

SUBROUTINE close_write(funit,i_vars,i_natom,i_delta_t,io_err)

  IMPLICIT NONE
  INTEGER,INTENT(IN)::funit,i_natom
  CHARACTER(80)::file_name
  INTEGER,DIMENSION(20),INTENT(IN)::i_vars
  DOUBLE PRECISION,INTENT(IN)::i_delta_t
  INTEGER,INTENT(OUT)::io_err

  INTEGER(KIND=4)::NTITLE, DCD_VARS(20),NATOM,HD(2)
  CHARACTER(4)::DCD_TYPE
  CHARACTER(80),ALLOCATABLE,DIMENSION(:)::TITLE
  REAL(KIND=4)::delta_t
  equivalence(DCD_VARS(10),delta_t)

  DCD_VARS=i_vars
  delta_t=i_delta_t

  INQUIRE(funit,name=file_name)
  CLOSE(funit)

  OPEN (unit=funit,file=TRIM(file_name),status='OLD',action='readwrite',form='unformatted',access='stream')
  WRITE(funit,pos=9) DCD_VARS(:)
  CLOSE (funit)

END SUBROUTINE close_write


  
