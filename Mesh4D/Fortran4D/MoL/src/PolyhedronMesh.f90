MODULE PolyhedronMesh

   USE TypesDef

   IMPLICIT NONE

   TYPE :: t_polyhedron
      REAL(dp), ALLOCATABLE :: vertici(:,:)
      REAL(dp), ALLOCATABLE :: vertici_new(:,:)
      INTEGER, ALLOCATABLE  :: facce(:,:)
   END TYPE t_polyhedron

CONTAINS

   SUBROUTINE MeshReader(vertices_file, vertices_new_file, facets_file, poly)
      CHARACTER(LEN=*), INTENT(IN)     :: vertices_file
      CHARACTER(LEN=*), INTENT(IN)     :: vertices_new_file
      CHARACTER(LEN=*), INTENT(IN)     :: facets_file
      TYPE(t_polyhedron), INTENT(OUT)  :: poly

      CALL read_real_matrix_3(vertices_file, poly%vertici)
      CALL read_real_matrix_3(vertices_new_file, poly%vertici_new)
      CALL read_integer_matrix_3(facets_file, poly%facce)

      IF (SIZE(poly%vertici, 1) /= SIZE(poly%vertici_new, 1)) THEN
         ERROR STOP 'ERRORE: il numero di vertici iniziali e finali non coincide.'
      END IF
   END SUBROUTINE MeshReader

   SUBROUTINE read_real_matrix_3(filename, a)
      CHARACTER(LEN=*), INTENT(IN)        :: filename
      REAL(dp), ALLOCATABLE, INTENT(OUT)  :: a(:,:)

      INTEGER                             :: unit
      INTEGER                             :: ios
      INTEGER                             :: n
      INTEGER                             :: i
      REAL(dp)                            :: tmp(3)

      n = count_rows(filename)
      ALLOCATE(a(n, 3))

      OPEN(NEWUNIT=unit, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=ios)
      IF (ios /= 0) ERROR STOP 'ERRORE: impossibile aprire un file di vertici.'

      DO i = 1, n
         READ(unit, *, IOSTAT=ios) tmp
         IF (ios /= 0) ERROR STOP 'ERRORE: impossibile leggere un file di vertici.'
         a(i,:) = tmp
      END DO

      CLOSE(unit)
   END SUBROUTINE read_real_matrix_3

   SUBROUTINE read_integer_matrix_3(filename, a)
      CHARACTER(LEN=*), INTENT(IN)       :: filename
      INTEGER, ALLOCATABLE, INTENT(OUT)  :: a(:,:)

      INTEGER                            :: unit
      INTEGER                            :: ios
      INTEGER                            :: n
      INTEGER                            :: i
      INTEGER                            :: tmp(3)

      n = count_rows(filename)
      ALLOCATE(a(n, 3))

      OPEN(NEWUNIT=unit, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=ios)
      IF (ios /= 0) ERROR STOP 'ERRORE: impossibile aprire un file di facce.'

      DO i = 1, n
         READ(unit, *, IOSTAT=ios) tmp
         IF (ios /= 0) ERROR STOP 'ERRORE: impossibile leggere un file di facce.'
         a(i,:) = tmp
      END DO

      CLOSE(unit)
   END SUBROUTINE read_integer_matrix_3

   INTEGER FUNCTION count_rows(filename) RESULT(n)
      CHARACTER(LEN=*), INTENT(IN) :: filename

      INTEGER                      :: unit
      INTEGER                      :: ios
      CHARACTER(LEN=1024)          :: line

      n = 0
      OPEN(NEWUNIT=unit, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=ios)
      IF (ios /= 0) ERROR STOP 'ERRORE: impossibile aprire un file mesh.'

      DO
         READ(unit, '(A)', IOSTAT=ios) line
         IF (ios /= 0) EXIT
         IF (LEN_TRIM(line) > 0) n = n + 1
      END DO

      CLOSE(unit)
   END FUNCTION count_rows

END MODULE PolyhedronMesh
