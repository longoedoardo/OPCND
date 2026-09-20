PROGRAM example_convex
   USE TypesDef, ONLY: dp
   USE PolyhedronMesh, ONLY: t_polyhedron, MeshReader
   USE OPC4D_Tensor_Module, ONLY: OPC4D_Tensor
   IMPLICIT NONE

   TYPE(t_polyhedron) :: poly
   REAL(dp), ALLOCATABLE :: XYZT(:,:), weights(:)
   INTEGER :: degree, c0, c1, rate
   REAL(dp) :: elapsed, integral

   degree = 1
   CALL MeshReader('convex_vertex.dat', 'convex_vertex_new.dat', 'convex_tri.dat', poly)

   WRITE(*,'(/,A)') '**************************************************************'
   WRITE(*,'(A)')   '                 OPTIMALPOLYCUBA4D Tensor'
   WRITE(*,'(A)')   '    Cubatura su Dominio Poliedrale Convesso in Movimento'
   WRITE(*,'(A,/)') '**************************************************************'
   WRITE(*,'(A,I0)') 'Ade:                   ', degree
   WRITE(*,'(A,I0)') 'Numero di vertici:     ', SIZE(poly%vertici,1)
   WRITE(*,'(A,I0)') 'Numero di facce:       ', SIZE(poly%facce,1)

   CALL SYSTEM_CLOCK(c0, rate)
   CALL OPC4D_Tensor(degree, poly%vertici, poly%vertici_new, poly%facce, 'GJCC', XYZT, weights)
   CALL SYSTEM_CLOCK(c1, rate)
   elapsed = REAL(c1 - c0, dp) / REAL(rate, dp)
   integral = SUM(weights)

   WRITE(*,'(A,I0)')       'Numero di nodi 4D:    ', SIZE(XYZT,1)
   WRITE(*,'(A,ES22.15)')  'Integrale numerico:   ', integral
   WRITE(*,'(A,ES12.6,A)') 'Tempo di calcolo:     ', elapsed, ' s'
END PROGRAM example_convex
