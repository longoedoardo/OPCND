PROGRAM example_concave

   USE TypesDef
   USE PolyhedronMesh
   USE OPC4D_MoL_Module

   IMPLICIT NONE

   TYPE(t_polyhedron)                 :: poly
   REAL(dp), ALLOCATABLE              :: XYZtau(:,:)
   REAL(dp), ALLOCATABLE              :: W(:)
   INTEGER                            :: ade
   INTEGER                            :: n_tau
   INTEGER                            :: c0
   INTEGER                            :: c1
   INTEGER                            :: rate
   REAL(dp)                           :: elapsed
   REAL(dp)                           :: integral

   ade = 1
   n_tau = 1000

   WRITE(*,'(/,A)') '**************************************************************'
   WRITE(*,'(A)')   '                       OPC4D_MoL'
   WRITE(*,'(A)')   '                - Metodo delle Linee -'
   WRITE(*,'(A)')   '    Cubatura su Dominio Poliedrale Concavo in Movimento'
   WRITE(*,'(A,/)') '**************************************************************'

   CALL MeshReader('concave_vertex.dat', 'concave_vertex_new.dat', 'concave_tri.dat', poly)

   WRITE(*,'(A,I0)') 'Ade:                   ', ade
   WRITE(*,'(A,I0)') 'Numero di nodi tau:    ', n_tau
   WRITE(*,'(A,I0)') 'Numero di vertici:     ', SIZE(poly%vertici, 1)
   WRITE(*,'(A,I0)') 'Numero di facce:       ', SIZE(poly%facce, 1)
   WRITE(*,'(A)')    'Funzione integranda:   f(x,y,z,tau) = 1'

   WRITE(*,'(/,A)') 'Inizio Cubatura...'
   CALL SYSTEM_CLOCK(c0, rate)
   CALL OPC4D_MoL(ade, n_tau, poly%vertici, poly%vertici_new, poly%facce, 'GJ', XYZtau, W)
   CALL SYSTEM_CLOCK(c1, rate)
   elapsed = REAL(c1 - c0, dp) / REAL(rate, dp)

   integral = SUM(W)

   WRITE(*,'(A)') 'Fine Cubatura...'
   WRITE(*,'(/,A,I0)')       'Numero di nodi 4D:    ', SIZE(XYZtau, 1)
   WRITE(*,'(A,ES22.15)')    'Integrale numerico:   ', integral
   WRITE(*,'(A,ES12.6,A)')   'Tempo di calcolo:     ', elapsed, ' s'

END PROGRAM example_concave
