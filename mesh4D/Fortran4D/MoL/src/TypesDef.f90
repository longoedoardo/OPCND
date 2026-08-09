MODULE TypesDef
    IMPLICIT NONE

!**********************************************************************
! Precisione numerica
!**********************************************************************
    INTEGER, PARAMETER          :: dp = KIND(0.0d0)

!**********************************************************************
! Costanti matematiche
!**********************************************************************
    REAL(dp), PARAMETER         :: PI = 3.14159265358979323846_dp

!**********************************************************************
! Parametri del metodo di cubatura spaziale
! ade    = grado algebrico di esattezza desiderato
! N      = numero di punti della griglia tensoriale di riferimento
! N_mom  = numero di monomi per grado totale <= ade in 3D,
!          pari a (ade+1)(ade+2)(ade+3)/6
!**********************************************************************
    INTEGER                     :: ade
    INTEGER                     :: N, N_mom
    INTEGER                     :: ind_curr(3)

!**********************************************************************
! Discretizzazione temporale (Clenshaw-Curtis su tau in [0,1])
!**********************************************************************
    INTEGER                     :: n_tau
    REAL(dp), ALLOCATABLE       :: tau_nodes(:)
    REAL(dp), ALLOCATABLE       :: tau_weights(:)
    REAL(dp)                    :: current_tau
    REAL(dp), ALLOCATABLE       :: tau_array(:)

!**********************************************************************
! Indici multi-indice di Chebyshev
! chebyshev_indices(k, 1:3) = [i, j, k]
!**********************************************************************
    INTEGER, ALLOCATABLE        :: chebyshev_indices(:,:)

!**********************************************************************
! Griglia di riferimento (Gauss-Chebyshev tensore su cubo di riferimento)
!**********************************************************************
    REAL(dp), ALLOCATABLE       :: XYZW_tens_ref(:,:)
    REAL(dp), ALLOCATABLE       :: XYZW_ref(:,:)
    REAL(dp), ALLOCATABLE       :: X(:,:)

!**********************************************************************
! Matrice di Vandermonde di Chebyshev e quantita' correlate
! V_ref(k, j)      = valore del j-esimo polinomio nel punto k
! coeffs(j)        = norma^2 del j-esimo polinomio di Chebyshev
! moments_ch(j)    = momento di Chebyshev j sul poliedro
!**********************************************************************
    REAL(dp), ALLOCATABLE       :: V_ref(:,:)
    REAL(dp), ALLOCATABLE       :: coeffs(:)
    REAL(dp), ALLOCATABLE       :: moments_ch(:)

!**********************************************************************
! Geometria del poliedro al tempo corrente
!**********************************************************************
    REAL(dp)                    :: bbox_tau(6)
    REAL(dp), ALLOCATABLE       :: vertices_tau(:,:)

!**********************************************************************
! Regola di cubatura riscalata sul poliedro
! XYZW_tens(k, 1:3) = coordinate del punto k (griglia riscalata)
! XYZW(k, 1:3)      = coordinate del punto k (regola finale)
! XYZW(k, 4)        = peso finale del punto k
!**********************************************************************
    REAL(dp), ALLOCATABLE       :: XYZW_tens(:,:)
    REAL(dp), ALLOCATABLE       :: XYZW(:,:)
    REAL(dp), ALLOCATABLE       :: W(:)

!**********************************************************************
! Valutazione della funzione integranda e risultato finale
!**********************************************************************
    REAL(dp), ALLOCATABLE       :: fXYZW(:)
    REAL(dp)                    :: I_4D

!**********************************************************************
! Indici
!**********************************************************************
    INTEGER                     :: kk

END MODULE TypesDef