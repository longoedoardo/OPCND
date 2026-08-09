MODULE TypesDef
    IMPLICIT NONE

    !**********************************************************************
    ! Definizione della precisione doppia (dp)
    !**********************************************************************
    INTEGER, PARAMETER :: dp = KIND(0.0d0)
 
    !**********************************************************************
    ! Costanti matematiche
    !**********************************************************************
    REAL(dp), PARAMETER         :: PI = 3.14159265358979323846_dp
 
    !**********************************************************************
    ! Parametri del metodo di cubatura
    ! ade = grado algebrico di esattezza desiderato
    ! N  = (ade+3)(ade+2)(ade+1)/6
    !**********************************************************************
    INTEGER                     :: ade
    INTEGER                     :: N, N_mom
    INTEGER                     :: ind_curr(3)

    REAL(dp), ALLOCATABLE       :: XYZW_tens_ref(:,:)
    REAL(dp), ALLOCATABLE       :: X(:,:)
    REAL(dp), ALLOCATABLE       :: XYZW_tens(:,:)
    REAL(dp), ALLOCATABLE       :: W(:)
    REAL(dp), ALLOCATABLE       :: fXYZW(:)
    REAL(dp)                    :: I_risultato
 
    !**********************************************************************
    ! Indici chebyshev_indices(k, 1:3) = [i, j, k]
    !**********************************************************************
    INTEGER, ALLOCATABLE        :: chebyshev_indices(:,:)
 
    !**********************************************************************
    ! Bounding box del poliedro
    !**********************************************************************
    REAL(dp), ALLOCATABLE       :: bbox(:)
 
    !**********************************************************************
    ! Regola di cubatura di riferimento (Gauss-Chebyshev tensore)
    !**********************************************************************
    REAL(dp), ALLOCATABLE       :: XYZW_ref(:,:)
 
    !**********************************************************************
    ! Regola di cubatura finale (riscalata sul poliedro)
    ! XYZW(k, 1:3) rappresenta le coordinate punto k
    ! XYZW(k, 4) rappresenta il peso finale punto k
    !**********************************************************************
    REAL(dp), ALLOCATABLE       :: XYZW(:,:)
 
    !**********************************************************************
    ! Matrice di Vandermonde di Chebyshev e coefficienti
    ! V_ref(k, j) rappresenta il valore del j-esimo polinomio nel punto k
    ! coeffs(j) e' la norma^2 del j-esimo polinomio di Chebyshev
    ! moments_ch(j) rappresenta il momento di Chebyshev j sul poliedro
    !**********************************************************************
    REAL(dp), ALLOCATABLE       :: V_ref(:,:)
    REAL(dp), ALLOCATABLE       :: coeffs(:)
    REAL(dp), ALLOCATABLE       :: moments_ch(:)

END MODULE TypesDef