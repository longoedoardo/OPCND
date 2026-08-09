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
    ! N = ((ade + 1) * (ade + 2) * (ade + 3) * (ade + 4)) / 24
    !**********************************************************************
    INTEGER                     :: ade
    INTEGER                     :: N, N_mom
    INTEGER                     :: ind_curr(4)

    REAL(dp), ALLOCATABLE       :: XYZWT_tens_ref(:,:)
    REAL(dp), ALLOCATABLE       :: X(:,:)
    REAL(dp), ALLOCATABLE       :: XYZWT_tens(:,:)
    REAL(dp), ALLOCATABLE       :: W(:)
    REAL(dp), ALLOCATABLE       :: fXYZWT(:)
    REAL(dp)                    :: I_4D
 
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
    REAL(dp), ALLOCATABLE       :: XYZWT_ref(:,:)
 
    !**********************************************************************
    ! Regola di cubatura finale (riscalata sul poliedro)
    ! XYZW(k, 1:3) rappresenta le coordinate punto k
    ! XYZW(k, 4) rappresenta il peso finale punto k
    !**********************************************************************
    REAL(dp), ALLOCATABLE       :: XYZWT(:,:)
 
    !**********************************************************************
    ! Matrice di Vandermonde di Chebyshev e coefficienti
    ! V_ref(k, j) rappresenta il valore del j-esimo polinomio nel punto k
    ! coeffs(j) e' la norma^2 del j-esimo polinomio di Chebyshev
    ! moments_ch(j) rappresenta il momento di Chebyshev j sul poliedro
    !**********************************************************************
    REAL(dp), ALLOCATABLE       :: V_ref(:,:)
    REAL(dp), ALLOCATABLE       :: coeffs(:)
    REAL(dp), ALLOCATABLE       :: moments_ch(:)

    !**********************************************************************
! Quadratura di riferimento sul prisma 4D (triangolo x tempo)
! XI_ref, ETA_ref = coordinate baricentriche sul triangolo di riferimento
! T_ref           = coordinata temporale locale in [0,1]
! W_ref           = peso di quadratura combinato (spazio x tempo)
!**********************************************************************
REAL(dp), ALLOCATABLE       :: XI_ref(:)
REAL(dp), ALLOCATABLE       :: ETA_ref(:)
REAL(dp), ALLOCATABLE       :: T_ref(:)
REAL(dp), ALLOCATABLE       :: W_ref(:)

!**********************************************************************
! Parametri della griglia di riferimento sul prisma
! nGP_tri   = punti di Gauss-Jacobi 1D per la base triangolare
! n1D_prism = punti di Gauss-Legendre 1D per l'asse temporale del prisma
! n_tri     = punti totali sul triangolo (= nGP_tri^2)
! n_tot_prism = punti totali della griglia tensoriale (= n_tri * n1D_prism)
!**********************************************************************
INTEGER                     :: nGP_tri
INTEGER                     :: n1D_prism
INTEGER                     :: n_tri
INTEGER                     :: n_tot_prism

END MODULE TypesDef