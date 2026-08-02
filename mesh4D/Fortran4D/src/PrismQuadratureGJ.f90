MODULE PrismQuadratureGJ
  USE TypesDef
  
  IMPLICIT NONE

CONTAINS

  SUBROUTINE reference_prism_quadrature(ade, dbox, XI, ETA, T, W_ref, scale_x, num_pts)

   !*******************************************************************************
   ! Genera la griglia di quadratura sul prisma di riferimento (Triangolo x [0,1])
   ! accoppiando una regola per triangoli basata su Gauss-Jacobi (via TriangleQuadraturePoints)
   ! con Gauss-Legendre per l'asse temporale.
   ! Adattato alla logica della versione MATLAB e sfruttando le subroutine del modulo.
   !*******************************************************************************

      IMPLICIT NONE

      ! Input variables
      INTEGER, INTENT(IN)                   :: ade
      REAL(dp), INTENT(IN)                  :: dbox(:, :)

      ! Output variables
      REAL(dp), ALLOCATABLE, INTENT(OUT)    :: XI(:), ETA(:), T(:), W_ref(:)
      REAL(dp), INTENT(OUT)                 :: scale_x
      INTEGER, INTENT(OUT)                  :: num_pts

      ! Local variables
      INTEGER                               :: nGP, n_spazio, n_tempo, i, j, idx
      REAL(dp), ALLOCATABLE                 :: P_std(:, :), W_std(:)
      REAL(dp), ALLOCATABLE                 :: t_time(:), wt_time(:)

      nGP = ade + 1

      ! 1. Quadratura sul triangolo di riferimento (0,0)-(1,0)-(0,1) 
      ! tramite TriangleQuadraturePoints (che usa internamente Gauss-Jacobi)
      n_spazio = nGP * nGP
      ALLOCATE(P_std(2, n_spazio))
      ALLOCATE(W_std(n_spazio))

      CALL TriangleQuadraturePoints(P_std, W_std, n_spazio, nGP)

      ! 2. Generazione nodi e pesi di Gauss-Legendre per l'asse temporale [0, 1]
      n_tempo = ade + 1
      ALLOCATE(t_time(n_tempo), wt_time(n_tempo))
      CALL gauleg(0.0_dp, 1.0_dp, t_time, wt_time, n_tempo)

      ! 3. Prodotto tensoriale Spazio (Triangolo) x Tempo
      num_pts = n_spazio * n_tempo

      IF (ALLOCATED(XI))    DEALLOCATE(XI)
      IF (ALLOCATED(ETA))   DEALLOCATE(ETA)
      IF (ALLOCATED(T))     DEALLOCATE(T)
      IF (ALLOCATED(W_ref)) DEALLOCATE(W_ref)

      ALLOCATE(XI(num_pts))
      ALLOCATE(ETA(num_pts))
      ALLOCATE(T(num_pts))
      ALLOCATE(W_ref(num_pts))

      idx = 1
      DO i = 1, n_spazio
         DO j = 1, n_tempo
            XI(idx)    = P_std(1, i)
            ETA(idx)   = P_std(2, i)
            T(idx)     = t_time(j)
            W_ref(idx) = W_std(i) * wt_time(j)
            idx        = idx + 1
         END DO
      END DO

      ! Fattore di scala fisico per la primitiva lungo X (derivante dalla dbox)
      scale_x = (dbox(2, 1) - dbox(1, 1)) / 2.0_dp

      DEALLOCATE(P_std, W_std, t_time, wt_time)

   END SUBROUTINE reference_prism_quadrature

   SUBROUTINE PrismQuad4D(V, XI_ref, ETA_ref, T_ref, W_ref, XYZTW, WV_CUB)

      IMPLICIT NONE

      REAL(dp), INTENT(IN)                  :: V(:, :)
      REAL(dp), INTENT(IN)                  :: XI_ref(:)
      REAL(dp), INTENT(IN)                  :: ETA_ref(:)
      REAL(dp), INTENT(IN)                  :: T_ref(:)
      REAL(dp), INTENT(IN)                  :: W_ref(:)

      REAL(dp), ALLOCATABLE, INTENT(OUT)    :: XYZTW(:, :)
      REAL(dp), ALLOCATABLE, INTENT(OUT)    :: WV_CUB(:)   ! ora e' la misura di flusso "con segno" (nx * dV)

      INTEGER                               :: num_pts, i
      REAL(dp), ALLOCATABLE                 :: X_base(:, :), X_top(:, :)
      REAL(dp), ALLOCATABLE                 :: g_xi(:, :), g_eta(:, :), g_tau(:, :)
      REAL(dp), ALLOCATABLE                 :: N1(:)        ! componente x del vettore normale generalizzato 4D

      num_pts = SIZE(W_ref)

      IF (ALLOCATED(XYZTW))  DEALLOCATE(XYZTW)
      IF (ALLOCATED(WV_CUB)) DEALLOCATE(WV_CUB)

      ALLOCATE(XYZTW(num_pts, 4))
      ALLOCATE(WV_CUB(num_pts))
      ALLOCATE(X_base(num_pts, 4), X_top(num_pts, 4))
      ALLOCATE(g_xi(num_pts, 4), g_eta(num_pts, 4), g_tau(num_pts, 4))
      ALLOCATE(N1(num_pts))

      DO i = 1, 4
         X_base(:, i) = V(1, i) + XI_ref * (V(2, i) - V(1, i)) + ETA_ref * (V(3, i) - V(1, i))
         X_top(:, i)  = V(4, i) + XI_ref * (V(5, i) - V(4, i)) + ETA_ref * (V(6, i) - V(4, i))
      END DO

      DO i = 1, 4
         XYZTW(:, i) = (1.0_dp - T_ref) * X_base(:, i) + T_ref * X_top(:, i)
      END DO

      ! Vettori tangenti locali (dipendono da T_ref: variano punto per punto se la
      ! mesh si deforma tra t=0 e t=1, cioe' se vertici_1 /= vertici_0)
      DO i = 1, 4
         g_xi(:, i)  = (1.0_dp - T_ref) * (V(2, i) - V(1, i)) + T_ref * (V(5, i) - V(4, i))
         g_eta(:, i) = (1.0_dp - T_ref) * (V(3, i) - V(1, i)) + T_ref * (V(6, i) - V(4, i))
         g_tau(:, i) = X_top(:, i) - X_base(:, i)
      END DO

      ! Componente x (indice 1) del prodotto vettoriale generalizzato 4D di
      ! (g_xi, g_eta, g_tau): usa solo le componenti y,z,tau (indici 2,3,4)
      N1 = g_xi(:,2) * (g_eta(:,3)*g_tau(:,4) - g_eta(:,4)*g_tau(:,3)) &
         - g_xi(:,3) * (g_eta(:,2)*g_tau(:,4) - g_eta(:,4)*g_tau(:,2)) &
         + g_xi(:,4) * (g_eta(:,2)*g_tau(:,3) - g_eta(:,3)*g_tau(:,2))

      ! Peso di flusso locale, con segno: sostituisce sia dV che la vecchia nx costante
      WV_CUB = W_ref * N1

      DEALLOCATE(X_base, X_top, g_xi, g_eta, g_tau, N1)

   END SUBROUTINE PrismQuad4D


  SUBROUTINE TriangleQuadraturePoints(IntGaussP,IntGaussW,nIntGP,nGP) 
    !**********************************************************************
    IMPLICIT NONE
    !**********************************************************************
    ! Argument list declaration                                               
    INTEGER                            :: nIntGP         ! Number of 2D integration points 
    REAL(dp)                           :: IntGaussP(2,nGP*nGP) ! Positions of 2D int. points 
    REAL(dp)                           :: IntGaussW(nGP*nGP)   ! Weights of 2D int. points  
    INTEGER                            :: nGP            ! Number of 1D Gausspoints        
    !**********************************************************************
    ! Local variable declaration                                              
    INTEGER                :: i,j                  ! Loop counters                   
    INTEGER                :: iIntGP               ! Loop counter                    
    REAL(dp)               :: tol                  ! Tolerance of 0.0                
    REAL(dp)               :: mu1(nGP)             ! 1D quadrature positions in y1   
    REAL(dp)               :: mu2(nGP)             ! 1D quadrature positions in y2   
    REAL(dp)               :: A1(nGP)              ! 1D quadrature weights for y1    
    REAL(dp)               :: A2(nGP)              ! 1D quadrature weights for y2    
    !**********************************************************************
    INTENT(IN)             :: nGP 
    INTENT(OUT)            :: IntGaussP, IntGaussW 
    !**********************************************************************
    tol = 1.0 / (10.0**(PRECISION(1.0)-2) )                                   
    !**********************************************************************                                                                         
    ! Quadrature points are defined by the conical product of 1D Gauss-Jacobi 
    ! formulas with nGP quadrature points. See Stround, p. 28ff for details.
    ! **********************************************************************
    nIntGP = nGP*nGP 
    ! 
    CALL gaujac(mu1,A1,nGP,1.0_dp,0.0_dp)     ! Get the Gauss-Jacobi positions and weights
    CALL gaujac(mu2,A2,nGP,0.0_dp,0.0_dp)     ! Get the Gauss-Jacobi positions and weights
    !
    mu1(:) = 0.5*mu1(:) + 0.5       ! Shift and rescale positions, because Stroud
    A1(:)  = 0.5**2*A1(:)           ! integrates over the interval [0,1] and
    mu2(:) = 0.5*mu2(:) + 0.5       ! the function gaujac of the num. recipes
    A2(:)  = 0.5**1*A2(:)           ! integrates over the interval [-1,1].
    !
    iIntGP = 1
    DO i = 1, nGP 
       DO j = 1, nGP
          intGaussP(1,iIntGP) = mu1(i)
          intGaussP(2,iIntGP) = mu2(j)*(1.-mu1(i))
          intGaussW(iIntGP)   = A1(i)*A2(j)      
          iIntGP              = iIntGP + 1
       ENDDO
    ENDDO
    !
    IF (     ((ABS(SUM(intGaussW(:)))-0.5).GT.tol)  ) THEN
       WRITE(*,*) '| Integration points calculated with conical product.'
       WRITE(*,*) '| Number of integration points is  ', nIntGP
       WRITE(*,*) '| SUM of Weights is ', SUM(intGaussW(:)), ' and must be 0.5! '
    END IF
    ! 
  END SUBROUTINE TriangleQuadraturePoints


  PURE ELEMENTAL FUNCTION gammln(xx)
  !**********************************************************************
  IMPLICIT NONE
  !**********************************************************************
  REAL(dp)             :: gammln,xx
  INTEGER          :: j
  REAL(dp)             :: ser,stp,tmp,x,y,cof(6)
  PARAMETER(stp=  2.5066282746310005  )
  PARAMETER(cof= (/           &
       76.18009172947146    , &
       -86.50532032941677   , &
       24.01409824083091    , &
       -1.231739572450155   , &
       .1208650973866179e-2 , &
       -.5395239384953e-5    /)       )
  !**********************************************************************
  INTENT(IN) :: xx
  !**********************************************************************

  x   = xx
  y   = x
  tmp = x+5.5d0
  tmp = (x+0.5d0)*LOG(tmp)-tmp
  ser = 1.000000000190015d0

  DO  j=1,6
     y  = y+1.d0
     ser= ser+cof(j)/y
  END DO

  gammln=tmp+LOG(stp*ser/x)

  RETURN

END FUNCTION gammln

PURE SUBROUTINE gauleg(x1,x2,x,w,n)
  !**********************************************************************
  IMPLICIT NONE
  !**********************************************************************
  INTEGER     ::  n
  REAL(dp)        :: x1,x2,x(n),w(n)
  REAL(dp)        :: EPS
  INTEGER     :: i,j,m
  REAL(dp)        :: p1,p2,p3,pp,xl,xm,z,z1
  !**********************************************************************
  INTENT(IN)  :: x1,x2,n
  INTENT(OUT) :: x,w
  !**********************************************************************
  PARAMETER (EPS=3.E-14)
  !**********************************************************************

  m  = (n+1)/2
  xm = 0.5*(x2+x1)
  xl = 0.5*(x2-x1)

  DO i=1,m
     z = COS(3.141592654*(i-.25)/(n+.5))
1 CONTINUE
     p1 = 1.
     p2 = 0.
     DO j = 1,n
        p3 = p2
        p2 = p1
        p1 = ((2.*j-1.)*z*p2-(j-1.)*p3)/j
     END DO
     pp = n*(z*p1-p2)/(z*z-1.)
     z1 = z
     z  = z1-p1/pp
     IF(ABS(z-z1).GT.EPS)GOTO 1
     x(i)    = xm-xl*z
     x(n+1-i)= xm+xl*z
     w(i)    = 2.*xl/((1.-z*z)*pp*pp)
     w(n+1-i)= w(i)
  END DO
  RETURN
END SUBROUTINE gauleg

SUBROUTINE gaujac(x,w,n,alf,bet)
  !**********************************************************************
  INTEGER ::  n,MAXIT
  REAL(dp)    ::  alf,bet,w(n),x(n)
  REAL(dp)    :: EPS
  INTEGER :: i,its,j,pr
  REAL(dp)    :: alfbet,an,bn,r1,r2,r3
  REAL(dp)    :: a,b,c,p1,p2,p3,pp,temp,z,z1
  REAL(dp)    :: test
  !**********************************************************************
  INTENT(IN)  :: n, alf, bet
  INTENT(OUT) :: x,w
  !**********************************************************************
  test = 1. 
  MAXIT=50
  pr = PRECISION(test)
  eps = 0.1**REAL(pr)
  DO i=1,n
     IF(i.EQ.1)THEN
        an=alf/n
        bn=bet/n
        r1=(1.+alf)*(2.78/(4.+n*n)+.768*an/n)
        r2=1.+1.48*an+.96*bn+.452*an*an+.83*an*bn
        z=1.-r1/r2
     ELSE IF(i.EQ.2)THEN
        r1=(4.1+alf)/((1.+alf)*(1.+.156*alf))
        r2=1.+.06*(n-8.)*(1.+.12*alf)/n
        r3=1.+.012*bet*(1.+.25*ABS(alf))/n
        z=z-(1.-z)*r1*r2*r3
     ELSE IF(i.EQ.3)THEN
        r1=(1.67+.28*alf)/(1.+.37*alf)
        r2=1.+.22*(n-8.)/n
        r3=1.+8.*bet/((6.28+bet)*n*n)
        z=z-(x(1)-z)*r1*r2*r3
     ELSE IF(i.EQ.n-1)THEN
        r1=(1.+.235*bet)/(.766+.119*bet)
        r2=1./(1.+.639*(n-4.)/(1.+.71*(n-4.)))
        r3=1./(1.+20.*alf/((7.5+alf)*n*n))
        z=z+(z-x(n-3))*r1*r2*r3
     ELSE IF(i.EQ.n)THEN
        r1=(1.+.37*bet)/(1.67+.28*bet)
        r2=1./(1.+.22*(n-8.)/n)
        r3=1./(1.+8.*alf/((6.28+alf)*n*n))
        z=z+(z-x(n-2))*r1*r2*r3
     ELSE
        z=3.*x(i-1)-3.*x(i-2)+x(i-3)
     ENDIF
     alfbet=alf+bet
     DO its=1,MAXIT
        temp=2.0+alfbet
        p1=(alf-bet+temp*z)/2.0
        p2=1.0
        DO j=2,n
           p3=p2
           p2=p1
           temp=2*j+alfbet
           a=2*j*(j+alfbet)*(temp-2.0)
           b=(temp-1.0)*(alf*alf-bet*bet+temp*(temp-2.0)*z)
           c=2.0*(j-1+alf)*(j-1+bet)*temp
           p1=(b*p2-c*p3)/a
        END DO
        pp=(n*(alf-bet-temp*z)*p1+2.0*(n+alf)*(n+bet)*p2)/(temp*(1.0-z*z))
        z1=z
        z=z1-p1/pp
        IF(ABS(z-z1).LE.EPS) GOTO 1
     END DO

     stop 'too many iterations in gaujac'
1    x(i)=z
     w(i)=EXP(gammln(alf+n)+gammln(bet+n)-gammln(n+1.0_dp)-gammln(n+alfbet+1.0_dp))*temp*2.0_dp**alfbet/(pp*p2)
  END DO
  RETURN
END SUBROUTINE gaujac

END MODULE PrismQuadratureGJ