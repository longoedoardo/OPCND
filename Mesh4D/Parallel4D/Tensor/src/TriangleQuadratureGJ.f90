MODULE TriangleQuadratureGJ

  USE TypesDef

  !**********************************************************************
  ! Quadratura di Gauss-Jacobi su un triangolo generico (3D), ottenuta
  ! mediante shifting affine dei punti/pesi calcolati sul triangolo di
  ! riferimento (0,0), (1,0), (0,1) tramite prodotto conico di formule
  ! 1D di Gauss-Jacobi.
  !**********************************************************************

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE TriangleQuadratureGJPoints(IntGaussP,IntGaussW,nIntGP,nGP) 

      IMPLICIT NONE

      !**********************************************************************
      ! Argomenti
      !**********************************************************************
      INTEGER, INTENT(IN) :: nGP
      INTEGER, INTENT(OUT) :: nIntGP
      REAL(dp), INTENT(OUT) :: IntGaussP(nGP*nGP, 2)
      REAL(dp), INTENT(OUT) :: IntGaussW(nGP*nGP)
      !**********************************************************************
      ! Variabili locali
      !**********************************************************************
      INTEGER :: i, j
      INTEGER :: iIntGP

      REAL(dp) :: tol
      REAL(dp) :: mu1(nGP)
      REAL(dp) :: mu2(nGP)
      REAL(dp) :: A1(nGP)
      REAL(dp) :: A2(nGP)
      tol = 1.0 / (10.0**(PRECISION(1.0)-2) )                                   
      !**********************************************************************

      nIntGP = nGP*nGP

      ! Calcolo dei punti e dei pesi di Gauss-Jacobi
      CALL gaujac(mu1, A1, nGP, 1.0_dp, 0.0_dp)
      CALL gaujac(mu2, A2, nGP, 0.0_dp, 0.0_dp)

      ! Trasformazione nell'intervallo [0,1]
      mu1(:) = 0.5_dp * mu1(:) + 0.5_dp
      A1(:)  = 0.5_dp**2 * A1(:)

      mu2(:) = 0.5_dp * mu2(:) + 0.5_dp
      A2(:)  = 0.5_dp * A2(:)

      ! Costruzione della quadratura sul triangolo di riferimento
      iIntGP = 1
      DO i = 1, nGP
         DO j = 1, nGP
            IntGaussP(iIntGP,1) = mu1(i)
            IntGaussP(iIntGP,2) = mu2(j) * (1.0_dp - mu1(i))
            IntGaussW(iIntGP) = A1(i) * A2(j)
            iIntGP = iIntGP + 1
         END DO
      END DO

      !**********************************************************************
      ! Controllo della somma dei pesi
      !**********************************************************************

      IF (ABS(SUM(IntGaussW) - 0.5_dp) > tol) THEN

         WRITE(*,*) '| Integration points calculated with conical product.'
         WRITE(*,*) '| Number of integration points is ', nIntGP
         WRITE(*,*) '| SUM of Weights is ', SUM(IntGaussW),' and must be 0.5!'
      END IF

   END SUBROUTINE TriangleQuadratureGJPoints


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
  INTEGER ::  MAXIT
  INTEGER, INTENT(IN) :: n
   REAL(dp), INTENT(IN) :: alf, bet
   REAL(dp), INTENT(OUT) :: x(n), w(n)
  REAL(dp)    :: EPS
  INTEGER :: i,its,j,pr
  REAL(dp)    :: alfbet,an,bn,r1,r2,r3
  REAL(dp)    :: a,b,c,p1,p2,p3,pp,temp,z,z1
  REAL(dp)    :: test
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
        z=3.*x(MAX(i-1,1))-3.*x(MAX(i-2,1))+x(MAX(i-3,1))
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
     w(i)=EXP(LOG_GAMMA(alf+n)+LOG_GAMMA(bet+n)-LOG_GAMMA(n+1.0_dp)-LOG_GAMMA(n+alfbet+1.0_dp))*temp*2.0_dp**alfbet/(pp*p2)
  END DO
  RETURN
END SUBROUTINE gaujac

END MODULE TriangleQuadratureGJ
