! from QE-6.0

subroutine ylmr2 (lmax2, ng, g, gg, ylm)
  !
  !     Real spherical harmonics ylm(G) up to l=lmax
  !     lmax2 = (lmax+1)^2 is the total number of spherical harmonics
  !     Numerical recursive algorithm based on the one given in Numerical
  !     Recipes but avoiding the calculation of factorials that generate
  !     overflow for lmax > 11
  !
  implicit none
  integer, parameter :: DP = 8
  real(DP), parameter :: PI = 4.d0*atan(1.d0)
  real(DP), parameter :: FPI = 4.d0*PI
  !
  integer, intent(in) :: lmax2, ng
  real(DP), intent(in) :: g (3, ng), gg (ng)
  !
  ! BEWARE: gg = g(1)^2 + g(2)^2 +g(3)^2  is not checked on input
  !         incorrect results will ensue if the above does not hold
  !
  real(DP), intent(out) :: ylm (ng,lmax2)
  !
  ! local variables
  !
  REAL(DP), PARAMETER :: eps = 1.0d-9
  REAL(DP), ALLOCATABLE :: cost (:), sent(:), phi (:), Q(:,:,:)
  REAL(DP) :: c, gmod
  INTEGER :: lmax, ig, l, m, lm

  IF (ng < 1 .or. lmax2 < 1) RETURN 

  DO lmax = 0, 25
    IF ((lmax+1)**2 == lmax2) go to 10
  ENDDO
  WRITE(*,*) 'Error ylmr: l > 25 or wrong number of Ylm required, lmax, lmax2 = ', lmax, lmax2
  STOP

10 CONTINUE
  !
  IF (lmax == 0) THEN
     ylm(:,1) =  sqrt (1.d0 / fpi)
     RETURN 
  ENDIF 

  !
  !  theta and phi are polar angles, cost = cos(theta)
  !
  allocate(cost(ng), sent(ng), phi(ng), Q(ng,0:lmax,0:lmax) )
  Q(:,:,:) = 0.d0 ! ffr
  !

  do ig = 1, ng
     gmod = sqrt (gg (ig) )
     if (gmod < eps) then
        cost(ig) = 0.d0
     else
        cost(ig) = g(3,ig)/gmod
     endif
     !
     !  beware the arc tan, it is defined modulo pi
     !
     if (g(1,ig) > eps) then
        phi(ig) = atan( g(2,ig)/g(1,ig) )
     else if (g(1,ig) < -eps) THEN
        phi(ig) = atan( g(2,ig)/g(1,ig) ) + pi
     else
        phi(ig) = sign( pi/2.d0,g(2,ig) )
     end if
     sent(ig) = sqrt(max(0d0,1.d0-cost(ig)**2))
     
     WRITE(*,'(1x,A,F18.10)') 'phi  = ', phi(ig)
     WRITE(*,'(1x,A,F18.10)') 'cost = ', cost(ig)
     WRITE(*,'(1x,A,F18.10)') 'sint = ', sent(ig)

  ENDDO
  !
  !  Q(:,l,m) are defined as sqrt ((l-m)!/(l+m)!) * P(:,l,m) where
  !  P(:,l,m) are the Legendre Polynomials (0 <= m <= l)
  !
  lm = 0
  
  DO l = 0, lmax
     
     WRITE(*,'(1x,A,I4)',advance='no') 'l = ', l

     c = sqrt (DBLE(2*l+1) / fpi)
     if ( l == 0 ) then
        do ig = 1, ng
           Q (ig,0,0) = 1.d0
        end do
     else if ( l == 1 ) then
        do ig = 1, ng
           Q (ig,1,0) = cost(ig)
           Q (ig,1,1) =-sent(ig)/sqrt(2.d0)
        end do
     else
        !
        !  recursion on l for Q(:,l,m)
        !
        do m = 0, l - 2
           do ig = 1, ng
              Q(ig,l,m) = cost(ig)*(2*l-1)/sqrt(DBLE(l*l-m*m))*Q(ig,l-1,m) &
                       - sqrt(DBLE((l-1)*(l-1)-m*m))/sqrt(DBLE(l*l-m*m))*Q(ig,l-2,m)
           end do
        end do
        do ig = 1, ng
           Q(ig,l,l-1) = cost(ig) * sqrt(DBLE(2*l-1)) * Q(ig,l-1,l-1)
        end do
        do ig = 1, ng
           Q(ig,l,l)   = - sqrt(DBLE(2*l-1))/sqrt(DBLE(2*l))*sent(ig)*Q(ig,l-1,l-1)
        end do
     end if
     !
     ! Y_lm, m = 0
     !
     lm = lm + 1
     WRITE(*,*) 'm = 0: lm = ', lm
     do ig = 1, ng
        ylm(ig, lm) = c * Q(ig,l,0)
     end do
     !
     do m = 1, l
        !
        ! Y_lm, m > 0
        !
        lm = lm + 1
        WRITE(*,*) 'm > 0: lm = ', lm
        do ig = 1, ng
           ylm(ig, lm) = c * sqrt(2.d0) * Q(ig,l,m) * cos (m*phi(ig))
        end do
        !
        ! Y_lm, m < 0
        !
        lm = lm + 1
        WRITE(*,*) 'm < 0: lm = ', lm
        DO ig = 1, ng
           ylm(ig, lm) = c * sqrt(2.d0) * Q(ig,l,m) * sin (m*phi(ig))
        ENDDO
     ENDDO 
  ENDDO 
  !
  DEALLOCATE(cost, sent, phi, Q)
  !
  RETURN

END SUBROUTINE  ylmr2

