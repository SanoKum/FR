! 1D shock tube using Flux Reconstruction Method 
! coded by Kumpei Sano. 2020 03 27

module lib4FR
    contains
    subroutine lagrange_poly(nSolPt, u_k, xi_k, x, output)
        implicit none

        integer,intent(in) :: nSolPt
        real(8),intent(in) :: u_k  (nSolPt)
        real(8),intent(in) :: xi_k(nSolPt)
        real(8),intent(in) :: x
        real(8),intent(out):: output
        real(8)            :: fai(nSolPt)
        integer            :: i , j


        do i=1,nSolPt
            fai(i) = 1.d0
            do j=1,nSolPt
                if (i.eq.j) cycle
                fai(i) = fai(i)*(x - xi_k(j))/(xi_k(i) - xi_k(j))
            end do
            !write(*,*) "fai(i)=",i,fai(i)
        end do

        output = dot_product(u_k(1:nSolPt),fai(1:nSolPt))

    end subroutine lagrange_poly

    subroutine calcA( rho, u, e, gam, A )
        implicit none

        real(8)             :: rho   , u , e , gam
        real(8)             :: H     , p
        real(8),intent(out) :: A(3,3)

        p = (gam-1.0)*(e-(rho*u**2)*0.5)
        H = (e + p)/rho

        A(1,1) = 0.0 ; A(1,2) = 1.0 ; A(1,3) = 0.0
        A(2,1) = (gam-3.0)*0.5*u**2 
        A(2,2) =  (3.0-gam)*u
        A(2,3) =  gam - 1.0
        A(3,1) = ((gam-1.0)*0.5*u**2 - H)*u
        A(3,2) = H - (gam-1.0)*u**2
        A(3,3) = gam*u

    end subroutine calcA

    subroutine calcA2( rho, u, H, gam, A )
        implicit none

        real(8)             :: rho   , u , e , gam
        real(8)             :: H     , p
        real(8),intent(out) :: A(3,3)

        A(1,1) = 0.0 ; A(1,2) = 1.0 ; A(1,3) = 0.0
        A(2,1) = (gam-3.0)*0.5*u**2 
        A(2,2) =  (3.0-gam)*u
        A(2,3) =  gam - 1.0
        A(3,1) = ((gam-1.0)*0.5*u**2 - H)*u
        A(3,2) = H - (gam-1.0)*u**2
        A(3,3) = gam*u

    end subroutine calcA2

    subroutine calcE( Q, gam, E)
        implicit none

        real(8)             :: Q(3) , gam
        real(8),intent(out) :: E(3)
        real(8)             :: rho , u , energy , p

        rho      = Q(1)
        u        = Q(2)/rho
        energy   = Q(3)

        p = (gam-1.0)*(energy - (rho*u**2)*0.5)

        E(1) = Q(2)
        E(2) = rho*u**2 + p
        E(3) = (energy+p)*u

    end subroutine calcE

    subroutine calcR( rho, u, H, gam, R )
        implicit none

        real(8)             :: rho , u , e , gam
        real(8)             :: H   , c
        real(8),intent(out) :: R(3,3)

        !p = (gam-1.0)*(e-(rho*u**2)*0.5)
        !H = (e + p)/rho
        c = (gam-1.0)*(H-0.5*u**2)
        if (c .le.0.0) then
            c = 0.0
        else
            c = dsqrt(c)
        end if

        R(1,1) = 1.0   ; R(1,2) = 1.0      ; R(1,3) = 1.0
        R(2,1) = u-c   ; R(2,2) = u        ; R(2,3) = u+c
        R(3,1) = H-c*u ; R(3,2) = 0.5*u**2 ; R(3,3) = H+c*u
    end subroutine calcR

    subroutine calcR_inv( rho, u, H, gam, R_inv )
        implicit none

        real(8)             :: rho , u , e , gam
        real(8)             :: H   ,  c
        real(8)             :: ita
        real(8),intent(out) :: R_inv(3,3)

        c = (gam-1.0)*(H-0.5*u**2)

        if (c .le.0.0) then
            c = 0.0
        else
            c = dsqrt(c)
        end if

        ita = (gam-1.0)/(c**2)

        R_inv(1,1) = (0.5*ita*u**2+u/c)*0.5   ; R_inv(1,2) = -(ita*u+1.d0/c)*0.5 ; R_inv(1,3) = ita*0.5
        R_inv(2,1) = (-ita*u**2)*0.5+1.0      ; R_inv(2,2) = ita*u               ; R_inv(2,3) = -ita
        R_inv(3,1) = ((ita*u**2)*0.5-u/c)*0.5 ; R_inv(3,2) = -(ita*u-1.d0/c)*0.5 ; R_inv(3,3) = ita*0.5
    end subroutine calcR_inv

    subroutine calcLambda( rho, u, H, gam, Lm )
        implicit none

        real(8)             :: rho , u , e , gam
        real(8)             :: H   , c
        real(8),intent(out) :: Lm(3,3)

        c = (gam-1.0)*(H-0.5*u**2)

        if (c .le.0.0) then
            c = 0.0
        else
            c = dsqrt(c)
        end if
 
        Lm(1:3,1:3) = 0.d0
        Lm(1,1) = dabs(u - c)
        Lm(2,2) = dabs(u)
        Lm(3,3) = dabs(u + c)

    end subroutine calcLambda

    subroutine calcRoeAverage(rL,rR,uL,uR,HL,HR , r2,u2,H2)
        implicit none

        real(8),intent(in)  :: rL , rR
        real(8),intent(in)  :: uL , uR
        real(8),intent(in)  :: HL , HR
        real(8),intent(out) :: r2 , u2 , H2

        r2 = dsqrt(rL*rR)
        u2 = (dsqrt(rL)*uL + dsqrt(rR)*uR)/(dsqrt(rL)+dsqrt(rR))
        H2 = (dsqrt(rL)*HL + dsqrt(rR)*HR)/(dsqrt(rL)+dsqrt(rR))

    end subroutine calcRoeAverage

    function calcH( rho, u, e, gam)
        implicit none

        real(8) :: rho , u , e , gam
        real(8) :: p 
        real(8) :: calcH

        p     = (gam-1.0)*(e-(rho*u**2)*0.5)
        calcH = (e + p)/rho

    end function calcH

    function calcC( rho, u, e, gam )
        implicit none

        real(8) :: rho , u , e , gam
        real(8) :: p   , c2, H
        real(8) :: calcC

        p = (gam-1.0)*(e-(rho*u**2)*0.5)
        H = (e + p)/rho
        c2= (gam-1.0)*(H-0.5*u**2)
        if (c2 .le.0.0) then
            calcC = 0.0
        else
            calcC = dsqrt(c2)
        end if
    
    end function calcC

end module

program main
    use lib4FR
    implicit none

    type Cell
        real(8)         :: center 
        integer         :: nSolPt
    end type Cell

    ! basic parameter
    character*2 :: method
    integer     :: nCell
    integer     :: nSolPt
    real(8)     :: DT
    integer     :: iSTEP , LastSTEP , outSTEP
    real(8)     :: xmax , xmin , dx
    real(8)     :: Cr   , umax

    ! polynomials
    real(8),allocatable    :: P1       (:) , P2        (:) , P3    (:)
    real(8),allocatable    :: P4       (:) , P5        (:) 
    real(8),allocatable    :: P1_dot   (:) , P2_dot    (:) , P3_dot(:)
    real(8),allocatable    :: P4_dot   (:) , P5_dot    (:) 
    real(8),allocatable    :: G_dg3_dot(:) , G_dg4_dot (:) , G_dg5_dot (:)
    real(8),allocatable    :: G_ga4_dot(:)
    real(8),allocatable    :: G_dotL   (:) 
    real(8),allocatable    :: G_dotR   (:) 

    ! Runge-Kutta 
    integer :: iRK
    integer :: nRK ! runge kutta order
    real(8),allocatable :: RK_K(:,:,:,:) !iCell,iSol,k,equation

    character*8  :: cstep
    character*20 :: fname 

    ! temp value
    integer :: i    , iCell , iSol
    integer :: irow 
    integer :: k    , l     , m
    real(8) :: fai_k_dot
    real(8) :: pai_pot , sig_pot
    real(8) :: diffL(3), diffR(3)
    real(8) :: x_in
    real(8) :: val
    real(8) :: fupw   (3)
    real(8) :: lambda (3,3)
    real(8) :: R      (3,3)
    real(8) :: R_inv  (3,3)
    real(8) :: Lm     (3,3)
    real(8) :: A      (3,3)
    real(8) :: uL , uR , uRoe
    real(8) :: rL , rR , rRoe
    real(8) :: HL , HR , HRoe
    real(8) :: eL , eR
    real(8) :: Q_L(3) , Q_R(3)
    real(8) :: E_L(3) , E_R(3)
    real(8) :: v3(3)

    ! mesh and variables
    type(Cell),allocatable :: cells(:)
    real(8),allocatable    :: x (:,:) ! iCell, iSol
    real(8),allocatable    :: xl(:)
    real(8),allocatable    :: D(:,:)

    real(8),allocatable    :: u    (:,:)
    real(8),allocatable    :: uN   (:,:)
    real(8),allocatable    :: rho  (:,:)
    real(8),allocatable    :: rhoN (:,:)
    real(8),allocatable    :: prs  (:,:)
    real(8),allocatable    :: prsN (:,:)
    real(8),allocatable    :: teng (:,:)
    real(8),allocatable    :: tengN(:,:)

    real(8),allocatable    :: f1     (:,:,:)
    real(8),allocatable    :: f1_dot (:,:,:)
    real(8),allocatable    :: f2     (:,:,:)
    real(8),allocatable    :: f2_dot (:,:,:)
    real(8),allocatable    :: f2_dotg(:,:,:)

    real(8),allocatable    :: Q      (:,:,:) , QN(:,:,:)
    real(8),allocatable    :: E      (:,:,:) , E2(:,:,:)
    real(8),allocatable    :: dEdx   (:,:,:) 
    real(8),allocatable    :: dEdx2  (:,:,:) 
    real(8),allocatable    :: dEdx2g (:,:,:) 

    ! shock parameter
    real(8),parameter     :: gam     = 1.4
    real(8),parameter     :: rhoL    = 1.0  
    real(8),parameter     :: prsL    = 1. 
    real(8),parameter     :: velL    = 0.0
    real(8),parameter     :: rhoR    = 0.125
    real(8),parameter     :: prsR    = 0.1
    real(8),parameter     :: velR    = 0.

    ! -----------------------
    ! *** initial setting ***
    ! -----------------------

    !method="DG"
    method="GA"

    if      (method.eq."DG") then ! DG method
        Cr       = 0.09162 ! Courant Number
        !LastSTEP = 764     ! wave travels a distance of 10 periods
        LastSTEP = 764*5   ! wave travels a distance of 50 periods
    else if (method.eq."GA") then
        Cr       = 0.14285 ! Courant Number
        !LastSTEP = 490     ! wave travels a distance of 10 periods
        LastSTEP = 490*5   ! wave travels a distance of 50 periods
    else
        write(*,*) "Error: unknown method"
        stop
    end if

!TEMP>
    !LastSTEP = 1457
    !LastSTEP = 1349
    LastSTEP = 123
!TEMP<

    nCell  = 100
    xmax   = 1.d0 ; xmin = 0.d0
    nSolPt = 4  

    umax = 1.0 ! extimated max velocity

    nRK = 4
    allocate(cells(nCell))

    allocate( x      (nCell,0:nSolPt+1) )
    allocate( u      (nCell,0:nSolPt+1) )
    allocate( uN     (nCell,0:nSolPt+1) )

    allocate( rho    (nCell,0:nSolPt+1) )
    allocate( rhoN   (nCell,0:nSolPt+1) )
    allocate( prs    (nCell,0:nSolPt+1) )
    allocate( prsN   (nCell,0:nSolPt+1) )
    allocate( teng   (nCell,0:nSolPt+1) )
    allocate( tengN  (nCell,0:nSolPt+1) )

    allocate( f1     (nCell,0:nSolPt+1,3) )
    allocate( f1_dot (nCell,0:nSolPt+1,3) )
    allocate( f2     (nCell,0:nSolPt+1,3) )
    allocate( f2_dot (nCell,0:nSolPt+1,3) )
    allocate( f2_dotg(nCell,0:nSolPt+1,3) )

    allocate( Q      (nCell,0:nSolPt+1,3) )
    allocate( QN     (nCell,0:nSolPt+1,3) )
    allocate( E      (nCell,0:nSolPt+1,3) )
    allocate( E2     (nCell,0:nSolPt+1,3) )
    allocate( dEdx   (nCell,0:nSolPt+1,3) )
    allocate( dEdx2  (nCell,0:nSolPt+1,3) )
    allocate( dEdx2g (nCell,0:nSolPt+1,3) )

    allocate( RK_K   (nCell,nSolPt,nRK,3))

    ! -----------------------
    ! *** set gauss point ***
    ! -----------------------
    allocate(xl(0:nSolPt+1))
    if (nSolPt.eq.4) then
        xl(0) = -1.0      ; xl(5) = 1.0
        xl(1) = -0.861136 ; xl(2) = -0.339981
        xl(4) =  0.861136 ; xl(3) =  0.339981
    else 
        write(*,*) "Error: only nSolPt = 2 is available"
        stop
    end if

    ! -------------
    ! *** set x ***
    ! -------------
    dx = (xmax - xmin)/nCell
    DT = CR*dx/umax*0.8! actually, devided by u
    do iCell=1,nCell
        cells(iCell)%nSolPt = nSolPt
        cells(iCell)%center = dx/2 + (iCell-1)*dx
    end do

    do iCell=1,nCell
        x(iCell,0:nSolPt+1) = cells(iCell)%center &
                            + dx/2*xl(0:nSolPt+1)
    end do

    ! -----------------------
    ! *** set polynomials ***
    ! -----------------------
    allocate( P1       (0:nSolPt+1) , P2       (0:nSolPt+1) , P3    (0:nSolPt+1))
    allocate( P4       (0:nSolPt+1) , P5       (0:nSolPt+1) )
    allocate( P1_dot   (0:nSolPt+1) , P2_dot   (0:nSolPt+1) , P3_dot(0:nSolPt+1))
    allocate( P4_dot   (0:nSolPt+1) , P5_dot   (0:nSolPt+1) )
    allocate( G_dg3_dot(0:nSolPt+1) , G_dg4_dot(0:nSolPt+1) , G_dg5_dot(0:nSolPt+1) )
    allocate( G_ga4_dot(0:nSolPt+1) )
    allocate( G_dotL(0:nSolPt+1) )
    allocate( G_dotR(0:nSolPt+1) )

    P4        = (35*xl**4 -30*xl**2 +3)/8
    P2_dot    = 3*xl
    P3_dot    = 0.5*(15*xl**2 -3)
    P4_dot    = 0.5*xl*(35*xl**2 -15)
    P5_dot    = (315*xl**4 -210*xl**2 +15)/8
    G_dg3_dot = (-1)**3*0.5*(P3_dot - P2_dot)
    G_dg4_dot = (-1)**4*0.5*(P4_dot - P3_dot)
    G_dg5_dot = (-1)**5*0.5*(P5_dot - P4_dot)
    G_ga4_dot = 4.d0/(2.d0*3 +1.d0)  *G_dg4_dot &
              + 3.d0/(2.d0*3.d0+1.d0)*G_dg3_dot 

    if      (method.eq."DG") then ! DG method
        G_dotL = G_dg4_dot
    else if (method.eq."GA") then ! GA method
        G_dotL = G_ga4_dot
    else
        write(*,*) "Error: unknown method"
        stop
    end if

    do i=0,nSolPt+1
        G_dotR(i) = - G_dotL(nSolPt+1-i)
    end do

    ! -----------------------------
    ! *** set derivative matrix ***
    ! -----------------------------
    allocate( D(nSolPt,nSolPt) )

    fai_k_dot = 0.d0

    do irow=1,nSolPt
        do k=1,nSolPt
            sig_pot = 0.d0
            do l=1,nSolPt
                if (l.eq.k) cycle
                pai_pot = 1.d0
                do m=1,nSolPt
                    if ( m.eq.l .or. m.eq.k )cycle
                    pai_pot = pai_pot*(xl(irow) - xl(m))/(xl(k) - xl(m))
                end do
                sig_pot = sig_pot + pai_pot/(xl(k)-xl(l))
            end do
            D(irow,k) = sig_pot
        end do
    end do

    ! ------------------------------------------------
    ! *** set initial value for shock tube problem ***
    ! ------------------------------------------------
    open(10,file="initVal.dat",status="replace")
    write(10,*) "# x , rho, u , teng , prs"
    do iCell=1,nCell
        do iSol=1,nSolPt
            if (iCell.le.nCell/2) then
                u   (iCell,iSol)= velL
                rho (iCell,iSol)= rhoL
                prs (iCell,iSol)= prsL
                teng(iCell,iSol)= prsL/(gam-1.0) + rhoL*velL**2*0.5
            else 
                u   (iCell,iSol)= velR
                rho (iCell,iSol)= rhoR
                prs (iCell,iSol)= prsR
                teng(iCell,iSol)= prsR/(gam-1.0) + rhoR*velR**2*0.5
            end if
            Q(iCell,iSol,1) = rho (iCell,iSol)
            Q(iCell,iSol,2) = rho (iCell,iSol)*u(iCell,iSol)
            Q(iCell,iSol,3) = teng(iCell,iSol)
            write(10,*) x(iCell,iSol),rho(iCell,iSol),u(iCell,iSol),teng(iCell,iSol),prs(iCell,iSol)
        end do
    end do
    close(10)


! -------------------------
! *** Start Calculation ***
! -------------------------
    do iSTEP=1,LastSTEP
        do iRK=1,nRK
        ! -----------------------------------------------------------
        ! *** Set Q and E at interface using Lagrange Polynomials ***
        ! -----------------------------------------------------------
            do iCell=1,nCell
                do iSol=1,nSolPt
                    call calcE(Q(iCell,iSol,1:3),gam,E(iCell,iSol,1:3))
                end do
                ! ------------
                ! *** left ***
                ! ------------
                x_in = -1.d0
                ! set Q
                call lagrange_poly(nSolPt, Q(iCell,1:nSolPt,1), xl(1:nSolPt), x_in, val)
                Q(iCell,0,1) = val
                call lagrange_poly(nSolPt, Q(iCell,1:nSolPt,2), xl(1:nSolPt), x_in, val)
                Q(iCell,0,2) = val
                call lagrange_poly(nSolPt, Q(iCell,1:nSolPt,3), xl(1:nSolPt), x_in, val)
                Q(iCell,0,3) = val

                ! set E
                call lagrange_poly(nSolPt, E(iCell,1:nSolPt,1), xl(1:nSolPt), x_in, val)
                E(iCell,0,1) = val
                call lagrange_poly(nSolPt, E(iCell,1:nSolPt,2), xl(1:nSolPt), x_in, val)
                E(iCell,0,2) = val
                call lagrange_poly(nSolPt, E(iCell,1:nSolPt,3), xl(1:nSolPt), x_in, val)
                E(iCell,0,3) = val

                ! -------------
                ! *** right ***
                ! -------------
                x_in =  1.d0
                ! set Q
                call lagrange_poly(nSolPt, Q(iCell,1:nSolPt,1), xl(1:nSolPt), x_in, val)
                Q(iCell,nSolPt+1,1) = val
                call lagrange_poly(nSolPt, Q(iCell,1:nSolPt,2), xl(1:nSolPt), x_in, val)
                Q(iCell,nSolPt+1,2) = val
                call lagrange_poly(nSolPt, Q(iCell,1:nSolPt,3), xl(1:nSolPt), x_in, val)
                Q(iCell,nSolPt+1,3) = val

                ! set E
                call lagrange_poly(nSolPt, E(iCell,1:nSolPt,1), xl(1:nSolPt), x_in, val)
                E(iCell,nSolPt+1,1) = val
                call lagrange_poly(nSolPt, E(iCell,1:nSolPt,2), xl(1:nSolPt), x_in, val)
                E(iCell,nSolPt+1,2) = val
                call lagrange_poly(nSolPt, E(iCell,1:nSolPt,3), xl(1:nSolPt), x_in, val)
                E(iCell,nSolPt+1,3) = val

                if (iRK.eq.1) then ! first step for Runge-Kutta
                    QN(iCell,0:nSolPt+1,1:3) = Q(iCell,0:nSolPt+1,1:3)
                end if
            end do

        ! ---------------------------------------------------
        ! *** set common flux at interface using Roe flux ***
        ! ---------------------------------------------------
            do iCell=1,nCell
                ! ----------------------
                ! *** left interface ***
                ! ----------------------
                if (iCell .eq. 1) then ! left boundary
                    E2(iCell, 0, 1:3) = E(iCell, 1, 1:3)
                else
                    Q_R(1:3) = Q(iCell, 0, 1:3)
                    E_R(1:3) = E(iCell, 0, 1:3)
                    rR       = Q(iCell, 0, 1)
                    uR       = Q(iCell, 0, 2)/rR
                    eR       = Q(iCell, 0, 3)
                    HR       = calcH(rR,uR,eR,gam)

                    Q_L(1:3) = Q(iCell-1, nSolPt+1, 1:3)
                    E_L(1:3) = E(iCell-1, nSolPt+1, 1:3)
                    rL       = Q(iCell-1, nSolPt+1, 1)
                    uL       = Q(iCell-1, nSolPt+1, 2)/rR
                    eL       = Q(iCell-1, nSolPt+1, 3)
                    HL       = calcH(rL,uL,eL,gam)

                    call calcRoeAverage(rL,rR,uL,uR,HL,HR , rRoe,uRoe,HRoe)
                    call calcR         (rRoe, uRoe, HRoe, gam, R )
                    call calcR_inv     (rRoe, uRoe, HRoe, gam, R_inv )
                    call calcLambda    (rRoe, uRoe, HRoe, gam, Lm )

                    A = matmul(R,Lm)
                    A = matmul(A,R_inv)
                    fupw(1:3) = 0.5*(E_L + E_R) -0.5*matmul(A,(Q_R-Q_L))
                    E2(iCell, 0, 1:3) = fupw(1:3)
                end if

                ! right interface
                if (iCell .eq. nCell) then ! right boundary
                    E2(iCell, nSolPt+1, 1:3) = E(iCell,nSolPt,1:3)
                else
                    Q_R(1:3) = Q(iCell+1, 0, 1:3)
                    E_R(1:3) = E(iCell+1, 0, 1:3)
                    rR       = Q(iCell+1, 0, 1)
                    uR       = Q(iCell+1, 0, 2)/rR
                    eR       = Q(iCell+1, 0, 3)
                    HR       = calcH(rR,uR,eR,gam)

                    Q_L(1:3) = Q(iCell, nSolPt+1, 1:3)
                    E_L(1:3) = E(iCell, nSolPt+1, 1:3)
                    rL       = Q(iCell, nSolPt+1, 1)
                    uL       = Q(iCell, nSolPt+1, 2)/rL
                    eL       = Q(iCell, nSolPt+1, 3)
                    HL       = calcH(rL,uL,eL,gam)

                    call calcRoeAverage(rL,rR,uL,uR,HL,HR , rRoe,uRoe,HRoe)
                    call calcR         (rRoe, uRoe, HRoe, gam, R )
                    call calcR_inv     (rRoe, uRoe, HRoe, gam, R_inv )
                    call calcLambda    (rRoe, uRoe, HRoe, gam, Lm )

                    A = matmul(R,Lm)
                    A = matmul(A,R_inv)
                    fupw(1:3) = 0.5*(E_L + E_R) -0.5*matmul(A,(Q_R-Q_L))
                    E2(iCell, nSolPt+1, 1:3) = fupw(1:3)

                end if

            end do

        ! ---------------------------------------------
        ! *** calc derivative of discontinuout flux ***
        ! ---------------------------------------------
            do iCell=1,nCell
                do iSol=1,nSolPt
                    dEdx(iCell,iSol,1) = dot_product(E(iCell,1:nSolPt,1),D(iSol,1:nSolPt))
                    dEdx(iCell,iSol,2) = dot_product(E(iCell,1:nSolPt,2),D(iSol,1:nSolPt))
                    dEdx(iCell,iSol,3) = dot_product(E(iCell,1:nSolPt,3),D(iSol,1:nSolPt))
                end do

            end do

        ! ----------------------------
        ! *** calc continuous flux ***
        ! ----------------------------
            do iCell=1,nCell
                diffL(1:3) = E2(iCell, 0       , 1:3) - E(iCell, 0        ,1:3)
                diffR(1:3) = E2(iCell, nSolPt+1, 1:3) - E(iCell, nSolPt+1 ,1:3)

                dEdx2 (iCell, 1:nSolPt, 1)  = dEdx(iCell,1:nSolPt,1)    &
                                            + G_dotL(1:nSolPt)*diffL(1) &
                                            + G_dotR(1:nSolPt)*diffR(1)


                dEdx2 (iCell, 1:nSolPt, 2)  = dEdx(iCell,1:nSolPt,2)    &
                                            + G_dotL(1:nSolPt)*diffL(2) &
                                            + G_dotR(1:nSolPt)*diffR(2)

                dEdx2 (iCell, 1:nSolPt, 3)  = dEdx(iCell,1:nSolPt,3)    &
                                            + G_dotL(1:nSolPt)*diffL(3) &
                                            + G_dotR(1:nSolPt)*diffR(3)

                dEdx2g(iCell, 1:nSolPt, 1:3)  = 2.d0/dx*dEdx2(iCell, 1:nSolPt, 1:3)


            end do

        ! -----------------------------
        ! *** time advance          *** 
        ! *** 4th order Runge-Kutta ***
        ! -----------------------------
            do iCell=1,nCell
                RK_K(iCell, 1:nSolPt, iRK, 1:3) = -dEdx2g(iCell, 1:nSolPt, 1:3)

                if      (iRK.eq.1) then
                    Q(iCell, 1:nSolPt, 1:3) = QN(iCell, 1:nSolPt, 1:3) &
                                            + DT*RK_K(iCell, 1:nSolPt, iRK, 1:3)*0.5
                else if (iRK.eq.2) then
                    Q(iCell, 1:nSolPt, 1:3) = QN(iCell, 1:nSolPt, 1:3) &
                                            + DT*RK_K (iCell, 1:nSolPt, iRK, 1:3)*0.5
                else if (iRK.eq.3) then
                    Q(iCell, 1:nSolPt, 1:3) = QN(iCell, 1:nSolPt, 1:3) &
                                            + DT*RK_K (iCell, 1:nSolPt, iRK, 1:3)
                else if (iRK.eq.4) then
                    Q(iCell, 1:nSolPt, 1:3) = QN(iCell, 1:nSolPt, 1:3)            &
                                            + DT*(   RK_K(iCell,1:nSolPt,1,1:3) &
                                                  +2*RK_K(iCell,1:nSolPt,2,1:3) &
                                                  +2*RK_K(iCell,1:nSolPt,3,1:3) &
                                                  +  RK_K(iCell,1:nSolPt,4,1:3) )/6
                else 
                    write(*,*) "Error: something wrong"
                    stop
                end if

            end do
        end do

    end do

    ! ------------------------
    ! *** output end value ***
    ! ------------------------
    open(10,file="endVal.dat",status="replace")
    write(10,*) "# x , rho, u , teng , prs"
    do iCell=1,nCell
        do iSol=1,nSolPt
            rho (iCell,iSol) = Q(iCell,iSol,1)
            u   (iCell,iSol) = Q(iCell,iSol,2)/rho(iCell,iSol)
            teng(iCell,iSol) = Q(iCell,iSol,3)
            prs (iCell,iSol) = (gam-1.0)*(teng(iCell,iSol)-0.5*rho(iCell,iSol)*u(iCell,iSol)**2)

            write(10,*) x(iCell,iSol),rho(iCell,iSol),u(iCell,iSol),teng(iCell,iSol),prs(iCell,iSol)
        end do
    end do
    close(10)



end program main

