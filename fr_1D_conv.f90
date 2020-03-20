! 1D convection using Flux Reconstruction Method 
! coded by Kumpei Sano. 2020 03 20

program main
    implicit none

    type Cell
        real(8)         :: center 
        integer         :: nSolPt
    end type Cell

    ! parameter
    integer :: nCell
    integer :: nSolPt
    real(8) :: DT
    integer :: iSTEP , LastSTEP , outSTEP
    real(8) :: a 
    real(8) :: xmax , xmin , dx
    real(8) :: Cr

    ! mesh and variables
    type(Cell),allocatable :: cells(:)
    real(8),allocatable    :: x (:,:) ! iCell, iSol
    real(8),allocatable    :: xl(:)
    real(8),allocatable    :: D(:,:)

    real(8),allocatable    :: u      (:,:)
    real(8),allocatable    :: uN     (:,:)
    real(8),allocatable    :: f1     (:,:)
    real(8),allocatable    :: f1_dot (:,:)
    real(8),allocatable    :: f2     (:,:)
    real(8),allocatable    :: f2_dot (:,:)
    real(8),allocatable    :: f2_dotg(:,:)

    ! polynomials
    real(8),allocatable    :: P1    (:) , P2    (:) , P3    (:)
    real(8),allocatable    :: P4    (:) , P5    (:) 
    real(8),allocatable    :: P1_dot(:) , P2_dot(:) , P3_dot(:)
    real(8),allocatable    :: P4_dot(:) , P5_dot(:) 
    real(8),allocatable    :: G_dg3_dot(:) , G_dg4_dot (:) , G_dg5_dot (:)
    real(8),allocatable    :: G_ga4_dot (:)
    real(8),allocatable    :: G_dotL(:) 
    real(8),allocatable    :: G_dotR(:) 

    ! Runge-Kutta 
    integer :: iRK
    integer :: nRK ! runge kutta order
    real(8),allocatable :: RK_K(:,:,:) !iCell,iSol,k

    character*8  :: cstep
    character*20 :: fname 

    ! temp value
    integer :: i    , iCell , iSol
    integer :: irow 
    integer :: k    , l     , m
    real(8) :: fai_k_dot
    real(8) :: pai_pot , sig_pot
    real(8) :: diffL, diffR
    real(8),allocatable :: u_temp(:) , x_temp(:) , f_temp(:)
    real(8) :: x_in
    real(8) :: val
    real(8) :: fupw

    ! -----------------------
    ! *** initial setting ***
    ! -----------------------
    nCell = 10
    xmax  = 1.d0 ; xmin = 0.d0
    nSolPt = 4  
    !LastSTEP = 764*5
    LastSTEP = 764 !Gdg
    !LastSTEP = 490 ! Gga
    !LastSTEP = 490*5
    !LastSTEP = 2
    !LastSTEP = 10   
    !LastSTEP = 2000
    outSTEP = 10
    a = 1.0
    Cr= 0.09162 ! Gdg
    !Cr= 0.14285 ! Gga
    nRK = 4
    allocate(cells(nCell))

    allocate(x      (nCell,0:nSolPt+1))
    allocate(u      (nCell,0:nSolPt+1))
    allocate(uN     (nCell,0:nSolPt+1))
    allocate(f1     (nCell,0:nSolPt+1))
    allocate(f1_dot (nCell,0:nSolPt+1))
    allocate(f2     (nCell,0:nSolPt+1))
    allocate(f2_dot (nCell,0:nSolPt+1))
    allocate(f2_dotg(nCell,0:nSolPt+1))

    allocate(RK_K(nCell,nSolPt,nRK))

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
    DT = CR*dx
    do iCell=1,nCell
        cells(iCell)%nSolPt    = nSolPt
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

    P4       (0:nSolPt+1) = (35*xl(0:nSolPt+1)**4 -30*xl(0:nSolPt+1)**2 +3)/8
    P3_dot   (0:nSolPt+1) = 0.5*(15*xl(0:nSolPt-1)**2 -3)
    P4_dot   (0:nSolPt+1) = 0.5*xl(0:nSolPt-1)*(35*xl(0:nSolPt-1)**2 -15)
    P5_dot   (0:nSolPt+1) = (315*xl(0:nSolPt-1)**4 -210*xl(0:nSolPt-1)**2 +15)/8
    G_dg3_dot(0:nSolPt+1) = (-1)**3*0.5*(P3_dot(0:nSolPt-1) - P2_dot(0:nSolPt-1))
    G_dg4_dot(0:nSolPt+1) = (-1)**4*0.5*(P4_dot(0:nSolPt-1) - P3_dot(0:nSolPt-1))
    G_dg5_dot(0:nSolPt+1) = (-1)**5*0.5*(P5_dot(0:nSolPt-1) - P4_dot(0:nSolPt-1))
    G_ga4_dot(0:nSolPt+1) = 4.d0/(2.d0*3 +1.d0)  *G_dg4_dot(0:nSolPt+1) &
                          + 3.d0/(2.d0*3.d0+1.d0)*G_dg3_dot(0:nSolPt+1) 

    G_dotL = G_dg4_dot
    !G_dotL = G_ga4_dot
    do i=0,nSolPt+1
        G_dotR(i) = - G_dotL(nSolPt+1-i)
    end do

    open(10,file="G_dot",status="replace")
    do i=0,nSolPt+1
        write(10,*) xl(i),G_dotL(i),G_dotR(i)
    end do
    close(10)

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

    ! set initial value
    open(10,file="initVal.dat",status="replace")
    do iCell=1,nCell
        do iSol=1,nSolPt
            !i = (iCell-1)*nSolPt + iSol
            u(iCell,iSol) = exp( -40.0*(x(iCell,iSol)-0.5)**2 )
            write(10,*) x(iCell,iSol),u(iCell,iSol)
        end do
    end do
    close(10)


    allocate( u_temp(1:nSolPt),x_temp(1:nSolPt),f_temp(1:nSolPt) )
! -------------------------
! *** Start Calculation ***
! -------------------------
    do iSTEP=1,LastSTEP
        do iRK=1,nRK
        ! set u(j+1/2,L) and u(j+1/2,R)
        !     f(j+1/2,L) and f(j+1/2,R)
            do iCell=1,nCell
                x_temp  (1:nSolPt) = xl (1:nSolPt)
                u_temp  (1:nSolPt) = u  (iCell,1:nSolPt)
                f_temp  (1:nSolPt) = a*u(iCell,1:nSolPt)
                f1(iCell,1:nSolPt) = f_temp(1:nSolPt)
                x_in = -1.d0
                call lagrange_poly(nSolPt, u_temp, x_temp, x_in, val)
                u (iCell,0) = val
                call lagrange_poly(nSolPt, f_temp, x_temp, x_in, val)
                f1(iCell,0) = val
                x_in =  1.d0
                call lagrange_poly(nSolPt, u_temp, x_temp, x_in, val)
                u (iCell,nSolPt+1) = val
                call lagrange_poly(nSolPt, f_temp, x_temp, x_in, val)
                f1(iCell,nSolPt+1) = val

                if (iRK.eq.1) uN(iCell,0:nSolPt+1) = u(iCell,0:nSolPt+1)
            end do

        ! set upwind flux at interface
            do iCell=1,nCell
                if (iCell .eq. 1) then
                    fupw               = a*u(nCell,nSolPt+1)
                    f2(iCell,0)        = fupw
                    fupw               = a*u(iCell,nSolPt+1)
                    f2(iCell,nSolPt+1) = fupw
                else
                    fupw               = a*u(iCell-1,nSolPt+1)
                    f2(iCell,0)        = fupw
                    fupw               = a*u(iCell  ,nSolPt+1)
                    f2(iCell,nSolPt+1) = fupw
                end if
            end do

        ! calc derivative of discontinuout flux at the nSolPt solution points
            do iCell=1,nCell
                do iSol=1,nSolPt
                    f1_dot(iCell,iSol) = dot_product(f1(iCell,1:nSolPt),D(iSol,1:nSolPt))
                end do
            end do

        ! calc continuous flux
            do iCell=1,nCell
                diffL = f2(iCell,0)        - f1(iCell,0)
                diffR = f2(iCell,nSolPt+1) - f1(iCell,nSolPt+1)
                f2_dot (iCell,1:nSolPt)  = f1_dot(iCell,1:nSolPt) &
                                         + diffL*G_dotL(1:nSolPt) &
                                         + diffR*G_dotR(1:nSolPt)
                f2_dotg(iCell,1:nSolPt)  = 2.d0/dx*f2_dot(iCell,1:nSolPt)
            end do

        ! time advance
        ! 4th order Runge-Kutta
            do iCell=1,nCell
                RK_K(iCell,1:nSolPt,iRK) = -f2_dotg(iCell,1:NSolPt)

                if      (iRK.eq.1) then
                    u(iCell,1:nSolPt) = uN(iCell,1:nSolPt) &
                                      + DT*RK_K (iCell,1:nSolPt,iRK)/2
                else if (iRK.eq.2) then
                    u(iCell,1:nSolPt) = uN(iCell,1:nSolPt) &
                                      + DT*RK_K (iCell,1:nSolPt,iRK)/2
                else if (iRK.eq.3) then
                    u(iCell,1:nSolPt) = uN(iCell,1:nSolPt) &
                                      + DT*RK_K (iCell,1:nSolPt,iRK)
                else if (iRK.eq.4) then
                    u(iCell,1:nSolPt) = uN(iCell,1:nSolPt)            &
                                      + DT*(   RK_K(iCell,1:nSolPt,1) &
                                            +2*RK_K(iCell,1:nSolPt,2) &
                                            +2*RK_K(iCell,1:nSolPt,3) &
                                            +  RK_K(iCell,1:nSolPt,4) )/6
                else 
                    write(*,*) "Error: something wrong"
                    stop
                end if
            end do
        end do

    end do

    ! end value
    open(10,file="endVal.dat",status="replace")
    do iCell=1,nCell
        do iSol=1,nSolPt
            write(10,*) x(iCell,iSol),u(iCell,iSol)
        end do
    end do
    close(10)


end program main

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
