   	program SWE_TANK

   	implicit none

	! Numerical solution to the Shallow Water Equations (Falconer, 1993)

	! Gabriel Thomas Scarlett, University of Edinburgh , U.K.

   	! Solves the Shallow Water Equations equatios in 1D space using:
	! Finite difference scheme
   	! Runge-Kutta fourth order time integration
	
	! This program addresses the problem of waves sloshing in a tank.
	! (see H. Lamb, 1932)

    ! WAVE I (h=2)

   	parameter(imax = 2001)
	double precision :: cf,S0,dx,g,dt,rho,tend,t,amp,Tp,om,kn,finish
	double precision :: hs(imax),h(imax),q(imax),zeta(imax)
	double precision :: znew(imax),qnew(imax),taubx
	integer :: imax

	call input (imax,S0,dx,cf,g,dt,rho,tend,zeta,hs,h,q,amp,Tp,om,kn)

     t=0.0d0			! time

500 continue

  	t=t+dt

    call calc (imax,hs,q,zeta,h,dx,cf,g,dt,rho,znew,qnew,taubx)

    call update (imax,q,zeta,h,hs,znew,qnew)

    open(4,file='output.dat',status='old')

  	write(4,*)t/Tp,zeta(2)/amp


	!100 format(f8.2,2(2x,f15.12))
     
  	if (t.lt.tend) goto 500	

		call cpu_time(finish)

        print '("Time = ",f6.3," seconds.")',finish

	stop

	end program SWE_TANK

!***************************************************************

	subroutine input (imax,S0,dx,cf,g,dt,rho,tend,zeta,hs,h,q,amp,Tp,om,kn)
    implicit none
	double precision :: pi,lx,g,S0,amp,h0,dx,x,cf,dt,rho,tend,Tp,kn,om

	double precision :: h(imax),hs(imax),q(imax),zeta(imax)
	integer i,imax
	
	!
    !   constants
	!

    pi=4.0d0*datan(1.0d0)

    !
	!  physical variables
	!

	g=9.81d0		! gravitational acceleration
    rho=1000.0d0	! density
	lx=1000.0d0		! domain length
	S0=0.0d0		! bed gradient
	cf=0.0d0		! bed friction coefficient
    h0=2.0d0
    amp=h0/100.0d0

    !
    !  numerical parameters
	!
     
	kn=2d0*pi/lx		! wave number
    om=kn*sqrt(g*h0)	! frequency
    Tp=2d0*pi/om		! period
	dx=lx/dble(imax-1)	! space step size
    dt=0.2d0			! time step
    tend=5.0d0*Tp		! duration of simulation
	
    !
    !   initial values of free surface elevation, still water depth
    !   total depth and flux
    !
	
	do i=1,imax
 	x=dble(i-1)*dx 					! distance along bed
  	zeta(i)=amp*cos(2.0d0*pi*x/lx) 	! initial free surface elevation
  	hs(i)=h0  						! still water level
  	h(i)=hs(i)+zeta(i)				! water depth		
  	q(i)=0.0d0	   					! local flux	
  	end do

  	end subroutine input

!***************************************************************
    subroutine calc (imax,hs,q,zeta,h,dx,cf,g,dt,rho,znew,qnew,taubx)
    implicit none
    double precision :: g,dt,dx,cf,u,rho,taubx,dzdt4,dqdt4,z1,ze,zw,h1,he,hw,q1,qe,qw
    double precision :: q(imax),zeta(imax),h(imax),hs(imax)
    double precision :: dzdt1(imax),dzdt2(imax),dzdt3(imax)
    double precision :: dqdt1(imax),dqdt2(imax),dqdt3(imax)
	double precision :: znew(imax),qnew(imax)
	integer i,imax

    
	!
	!  first step in RK4
	!
    
	do i=2,imax-1
	dzdt1(i)=-(q(i+1)-q(i-1))/2.0/dx
  	u=q(i)/h(i)						
  	taubx=cf*rho*u*dabs(u)			! bed shear stress
  	dqdt1(i)=-(q(i+1)*q(i+1)/h(i+1)-q(i-1)*q(i-1)/h(i-1))/2.0d0/dx &	
           -g*h(i)*(zeta(i+1)-zeta(i-1))/2.0d0/dx  &
              -taubx/rho
       end do
       dzdt1(1)=2.0d0*dzdt1(2)-dzdt1(3)
       dqdt1(1)=0.0d0
       dzdt1(imax)=2.0d0*dzdt1(imax-1)-dzdt1(imax-2)
       dqdt1(imax)=0.0d0
	!
	!  second step in RK4
	!
       do i=2,imax-1
       z1=zeta(i)+dt/2.0d0*dzdt1(i)
       ze=zeta(i+1)+dt/2.0d0*dzdt1(i+1)
       zw=zeta(i-1)+dt/2.0d0*dzdt1(i-1)
       h1=hs(i)+z1
       he=hs(i+1)+ze
       hw=hs(i-1)+zw
       q1=q(i)+dt/2.0d0*dqdt1(i)
       qe=q(i+1)+dt/2.0d0*dqdt1(i+1)
       qw=q(i-1)+dt/2.0d0*dqdt1(i-1)
       dzdt2(i)=-(qe-qw)/2.0/dx
       u=q1/h1						
       taubx=cf*rho*u*dabs(u)			! bed shear stress
       dqdt2(i)=-(qe*qe/he-qw*qw/hw)/2.0d0/dx &	
           -g*h1*(ze-zw)/2.0d0/dx  &
              -taubx/rho
       end do
       dzdt2(1)=2.0d0*dzdt2(2)-dzdt2(3)
       dqdt2(1)=0.0d0
       dzdt2(imax)=2.0d0*dzdt2(imax-1)-dzdt2(imax-2)
       dqdt2(imax)=0.0d0
   	!
	!  third step in RK4
	!
       do i=2,imax-1
       z1=zeta(i)+dt/2.0d0*dzdt2(i)
       ze=zeta(i+1)+dt/2.0d0*dzdt2(i+1)
       zw=zeta(i-1)+dt/2.0d0*dzdt2(i-1)
       h1=hs(i)+z1
       he=hs(i+1)+ze
       hw=hs(i-1)+zw
       q1=q(i)+dt/2.0d0*dqdt2(i)
       qe=q(i+1)+dt/2.0d0*dqdt2(i+1)
       qw=q(i-1)+dt/2.0d0*dqdt2(i-1)
       dzdt3(i)=-(qe-qw)/2.0/dx
       u=q1/h1						
       taubx=cf*rho*u*dabs(u)			! bed shear stress
       dqdt3(i)=-(qe*qe/he-qw*qw/hw)/2.0d0/dx &	
           -g*h1*(ze-zw)/2.0d0/dx  &
              -taubx/rho
       end do
       dzdt3(1)=2.0d0*dzdt3(2)-dzdt3(3)
       dqdt3(1)=0.0d0
       dzdt3(imax)=2.0d0*dzdt3(imax-1)-dzdt3(imax-2)
       dqdt3(imax)=0.0d0
  	!
	!  fourth step in RK4
	!
       do i=2,imax-1
       z1=zeta(i)+dt*dzdt3(i)
       ze=zeta(i+1)+dt*dzdt3(i+1)
       zw=zeta(i-1)+dt*dzdt3(i-1)
       h1=hs(i)+z1
       he=hs(i+1)+ze
       hw=hs(i-1)+zw
       q1=q(i)+dt*dqdt3(i)
       qe=q(i+1)+dt*dqdt3(i+1)
       qw=q(i-1)+dt*dqdt3(i-1)
       dzdt4=-(qe-qw)/2.0/dx
       u=q1/h1						
       taubx=cf*rho*u*dabs(u)			! bed shear stress
       dqdt4=-(qe*qe/he-qw*qw/hw)/2.0d0/dx &	
           -g*h1*(ze-zw)/2.0d0/dx  &
              -taubx/rho
       znew(i)=zeta(i)+dt*(dzdt1(i)+2.0d0*dzdt2(i) &
                      +2.0d0*dzdt3(i)+dzdt4)/6.0d0 
       qnew(i)=q(i)+dt*(dqdt1(i)+2.0d0*dqdt2(i) &
                      +2.0d0*dqdt3(i)+dqdt4)/6.0d0      
  	end do
    !
    ! boundary conditions
    !
  	znew(1)=2.0d0*znew(2)-znew(3)
  	qnew(1)=0.0d0
  	znew(imax)=2.0d0*znew(imax-1)-znew(imax-2)
  	qnew(imax)=0.0d0


    end subroutine calc

!******************************************************************

    subroutine update (imax,q,zeta,h,hs,znew,qnew)
    implicit none
    double precision :: q(imax),zeta(imax),h(imax),hs(imax)
	double precision :: znew(imax),qnew(imax)
	integer i,imax
  	do i=1,imax
    zeta(i)=znew(i)
    h(i)=hs(i)+zeta(i)
    q(i)=qnew(i)
  	end do

    end subroutine update

!******************************************************************


	subroutine out (imax,q,h)
    implicit none
    double precision :: q(imax),h(imax)
    integer :: imax
    write(*,100)h,q/h					! prints the horizontal velocity
    100 format(f8.2,2(2x,f15.12))
    
	end subroutine out