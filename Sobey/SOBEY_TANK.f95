   	program SOBEY_TANK

   	implicit none
	
	! Numerical solution to phase-resolving depth averaged wave evolution equations (R. J. Sobey, 2014)

	! Gabriel Thomas Scarlett, University of Edinburgh , U.K.

   	! Solves Sobey's equatios in 1D space using:
	! Finite difference scheme
    ! Tri-diagonal matrix algebra with Thomas algorithm
   	! Runge-Kutta fourth order time integration
	
	! This program addresses the problem of waves sloshing in a tank.
	! (see H. Lamb, 1932)

    ! WAVE I (h=2)

   	parameter(imax = 2001)
    double precision :: Ia,Ip10,Ip11,Ip2,tp,amp
	double precision :: cf,dx,g,dt,rho,tend,t,finish
	double precision :: hs(imax),h(imax),q(-1:imax+2),zeta(imax)
	double precision :: znew(imax),qnew(-1:imax+2),taubx,x(imax)
	integer :: imax,i

	call input (imax,dx,cf,g,dt,rho,tend,zeta,hs,h,q,x,Tp,amp)
    call shape(Ia,Ip10,Ip11,Ip2)

    t=0.0d0			! time

500 continue

  	t=t+dt

    call calc (imax,hs,q,zeta,h,dx,cf,g,dt,rho,znew,qnew,taubx,Ia,Ip10,Ip11,Ip2)
    call update (imax,q,zeta,h,hs,znew,qnew)

        open(4,file='tank_W2.dat',status='old')
	write(4,*)t/Tp,zeta(2)/amp
    write(*,*)t,zeta(500)

     
  	if (t.lt.tend) goto 500
    	
      
	call cpu_time(finish)

    print '("Time = ",f6.3," seconds.")',finish
    close(4)

	stop

	end program SOBEY_TANK

!***************************************************************

	subroutine input (imax,dx,cf,g,dt,rho,tend,zeta,hs,h,q,x,tp,amp)
    implicit none
    double precision :: kn,amp,om,Tp,x
	double precision :: pi,lx,g,h0,dx,cf,dt,rho,tend
	double precision :: h(imax),hs(imax),q(-1:imax+2),zeta(imax)
	integer i,imax

    !   constants

    pi=4.0d0*datan(1.0d0)
	
	!
    !  physical variables
	!
	
	g=9.81d0						! gravitational acceleration
    rho=1000.0d0					! density
	lx=1000.0d0						! domain length
	cf=0.0d0						! bed friction coefficient
    h0=2.0d0						! initial water depth
    amp=h0/100.0d0  				! initial disturbance

    !
    !  numerical parameters
    !
	
	kn=2.0d0*pi/lx					! wave number
    om=kn*sqrt(g*h0)				! frequency
    Tp=2.0d0*pi/om					! period
	dx=lx/dble(imax-1)				! space step size
    dt=0.02d0						! time step
    tend=5.0d0*Tp					! duration of simulation

    !
    !   initial values of free surface elevation, still water depth
    !   total depth and flux
    !
	 
	do i=1,imax
 	x=dble(i-1)*dx 					! distance along bed
  	zeta(i)=amp*cos(2.0d0*pi*x/lx) 	! free surface elevation
  	hs(i)=h0  						! still water level
  	h(i)=hs(i)+zeta(i)				! water depth (still water level + free surface elevation)	
  	q(i)=0.0d0	   	    			! local flux	
  	end do
    
	!
	!   Boundary conditions
	!
     
        q(0)=0.0d0
        q(-1)=0.0d0
        q(imax+1)=0.0d0
        q(imax+2)=0.0d0
        
  	end subroutine input

!***************************************************************
	
	subroutine shape(Ia,Ip10,Ip11,Ip2)
    implicit none
    double precision :: Ia,Ip10,Ip11,Ip2


    ! MODULE HERE TO ACCESS pms1d.95 & splineIa,Ip10,Ip11,I2.f95
    ! to determine shape factors

    !For now assume wave(I) (Sobey,p.6)

	!$$$$$$     Ia =1.024d0
	!$$$$$$     Ip10=1.045d0
	!$$$$$$     Ip11=-0.563d0
	!$$$$$$     Ip2=0.910d0


    Ia =1.0d0
    Ip10=1.045d0
    Ip11=0.0d0
    Ip2=0.0d0
    end subroutine shape
    

!***************************************************************

    subroutine calc (imax,hs,q,zeta,h,dx,cf,g,dt,rho,znew,qnew,taubx,Ia,Ip10,Ip11,Ip2)

    implicit none
    double precision :: g,dt,dx,cf,u,rho,taubx,z1,ze,zw,h1,he,hw,q1,qe,qw
    double precision :: u1,qee,qww
    double precision :: factor,Ia,Ip10,Ip11,Ip2, part1,part2,part3,t1,t2,t3
    double precision :: q(-1:imax+2),zeta(imax),h(imax),hs(imax)
    double precision :: F(imax),a(imax),b(imax),c(imax)
    double precision :: dqdt1(-1:imax+2),dqdt2(-1:imax+2),dqdt3(-1:imax+2),dqdt4(-1:imax+2)
    double precision :: dzdt1(imax),dzdt2(imax),dzdt3(imax),dzdt4(imax)
	double precision :: znew(imax),qnew(-1:imax+2)
	integer i,imax
	!
	!  first step in RK4
	!


	do i=2, imax-1
	
	!
	! continuity equation
	!
	
    dzdt1(i)=-(q(i+1)-q(i-1))/2.0d0/dx 
	
	!
	! momentum equation
	!
	
     u=q(i)/h(i)						
     taubx=cf*rho*u*dabs(u)

     part1=-g*h(i)*(zeta(i+1)-zeta(i-1))/2.0d0/dx
     part2=-Ia*(q(i+1)*q(i+1)/h(i+1)-q(i-1)*q(i-1)/h(i-1))/2.0d0/dx
     
     t1=(q(i+1)-q(i-1))/2.0d0/dx*(q(i+1)/h(i+1)-2.0d0*q(i)/h(i)+q(i-1)/h(i-1))/dx/dx
     t2=2.0d0*(q(i+1)/h(i+1)-q(i-1)/h(i-1))/2.0d0/dx*(q(i+1)-2.0d0*q(i)+q(i-1))/dx/dx
     t3=q(i)/h(i)*(-q(i-2)+2.0d0*q(i-1)-2.0d0*q(i+1)+q(i+2))/2.0d0/dx/dx/dx
     
     part3=Ip2*hs(i)*hs(i)/2.0d0*(t1+t2+t3)
     F(i)=part1+part2+part3-taubx/rho       
        
    a(i) = -(Ip10*hs(i)*hs(i))/2.0d0 +Ip11*hs(i)*(zeta(i))/dx/dx
    b(i) = 1.0d0+2.0d0*(Ip10*hs(i)*hs(i))/2.0d0 +Ip11*hs(i)*(zeta(i))/dx/dx
    c(i) = -(Ip10*hs(i)*hs(i))/2.0d0 +Ip11*hs(i)*(zeta(i))/dx/dx
	end do
     
	do i=3,imax-1
	factor=a(i)/b(i-1)
	b(i)=b(i)-factor*c(i-1)
	F(i)=F(i)-factor*F(i-1)
	end do
	
	!
    ! Backward substitution
	!
	
	dqdt1(imax-1)=F(imax-1)/b(imax-1)
	dqdt1(imax-2)=(F(imax-2)-c(imax-2)*dqdt1(imax-1))/b(imax-2)
    
	do i=imax-2,2,-1
  	dqdt1(i)=(F(i)-c(i)*dqdt1(i+1))/b(i)
	end do

    !
    ! boundary conditions
    !
	
       dzdt1(1)=4.0d0*dzdt1(2)-6.0d0*dzdt1(3)+4.0d0*dzdt1(4)-dzdt1(5)
       dzdt1(imax)=4.0d0*dzdt1(imax-1)-6.0d0*dzdt1(imax-2)+4.0d0*dzdt1(imax-3)-dzdt1(imax-4)
       
       dqdt1(1)=0.0d0
       dqdt1(0)=0.0d0
       dqdt1(-1)=0.0d0
       dqdt1(imax)=0.0d0
       dqdt1(imax+1)=0.0d0
       dqdt1(imax+2)=0.0d0
	   
	!
	!  second step in RK4
	!

        do i=2,imax-1

        z1=zeta(i)+0.5d0*dt*dzdt1(i)
        ze=zeta(i+1)+0.5d0*dt*dzdt1(i+1)
        zw=zeta(i-1)+0.5d0*dt*dzdt1(i-1)

        h1=hs(i)+z1
        he=hs(i+1)+ze
        hw=hs(i-1)+zw
        
        q1=q(i)+0.5d0*dt*dqdt1(i)
        qe=q(i+1)+0.5d0*dt*dqdt1(i+1)
        qw=q(i-1)+0.5d0*dt*dqdt1(i-1)
        qee=q(i+2)+0.5d0*dt*dqdt1(i+2)
        qww=q(i-2)+0.5d0*dt*dqdt1(i-2)

	!
	! continuity equation
	!
	
        dzdt2(i)=-(qe-qw)/2.0d0/dx 
		
	!
	! momentum equation
	!
	
     u1=q1/h1     						
     taubx=cf*rho*u1*dabs(u1)
     part1=-g*h1*(ze-zw)/2.0d0/dx
     part2=-Ia*(qe*qe/he-qw*qw/hw)/2.0d0/dx
     
	!$$$$$$      t1=(ze-zw)/2.0d0/dx*(qe/he-2.0d0*q1/h1+qw/hw)/dx/dx
	!$$$$$$      t2=(qe/he-qw/hw)/2.0d0/dx*(ze-2.0d0*z1+zw)/dx/dx
	!$$$$$$      t3=q1/h1*(-zww+2.0d0*zw-2.0d0*ze+zee)/2.0d0/dx/dx/dx

     t1=(qe-qw)/2.0d0/dx*(qe/he-2.0d0*q1/h1+qw/hw)/dx/dx
     t2=2.0d0*(qe/he-qw/hw)/2.0d0/dx*(qe-2.0d0*q1+qw)/dx/dx
     t3=q1/h1*(-qww+2.0d0*qw-2.0d0*qe+qee)/2.0d0/dx/dx/dx
     
     part3=Ip2*hs(i)*hs(i)/2.0d0*(t1+t2+t3)
     F(i)=part1+part2+part3-taubx/rho       
  
	!$$$$$$     a(i) = -(Ip10*hs(i)*hs(i))/2.0d0 +Ip11*h1*z1/dx/dx
	!$$$$$$     b(i) = 1.0d0+2.0d0*(Ip10*h1*h1)/2.0d0 +Ip11*h1*z1/dx/dx
	!$$$$$$     c(i) = -(Ip10*h1*h1)/2.0d0 +Ip11*h1*z1/dx/dx

    a(i) = -(Ip10*hs(i)*hs(i))/2.0d0 +Ip11*hs(i)*z1/dx/dx
    b(i) = 1.0d0+2.0d0*(Ip10*hs(i)*hs(i))/2.0d0 +Ip11*hs(i)*z1/dx/dx
    c(i) = -(Ip10*hs(i)*hs(i))/2.0d0 +Ip11*hs(i)*z1/dx/dx
	end do
 
  
	do i=3,imax-1
	factor=a(i)/b(i-1)
	b(i)=b(i)-factor*c(i-1)
	F(i)=F(i)-factor*F(i-1)
	end do
	
	!
    ! Backward substitution
	!
	
	dqdt2(imax-1)=F(imax-1)/b(imax-1)
	dqdt2(imax-2)=(F(imax-2)-c(imax-2)*dqdt2(imax-1))/b(imax-2)
    
	do i=imax-2,2,-1
  	dqdt2(i)=(F(i)-c(i)*dqdt2(i+1))/b(i)
	end do
	
    !
    ! boundary conditions
    !


       dzdt2(1)=4.0d0*dzdt2(2)-6.0d0*dzdt2(3)+4.0d0*dzdt2(4)-dzdt2(5)
       dzdt2(imax)=4.0d0*dzdt2(imax-1)-6.0d0*dzdt2(imax-2)+4.0d0*dzdt2(imax-3)-dzdt2(imax-4)
       
       dqdt2(1)=0.0d0
       dqdt2(0)=0.0d0
       dqdt2(-1)=0.0d0
       dqdt2(imax)=0.0d0
       dqdt2(imax+1)=0.0d0
       dqdt2(imax+2)=0.0d0

	!
	!  third step in RK4
	!

        do i=2,imax-1

        z1=zeta(i)+0.5d0*dt*dzdt2(i)
        ze=zeta(i+1)+0.5d0*dt*dzdt2(i+1)
        zw=zeta(i-1)+0.5d0*dt*dzdt2(i-1)

        h1=hs(i)+z1
        he=hs(i+1)+ze
        hw=hs(i-1)+zw
        
        q1=q(i)+0.5d0*dt*dqdt2(i)
        qe=q(i+1)+0.5d0*dt*dqdt2(i+1)
        qw=q(i-1)+0.5d0*dt*dqdt2(i-1)
        qee=q(i+2)+0.5d0*dt*dqdt2(i+2)
        qww=q(i-2)+0.5d0*dt*dqdt2(i-2)
		
	!
	! continuity equation
	!
	
        dzdt3(i)=-(qe-qw)/2.0d0/dx 
		
	!
	! momentum equation
	!
	
     u1=q1/h1     						
     taubx=cf*rho*u1*dabs(u1)
    
     part1=-g*h1*(ze-zw)/2.0d0/dx
     part2=-Ia*(qe*qe/he-qw*qw/hw)/2.0d0/dx
     t1=(qe-qw)/2.0d0/dx*(qe/he-2.0d0*q1/h1+qw/hw)/dx/dx
     t2=2.0d0*(qe/he-qw/hw)/2.0d0/dx*(qe-2.0d0*q1+qw)/dx/dx
     t3=q1/h1*(-qww+2.0d0*qw-2.0d0*qe+qee)/2.0d0/dx/dx/dx
     part3=Ip2*hs(i)*hs(i)/2.0d0*(t1+t2+t3)
     F(i)=part1+part2+part3-taubx/rho       

    a(i) = -(Ip10*hs(i)*hs(i))/2.0d0 +Ip11*hs(i)*z1/dx/dx
    b(i) = 1.0d0+2.0d0*(Ip10*hs(i)*hs(i))/2.0d0 +Ip11*hs(i)*z1/dx/dx
    c(i) = -(Ip10*hs(i)*hs(i))/2.0d0 +Ip11*hs(i)*z1/dx/dx
	end do
     
	do i=3,imax-1
	factor=a(i)/b(i-1)
	b(i)=b(i)-factor*c(i-1)
	F(i)=F(i)-factor*F(i-1)
	end do
	
	!
    ! Backward substitution
	!
	
	dqdt3(imax-1)=F(imax-1)/b(imax-1)
	dqdt3(imax-2)=(F(imax-2)-c(imax-2)*dqdt3(imax-1))/b(imax-2)
    
	do i=imax-2,2,-1
  	dqdt3(i)=(F(i)-c(i)*dqdt3(i+1))/b(i)
	end do
	
    !
    ! boundary conditions
    !
	
       dzdt3(1)=4.0d0*dzdt3(2)-6.0d0*dzdt3(3)+4.0d0*dzdt3(4)-dzdt3(5)
       dzdt3(imax)=4.0d0*dzdt3(imax-1)-6.0d0*dzdt3(imax-2)+4.0d0*dzdt3(imax-3)-dzdt3(imax-4)
       
       dqdt3(1)=0.0d0
       dqdt3(0)=0.0d0
       dqdt3(-1)=0.0d0
       dqdt3(imax)=0.0d0
       dqdt3(imax+1)=0.0d0
       dqdt3(imax+2)=0.0d0
	   
	!
	! fourth step in RK4
	!

        do i=2,imax-1

        z1=zeta(i)+0.5d0*dt*dzdt3(i)
        ze=zeta(i+1)+0.5d0*dt*dzdt3(i+1)
        zw=zeta(i-1)+0.5d0*dt*dzdt3(i-1)

        h1=hs(i)+z1
        he=hs(i+1)+ze
        hw=hs(i-1)+zw
        
        q1=q(i)+0.5d0*dt*dqdt3(i)
        qe=q(i+1)+0.5d0*dt*dqdt3(i+1)
        qw=q(i-1)+0.5d0*dt*dqdt3(i-1)
        qee=q(i+2)+0.5d0*dt*dqdt3(i+2)
        qww=q(i-2)+0.5d0*dt*dqdt3(i-2)
		
	!
	! continuity equation
	!
	
        dzdt4(i)=-(qe-qw)/2.0d0/dx 
		
	!
	! momentum equation
	!
	
     u1=q1/h1     						
     taubx=cf*rho*u1*dabs(u1)
    
     part1=-g*h1*(ze-zw)/2.0d0/dx
     part2=-Ia*(qe*qe/he-qw*qw/hw)/2.0d0/dx
     t1=(qe-qw)/2.0d0/dx*(qe/he-2.0d0*q1/h1+qw/hw)/dx/dx
     t2=2.0d0*(qe/he-qw/hw)/2.0d0/dx*(qe-2.0d0*q1+qw)/dx/dx
     t3=q1/h1*(-qww+2.0d0*qw-2.0d0*qe+qee)/2.0d0/dx/dx/dx
     part3=Ip2*hs(i)*hs(i)/2.0d0*(t1+t2+t3)
     F(i)=part1+part2+part3-taubx/rho  
          
     a(i) = -(Ip10*hs(i)*hs(i))/2.0d0 +Ip11*hs(i)*z1/dx/dx
     b(i) = 1.0d0+2.0d0*(Ip10*hs(i)*hs(i))/2.0d0 +Ip11*hs(i)*z1/dx/dx
     c(i) = -(Ip10*hs(i)*hs(i))/2.0d0 +Ip11*hs(i)*z1/dx/dx
	end do
     
	do i=3,imax-1
	factor=a(i)/b(i-1)
	b(i)=b(i)-factor*c(i-1)
	F(i)=F(i)-factor*F(i-1)
	end do
	
	!
    ! Backward substitution
	!
	
	dqdt4(imax-1)=F(imax-1)/b(imax-1)
	dqdt4(imax-2)=(F(imax-2)-c(imax-2)*dqdt4(imax-1))/b(imax-2)
    
	do i=imax-2,2,-1
  	dqdt4(i)=(F(i)-c(i)*dqdt4(i+1))/b(i)
	end do
	
    !
    ! boundary conditions
    !
	
       dzdt4(1)=4.0d0*dzdt4(2)-6.0d0*dzdt4(3)+4.0d0*dzdt4(4)-dzdt4(5)
       dzdt4(imax)=4.0d0*dzdt4(imax-1)-6.0d0*dzdt4(imax-2)+4.0d0*dzdt4(imax-3)-dzdt4(imax-4)
       
       dqdt4(1)=0.0d0
       dqdt4(0)=0.0d0
       dqdt4(-1)=0.0d0
       dqdt4(imax)=0.0d0
       dqdt4(imax+1)=0.0d0
       dqdt4(imax+2)=0.0d0

    !
    ! RK4 algorithm
    !
	
       do i=1,imax
       znew(i)=zeta(i)+dt*(dzdt1(i)+2.0d0*dzdt2(i)+2.0d0*dzdt3(i)+dzdt4(i))/6.0d0      
       qnew(i)=q(i)+dt*(dqdt1(i)+2.0d0*dqdt2(i)+2.0d0*dqdt3(i)+dqdt4(i))/6.0d0 
       end do

	!
	! boundary conditions for dependent variables, zeta and q
	!
	
  		qnew(0)=0.0d0
        qnew(-1)=0.0d0

  		qnew(imax+1)=0.0d0
  		qnew(imax+2)=0.0d0



    end subroutine calc

!******************************************************************

    subroutine update (imax,q,zeta,h,hs,znew,qnew)
    implicit none
    double precision :: q(-1:imax+2),zeta(imax),h(imax),hs(imax)
	double precision :: znew(imax),qnew(-1:imax+2)
	integer i,imax

        do i=-1,imax+2
        q(i)=qnew(i)
        end do

  		do i=1,imax
    	zeta(i)=znew(i)
        h(i)=hs(i)+zeta(i)
  		end do

    end subroutine update

