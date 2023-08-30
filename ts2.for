       dimension x(25000000), y(25000000), t(25000000), l(25000000)
       dimension x1(25000000), y1(25000000), xd1(25000000),
     * yd1(25000000)
       dimension x2(25000000), x3(25000000), ld(25000000), y3(25000000)
       dimension z1(25000000), yd2(25000000), yd3(25000000),
     * xd2(25000000),xd3(25000000),y2(25000000)
       dimension zlp(25000000)
       real t,x1,l,y1,xd1,yd1,ld,x1delay,tau1,beta,delta,xi,h,x2,y3,
     * x,y,x3,z1,yd2,yd3,zlp
       real sigma1,sigma2,sigma3,g,gamma1,alpha,y2,rho
       integer n,kn,i,j,kdnp,m,kn1,ik,index,cont

       open(51,file='delayseriesme1.dat')

       x1(1)=0.3
       x2(1)=0.3
       x3(1)=0.3

       x1delay=0.25
       x2delay=0.25
       x3delay=0.25

       sigma1 = 0.52
       sigma2 = 0.72
       sigma3 = 0.41

       xi= 0.50
       beta = 2.50
       delta = 1.39
       tau1 = 0.019
       
       h=0.01

       do 100 n=1,10000000
       t(n)=n*h
       l(n)=t(n)-tau1
       kn=aint(l(n)/h,kind(l(n)/h))

       if(l(n) .le. 0.0)then
       y1(n)=x1delay
       else
       y1(n)=((x1(kn+1)-x1(kn))*(l(n)-kn*h)/h)+x1(kn)
       endif

       if(l(n) .le. 0.0)then
       y2(n)=x2delay
       else
       y2(n)=((x2(kn+1)-x2(kn))*(l(n)-kn*h)/h)+x2(kn)
       endif

       if(l(n) .le. 0.0)then
       y3(n)=x3delay
       else
       y3(n)=((x3(kn+1)-x3(kn))*(l(n)-kn*h)/h)+x3(kn)
       endif

       xd1(n+1)=x1(n)+(x1(n)*((1-sigma1)*y1(n)+(1-sigma1)*y2(n)
     *-sigma1*y3(n)+(sigma1-xi)))*h
       xd2(n+1)=x2(n)+(x2(n)*((1-sigma2)*y1(n)+(1-sigma2)*y2(n)
     *-(delta+sigma2)*y3(n)+(sigma2-xi)))*h
       xd3(n+1)=x3(n)+(x3(n)*((beta-sigma3)*y1(n)+(beta-delta-sigma3)
     **y2(n)-sigma3*y3(n)+(sigma3-xi)))*h

       t(n+1)=(n+1)*h
       ld(n+1)=t(n+1)-tau1
       kdnp=aint(ld(n+1)/h,kind(ld(n+1)/h))

       if(ld(n+1) .le. 0.0) then
       yd1(n+1)=x1delay
       endif
       if(ld(n+1) .gt. 0.0 .and. ld(n+1) .le. t(n)) then
       yd1(n+1)=((x1(kdnp+1)-x1(kdnp))*(ld(n+1)-kdnp*h)/h)+x1(kdnp)
       endif
       if(ld(n+1) .gt. t(n) .and. ld(n+1) .lt. t(n+1)) then
       yd1(n+1)=((xd1(n+1)-x1(n))*(ld(n+1)-t(n)*h)/h)+x1(n)
       endif
       if(ld(n+1) .eq. t(n+1)) then
       yd1(n+1)=xd1(n+1)
       endif

       if(ld(n+1) .le. 0.0) then
       yd2(n+1)=x2delay
       endif
       if(ld(n+1) .gt. 0.0 .and. ld(n+1) .le. t(n)) then
       yd2(n+1)=((x2(kdnp+1)-x2(kdnp))*(ld(n+1)-kdnp*h)/h)+x2(kdnp)
       endif
       if(ld(n+1) .gt. t(n) .and. ld(n+1) .lt. t(n+1)) then
       yd2(n+1)= ((xd2(n+1)-x2(n))*(ld(n+1)-t(n)*h)/h)+x2(n)
       endif
       if(ld(n+1) .eq. t(n+1)) then
       yd2(n+1)=xd2(n+1)
       endif

       if(ld(n+1) .le. 0.0) then
       yd3(n+1)=x3delay
       endif
       if(ld(n+1) .gt. 0.0 .and. ld(n+1) .le. t(n)) then
       yd3(n+1)=((x3(kdnp+1)-x3(kdnp))*(ld(n+1)-kdnp*h)/h)+x3(kdnp)
       endif
       if(ld(n+1) .gt. t(n) .and. ld(n+1) .lt. t(n+1)) then
       yd3(n+1)= ((xd3(n+1)-x3(n))*(ld(n+1)-t(n)*h)/h)+x3(n)
       endif
       if(ld(n+1) .eq. t(n+1)) then
       yd3(n+1)=xd3(n+1)
       endif


       x1(n+1)=x1(n)+(x1(n)*((1-sigma1)*y1(n)+(1-sigma1)*y2(n)-sigma1
     **y3(n)+(sigma1-xi)))*h/2.0+(xd1(n+1)*((1-sigma1)*yd1(n+1)
     *+(1-sigma1)*yd2(n+1)-sigma1*yd3(n+1)+(sigma1-xi)))*h/2.0
       x2(n+1)=x2(n)+(x2(n)*((1-sigma2)*y1(n)+(1-sigma2)*y2(n)
     *-(delta+sigma2)*y3(n)+(sigma2-xi)))*h/2.0+(xd2(n+1)*((1-sigma2)
     **yd1(n+1)+(1-sigma2)*yd2(n+1)-(sigma2+delta)*yd3(n+1)
     *+(sigma2-xi)))*h/2.0
       x3(n+1)=x3(n)+(x3(n)*((beta-sigma3)*y1(n)+(beta-delta-sigma3)
     **y2(n)-sigma3*y3(n)+(sigma3-xi)))*h/2.0+(xd3(n+1)*((beta-sigma3)
     **yd1(n+1)+(beta-delta-sigma3)*yd2(n+1)-sigma3*yd3(n+1)
     *+(sigma3-xi)))*h/2.0


       if(n .gt. 9800000)then
       write(*,*)t(n+1),x1(n+1),x2(n+1),x3(n+1), x1(n+1)+x2(n+1)+x3(n+1)
       write(51,*)t(n+1),x1(n+1),x2(n+1),x3(n+1),x1(n+1)+x2(n+1)+x3(n+1)
     * ,yd1(n+1),yd2(n+1),yd3(n+1),yd1(n+1)+yd2(n+1)+yd3(n+1),
     *x1(n+1)/(x1(n+1)+x2(n+1)+x3(n+1)),x2(n+1)/(x1(n+1)+x2(n+1)
     * +x3(n+1)), x3(n+1)/(x1(n+1)+x2(n+1)+x3(n+1))
       endif

100    continue
       stop
       end
