//------------------------------------------------------------------------------
// FEDERAL UNIVERSITY OF UBERLANDIA
// Faculty of Electrical Engineering
// Biomedical Engineering Lab
// Uberlandia, Brazil
//------------------------------------------------------------------------------
// Author: Andrei Nakagawa, MSc
// Contact: andrei.ufu@gmail.com
//------------------------------------------------------------------------------
//Simulation parameters
t0 = 0;
tf = 50;
dt = 0.01;
t = t0:dt:tf;
//------------------------------------------------------------------------------
//Continuous-time model
f = [0 1; 0 0];
g = [0;1];
h = [1 0; 0 1];
contSys = syslin('c',f,g,h);
//------------------------------------------------------------------------------
//Discrete-time model
discSys = dscr(contSys,dt);
//------------------------------------------------------------------------------
//Running a simulation
u=1;
xk = []
yk = []
yn = []
x0 = [0;0];
for k=1:length(t)
    x = discSys.A * x0 + discSys.B  * u;
    y = discSys.C * x;
    ypn = y(1) + rand(1,'uniform');
    yvn = y(2) + rand(1,'uniform');
    ynoise = [ypn;yvn];
    x0 = x;
    xk = [xk x];
    yk = [yk y];
    yn = [yn ynoise];
end
//------------------------------------------------------------------------------
//State estimation
//Kalman filter - using scilab function
q = 1;
r = [0.001 0; 0 0.01];
p0 = [1 0; 0 1];
x0 = [0;0];
x1 = x0;
p1 = p0;
for k=1:length(t)-1
    [x1(:,k+1),p1,x,p] = kalm(yn(:,k),x1(:,k),p1,discSys.A,discSys.B,discSys.C,q,r);
end
//------------------------------------------------------------------------------
//State estimation
//Kalman filter - using scilab algorithm
q = 1;
r = [0.001 0; 0 0.01];
p0 = [1 0; 0 1];
x0 = [0;0];
xk1 = [x0];
//kalman
for k=1:length(t)-1
    //measurement update (correction)
    E = yn(:,k) - discSys.C * x0;
    S = discSys.C * p0 * discSys.C' + r;
    K = p0*discSys.C' * inv(S);
    x = x0 + K*E;
    p = p0 - K*discSys.C*p0;
    //time update (prediction)
    xp = discSys.A * x;
    pp = discSys.A * p * discSys.A' + discSys.B * q * discSys.B';
    //stores the state estimate
    xk1 = [xk1 xp];
    x0 = xp;
    p0 = pp;
end
//------------------------------------------------------------------------------
figure();
plot(t,xk(1,:),'b');
plot(t,x1(1,:),'r');
plot(t,xk1(1,:),'g');
title('Position');
figure();
plot(t,xk(2,:),'b');
plot(t,x1(2,:),'r');
plot(t,xk1(2,:),'g');
title('Velocity');
