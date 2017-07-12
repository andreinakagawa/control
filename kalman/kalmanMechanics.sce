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
tf = 100;
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
    ynoise = y + rand(1,'uniform');
    x0 = x;
    xk = [xk x];
    yk = [yk y];
    yn = [yn ynoise];
end
//------------------------------------------------------------------------------
//State estimation
//Kalman filter - using scilab function
q = 1;
r = [0.01 0; 0 0.01];
p0 = [100 0; 0 100];
x0 = [0;0];
x1 = x0;
p1 = p0;
for k=1:length(t)-1
    [x1(:,k+1),p1,x,p] = kalm(yn(:,k),x1(:,k),p1,discSys.A,discSys.B,discSys.C,q,r);
end
//------------------------------------------------------------------------------
//State estimation
//Kalman filter - using scilab algorithm
//------------------------------------------------------------------------------
figure();
plot(t,xk(1,:),'b');
plot(t,x1(1,:),'r');
title('Position');
figure();
plot(t,xk(2,:),'b');
plot(t,x1(2,:),'r');
title('Velocity');
