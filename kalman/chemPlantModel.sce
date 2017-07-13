//    FEDERAL UNIVERSITY OF UBERLANDIA
//    Biomedical Engineering Lab
//    Uberlandia, Brazil
//--------------------------------------------------    
//    Author: Andrei Nakagawa, MSc
//    contact: andrei.ufu@gmail.com
//--------------------------------------------------
//    Control II
//    Exercise I
//    Finding steady-state values
//--------------------------------------------------
//--------------------------------------------------
clear;
clc;
//--------------------------------------------------
//Model equations
function [f] = modelEq(x,Ca0,rho,Cp,dH,Qk,k0,Ear,Vr,u)
    Ca = x(1); //Concentration of A
    Cb = x(2); //Concentration of B
    T = x(3); //Reactor jacket temperature
    //Retrieving the values of K
    K = tempEq(T,k0,Ear);            
    f(1) = u(1) * (Ca0 - Ca) - K(1)*Ca - K(3)*Ca^2;
    f(2) = -Cb*u(1) + K(1)*Ca - K(2)*Cb;
    f(3) = -1/(rho*Cp)*(K(1)*Ca*dH(1) + K(2)*Cb*dH(2) + K(3)*Ca*Ca*dH(3)) + u(1)*(u(2)-T) + Qk/(rho*Cp*Vr);
endfunction
//--------------------------------------------------
//Temperature
function [k] = tempEq(T,k0,Ear)
    for i=1:3
        k(i) = k0(i) * exp(-Ear(i) / (T+273.15));
    end
endfunction
//--------------------------------------------------
function [f] = integModel(t,x,Ca0,rho,Cp,dH,Qk,k0,Ear,Vr,u)
    f = modelEq(x,Ca0,rho,Cp,dH,Qk,k0,Ear,Vr,u);
endfunction
//--------------------------------------------------
//Linear model
function [A,B]=linearModel(x,Ca0,rho,Cp,dH,Qk,k0,Ea_R,Vr,u)
    Ca=x(1); Cb=x(2); T=x(3);
    A(1,1)=-2*Ca*k0(3)*exp(-Ea_R(3)/(T + 273.15)) - k0(1)*exp(-Ea_R(1)/(T + 273.15)) - u(1);
    A(1,2)=0;
    A(1,3)=-Ca^2*Ea_R(3)*k0(3)*exp(-Ea_R(3)/(T + 273.15))/(T + 273.15)^2 - Ca*Ea_R(1)*k0(1)*exp(-Ea_R(1)/(T + 273.15))/(T + 273.15)^2;
    A(2,1)=k0(1)*exp(-Ea_R(1)/(T + 273.15));
    A(2,2)=-k0(2)*exp(-Ea_R(2)/(T + 273.15)) - u(1);
    A(2,3)=Ca*Ea_R(1)*k0(1)*exp(-Ea_R(1)/(T + 273.15))/(T + 273.15)^2 - Cb*Ea_R(2)*k0(2)*exp(-Ea_R(2)/(T + 273.15))/(T + 273.15)^2;
    A(3,1)=-(2*Ca*dH(3)*k0(3)*exp(-Ea_R(3)/(T + 273.15)) + dH(1)*k0(1)*exp(-Ea_R(1)/(T + 273.15)))/(Cp*rho);                                           
    A(3,2)=-dH(2)*k0(2)*exp(-Ea_R(2)/(T + 273.15))/(Cp*rho);
    A(3,3)= -u(1) - (Ca^2*Ea_R(3)*dH(3)*k0(3)*exp(-Ea_R(3)/(T + 273.15))/(T + 273.15)^2 + Ca*Ea_R(1)*dH(1)*k0(1)*exp(-Ea_R(1)/(T + 273.15))/(T + 273.15)^2 + Cb*Ea_R(2)*dH(2)*k0(2)*exp(-Ea_R(2)/(T + 273.15))/(T + 273.15)^2)/(Cp*rho);
    B=[-Ca + Ca0,0; -Cb,0; -T + u(2),u(1)];
endfunction
//--------------------------------------------------
//Parameters
rho = 0.9342;
Cp = 3.01;
Qk = -4496;
Vr = 0.01 * 10^3;
H1 = 4.2;
H2 = -11;
H3 = -41.85;
Ca0 = 5.1;
Ear = [9758.3,9758.3,8560];
k0 = [1.287e12,1.287e12,9.043e9];
dH = [H1;H2;H3];
//Input
u10 = 18.83;
u20 = 130;
u0 = [u10;u20];
//--------------------------------------------------
//Finding all the stationary points
u1 = 5:0.1:160;
x0 = [1;1;130];
x2 = [];
for i=1:length(u1)
    xss = fsolve(x0,list(modelEq,Ca0,rho,Cp,dH,Qk,k0,Ear,Vr,[u1(i);u20]));
    x2 = [x2 xss(2)];
end
//--------------------------------------------------
//Linear model
//Stationary point
uss = [18.83;130];
xss = fsolve(x0,list(modelEq,Ca0,rho,Cp,dH,Qk,k0,Ear,Vr,uss));
[Ac,Bc] = linearModel(xss,Ca0,rho,Cp,dH,Qk,k0,Ear,Vr,uss);
C = eye(3,3);
contSys = syslin('c',Ac,Bc,C);
//--------------------------------------------------
//Discrete model
dt = 0.01;
discSys = dscr(contSys,dt);
//--------------------------------------------------
//Simulation parameters
t0=0; tf=1;
t=t0:dt:tf;
//--------------------------------------------------
//--------------------------------------------------
//Adding noise
rand('seed',2);rand('normal');
//state uncertainty
w = rand(3,length(t))/1000;
for i=1:3
    w(i,:) = w(i,:) - mean(w(i,:));
end
//measurement noise
v = rand(3,length(t))/100;
for i=1:3
    v(i,:) = v(i,:) - mean(v(i,:));
end
//input noise covariance matrix
q = [mean(w(1,:)*w(1,:)') 0 0; 0 mean(w(2,:)*w(2,:)') 0; 0 0 mean(w(3,:)*w(3,:)')];
//output noise covariance matrix
r = [mean(v(1,:)*v(1,:)') 0 0; 0 mean(v(2,:)*v(2,:)') 0; 0 0 mean(v(3,:)*v(3,:)')];
//--------------------------------------------------
//Simulation
amp = 0.8;
freq=1;
usine = amp * sin(2*%pi*freq*t);
Xk = []; //without noise
Xkn = []; //with noise
uk = [1;0];
xk = [0;0;0];
xkn = [0;0;0];
Yk = []
//Bd = zeros(3,3);
//Bd(:,1:2) = discSys.B;
//discSys.B = Bd;
Xkalman = [];
q = [0 0; 0 0];
r = [100 0 0; 0 100 0; 0 0 100];
p0 = chol([100 0 0; 0 100 0; 0 0 100]);
p1 = p0;
x1 = x0;
for k=1:length(t)
    //model without noise
    xk = (discSys.A * xk) + (discSys.B * uk);
    Xk = [Xk xk];
    //model with noise
    xkn = (discSys.A * xkn) + (discSys.B * uk);
    Xkn = [Xkn xkn];
    //output with noise
    yk = discSys.C*xkn + v(:,k);
    Yk = [Yk yk];
    //kalman filter
    [x1,p1] = kalm(yk,x1,p1,discSys.A,discSys.B,discSys.C,q,r);
    Xkalman = [Xkalman x1];
end
//------------------------------------------------------------------------------
figure();
plot(t,Xk(1,:),'r');
plot(t,Xk(2,:),'g');
plot(t,Xk(3,:),'b');
//------------------------------------------------------------------------------
figure();
plot(t,Yk(1,:),'r');
plot(t,Yk(2,:),'g');
plot(t,Yk(3,:),'b');
//------------------------------------------------------------------------------
figure();
plot(t,Xkalman(1,:),'r');
plot(t,Xkalman(2,:),'g');
plot(t,Xkalman(3,:),'b');
//------------------------------------------------------------------------------
