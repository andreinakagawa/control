//------------------------------------------------------------------------------
//FEDERAL UNIVERSITY OF UBERLANDIA
//Discipline: Control II
//Project 1:
//Botan et al.: "Discrete Time Linear Quadratic Optimal Control for an Electrical Servo Drive System"
//Authors: Andrei Nakagawa e Henrique Oyama
//Last update: May, 2017
//Uberlandia, Brazil
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
clear;
clc;
//------------------------------------------------------------------------------
//Dynamic model for integration
function dxdt=integModel(t,x,xd,Ac,Bc,K)
  u=-Kcont*(x-xd);
  dxdt=Ac*x+Bc*u;
endfunction
//------------------------------------------------------------------------------
//Model of the DC motor
function [Ac,Bc]=dcMotorModel(r,pho,J,Cm,Ce,L,Re)
  Ac(1,1)=0;
  Ac(1,2)=r;
  Ac(1,3)=0;
  Ac(2,1)=0;
  Ac(2,2)=-pho/J;
  Ac(2,3)=Cm/J;
  Ac(3,1)=0;
  Ac(3,2)=-Ce/L;
  Ac(3,3)=-Re/L;
  Bc=[0; 0; 1/L];
endfunction
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//Parameters
Re=3.1// Ohm
L=0.16//H
Ce=0.58//Vs/rad
Cm=0.58//Nm/A
J=0.028//Nms^2/rad
r=1/15
tau=0.003//s
tf=0.9//s
pho=0.001008//
m=0.78//Nm
//------------------------------------------------------------------------------
//Simulation parameters
x0 = [0;0;0]; //initial state
t0 = 0; //initial time
tf = 10; //final time
dt = 0.002; //step
t = t0:dt:tf; //time vector
//------------------------------------------------------------------------------
//Model matrices
[A,B] = dcMotorModel(r,pho,J,Cm,Ce,L,Re);
//------------------------------------------------------------------------------
//Desired setpoints or reference trajectory
xd = [2;0;0];
//Weight matrices
Q=diag([10,0,3.1]);
R=0.2;
//Solving P -- Riccati equation
P=riccati(A,B*inv(R)*B',Q,'c','eigen');
//Kalman gain
Kcont=inv(R)*B'*P;
//------------------------------------------------------------------------------
//Integrating the model
xint = ode(x0,t0,t,list(integModel,xd,A,B,Kcont));
//------------------------------------------------------------------------------
//Recovering the control policy over time
uc = []
for i=1:length(t)
    uc = [uc -Kcont*(xint(:,i)-xd)];
end
//------------------------------------------------------------------------------
//Plots
figure();
plot(t,xint(1,:),'r');
plot(t,xint(2,:),'b');
plot(t,xint(3,:),'g');
plot(t,uc,'k');
xlabel('Time (s)');
legend({'Angular position (x1)', 'Angular velocity (x2)', 'Current (x3)', 'Voltage (u)'},-1);
//------------------------------------------------------------------------------
