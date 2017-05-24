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
function dxdt=integModel(t,x,Ac,Bc)
  //u=-Kcont*(x-xd);
  if(t >= 1 & t < 2)
      u = 1;
  else
      u = 0;
   end
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
u = 2; //V
x0 = [0;0;0]; //initial state
t0 = 0;
tf = 5;
dt = 0.01;
t = t0:dt:tf;
[A,B] = dcMotorModel(r,pho,J,Cm,Ce,L,Re);
xint = ode(x0,t0,t,list(integModel,A,B));
figure();
plot(t,xint(1,:),'r');
plot(t,xint(2,:),'b');
plot(t,xint(3,:),'g');
xlabel('Time (s)');
legend({'Angular position', 'Angular velocity', 'Current'},-1);
