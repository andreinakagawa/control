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
  u=-K*(x-xd);
  dxdt=Ac*x+Bc*u;
endfunction
//------------------------------------------------------------------------------
//Continuous-time model of the DC motor
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
tf = 5; //final time
dt = 0.002; //step
t = t0:dt:tf; //time vector
//------------------------------------------------------------------------------
//Model matrices
[A,B] = dcMotorModel(r,pho,J,Cm,Ce,L,Re);
//------------------------------------------------------------------------------
//Creating the discrete model
C = eye(3,3); //output matrix
//Continuous linear system
contSys = syslin('c',A,B,C);
//Discrete linear system
discSys = dscr(contSys,dt);
//------------------------------------------------------------------------------
//Desired setpoints or reference trajectory
xd = [2;0;0];
//Weight matrices
Qd=diag([10,0,3.1]);
Rd=0.2;
//Discrete riccati
F = discSys(2); //A
Bd = discSys(3); //B
G1=Bd;
G2=Rd;
G = G1/G2 * G1';
H = Qd;
Pdisc = ricc(F,G,H,'disc'); //solution to riccati
Kdisc = inv(Bd'*Pdisc*Bd + Rd)*(Bd'*Pdisc*F); //Kalman gain
//------------------------------------------------------------------------------
//Integrating the discrete model
//------------------------------------------------------------------------------
xint = []; //vector of states
uint = []; //vector of input (control)
x=[0;0;0]; //initial states
//For each time step, integrate the discrete-time model
for k=1:length(t)
    //Estimating the input given the Kalman gain 
    //calculated previously
   ud = -Kdisc*(x-xd);
   //Find the new states
   //Ax + Bu with A and B coming from the discrete model
   x = discSys(2)*x + discSys(3)*ud;
   //Saves the new states
   xint = [xint x];
   //Saves the new control
   uint = [uint ud];
end
//------------------------------------------------------------------------------
//Plots
figure();
plot(t,xint(1,:),'r');
plot(t,xint(2,:),'b');
plot(t,xint(3,:),'g');
plot(t,uint,'k');
title('Discrete-time Optimal Control');
xlabel('Time (s)');
legend({'Angular position (x1)', 'Angular velocity (x2)', 'Current (x3)', 'Voltage (u)'},-1);
xs2jpg(gcf(), 'simulation_discrete_LTI.jpg'); // Export to a JPG file
//------------------------------------------------------------------------------
