//------------------------------------------------------------------------------
//FEDERAL UNIVERSITY OF UBERLANDIA
//Discipline: Control II
//Project 1:
//Botan et al.: "Discrete Time Linear Quadratic Optimal Control for an Electrical Servo Drive System"
//Authors: Andrei Nakagawa and Henrique Oyama
//Last update: May, 2017
//Uberlandia, Brazil
//------------------------------------------------------------------------------
//Description: Time-varying LQR
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
clear;
clc;
//------------------------------------------------------------------------------
//Solution to the riccati differential equation
function dxdt = mricc(t,x,A,B,Q,R)
    dxdt = A'*x + x*A - x*B*inv(R)*B'*x + Q;
endfunction
//------------------------------------------------------------------------------
//Function that integrates the time-varying optimal control
//with the model
function dxdt = optControlTV(t,x,xd,A,B,Q,R)
    //Retrieves the current state
    xs = x(1,:)'; 
    //Retrieves the current Riccati solution
    pr = x(2:4,:); 
    //Time-varying gain
    K = inv(R)*B'*P; 
    //Control    
    u = -K*(xs - xd); 
    //New states
    xs = A*xs + B*u; 
    //New Riccati solution
    P = A'*pr + pr*A  - pr*B*inv(R)*B'*pr + Q;
    //Output of the function
    dxdt = zeros(5,3);
    dxdt(1,:) = xs'; //states
    dxdt(2:4,:) = P; //riccati
    dxdt(5,1) = u; //input
endfunction
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
dt = 0.01; //step
t0 = 0; //initial time
tf = 50; //final time
t = t0:dt:tf; //time vector
//------------------------------------------------------------------------------
//Model matrices
[A,B] = dcMotorModel(r,pho,J,Cm,Ce,L,Re);
//------------------------------------------------------------------------------
//Desired setpoints or reference trajectory
xd = [2;0;0];
//Weight matrices
Q=diag([10,0,3.1]);
R=0.3;
//Solving P -- Riccati equation
P=riccati(A,B*inv(R)*B',Q,'c','eigen');
//Kalman gain
Kcont=inv(R)*B'*P;
//------------------------------------------------------------------------------
//Integrating the model
//x will have the states, the riccati solution (matrix) and u
x0 = zeros(5,3);
x0(1,:) = [0,0,0]; //initial states
x0(2:4,:) = diag([60,0,2]); //riccati
x0(5,:) = 0; //input
xint = ode(x0,t0,t,list(optControlTV,xd,A,B,Q,R));
//------------------------------------------------------------------------------
//Retrieves only the states from the output of the integration
xs = xint(1,:);
ids = 1:3:length(xs);
//Retrieves only the inputs from the output of the integration
uo = xint(5,:);
uid = 1:3:length(uo);
//------------------------------------------------------------------------------
//Plots
figure();
plot(t,xs(ids),'r'); //position
plot(t,xs(ids+1),'b'); //angular velocity
plot(t,xs(ids+2),'g'); //current
plot(t,uo(ids),'k'); //voltage
legend({'Angular position (x1)', 'Angular velocity (x2)', 'Current (x3)', 'Voltage (u)'},-1);
xlabel('Time (s)');
title('Simulation results');
xs2jpg(gcf(), 'simulationResultsTV.jpg'); // Export to a JPG file
//------------------------------------------------------------------------------
