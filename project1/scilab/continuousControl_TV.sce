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
//and the time-varying gain
function dsk = tvGain(t,x,A,B,Q,R)
    S = x(1:3,:); //Retrieves the current riccati solution    
    K = inv(R)*B'*S; //new gain
    S = A'*S + S*A - S*B*inv(R)*B'*S + Q; //new riccati solution
    dsk = zeros(4,3); //solution of the function
    dsk(1:3,:) = S; //riccati
    dsk(4,:) = K; //gain
endfunction
//------------------------------------------------------------------------------
//Solution to the riccati differential equation
function dxdt = mricc(t,x,A,B,Q,R)
    dxdt = A'*x + x*A - x*B*inv(R)*B'*x + Q;
endfunction
//------------------------------------------------------------------------------
//Dynamic model for integration
function dxdt=integModel(t,x,xd,Ac,Bc,K)
  u=-K*(x-xd);
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
tf = 20; //final time
t = t0:dt:tf; //time vector
//------------------------------------------------------------------------------
//Model matrices
[A,B] = dcMotorModel(r,pho,J,Cm,Ce,L,Re);
//------------------------------------------------------------------------------
//Control parameters
Q = diag([10,0,3.1]); //weights the states
R = 0.6; //weights the controls
ST = diag([600,150,0]); //desired matrix at time T
//Calculating the time-varying gain
//initial conditions
sk0 = zeros(4,3);
sk0(1:3,:) = ST; //riccati
//Integrates first to find K
skint = ode(sk0,t0,t,list(tvGain,A,B,Q,R));
//------------------------------------------------------------------------------
//Second integration -- applying optimal control with time-varying gain
//Initial states
x0 = [0;0;0];
//Desired setpoints or reference trajectory
xd = [2;0;0];
//counter for looping through the gains
is=1;
//stores all the states
X = [];
//stores all the inputs
U = [];
//Integrating the dc motor model step-by-step
//in each step, a new gain will be used to update the control law
for i=1:length(t)-1
    //Retrieves the current gain
    K = skint(4,is:is+2);
    //Integrates the model from the current time until the next timestep
    //the function will apply optimal control given the time-varying gain
    xint = ode(x0,t(i),t(i+1),list(integModel,xd,A,B,K));
    //The initial state of the next integration will be the current
    //states
    x0 = xint;
    //Updates the counter for looping through the gains
    is = is+3;
    //Stores the input
    U = [U -K*(x0-xd)];
    //Stores the states
    X = [X xint];
end
//------------------------------------------------------------------------------
////Plots
//corrects the time vector since the last time step was not calculated
t = t(1:$-1); 
figure();
plot(t,X(1,:),'r'); //position
plot(t,X(2,:),'b'); //angular velocity
plot(t,X(3,:),'g'); //current
plot(t,U,'k'); //voltage
legend({'Angular position (x1)', 'Angular velocity (x2)', 'Current (x3)', 'Voltage (u)'},-1);
xlabel('Time (s)');
title('Simulation results');
xs2jpg(gcf(), 'simulationResultsTV.jpg'); // Export to a JPG file
////------------------------------------------------------------------------------
