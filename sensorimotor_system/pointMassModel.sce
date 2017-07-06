//------------------------------------------------------------------------------
// FEDERAL UNIVERSITY OF UBERLANDIA
// Faculty of Electrical Engineering
// Biomedical Engineering Lab
// Uberlandia, Brazil
//------------------------------------------------------------------------------
// Author: Andrei Nakagawa, MSc
// Contact: andrei.ufu@gmail.com
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//States
//Position and velocity in X
//Position and velocity in Y
//Inputs
//Force in X and Y
function [Ac,Bc,Cc] = pointMassModel(m)
    Ac = [0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0];
    Bc = [0 0; 1/m 0; 0 0; 0 1/m];
    Cc = eye(size(Ac,1),size(Ac,2))
endfunction
//------------------------------------------------------------------------------
[A,B,C] = pointMassModel(1);
//------------------------------------------------------------------------------
t0=0;
tf=5;
dt = 0.01;
t = t0:dt:tf;
contSys = syslin('c',A,B,C);
discSys = dscr(contSys,dt);
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//Weight matrices
Qd=diag([10,0,10,0]);
Rd=diag([0.2,0.2]);
//Discrete riccati
Ad = discSys(2); //A
Bd = discSys(3); //B
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//Differential Riccati Equation - Discrete-time
//Calculating the solution to riccati for each instant in time
//and then finding the time-varying gain for each time step
//------------------------------------------------------------------------------
Sdisc = [];
Kdisc = [];
S0 = diag([120,0,120,0]); //Estimate for the Riccati matrix 
for k=1:length(t)-1
    //Calculating the time-varying gain
    K = inv(Bd'*S0*Bd + Rd)*(Bd'*S0*Ad);
    //New riccati solution
    S0 = Ad'*S0*Ad - Ad'*S0*Bd*((Rd + Bd'*S0*Bd)^-1)*Bd'*S0*Ad + Qd;
    //Stores the riccati solution
    Sdisc = [Sdisc S0];    
    //Stores the gain
    Kdisc = [Kdisc K];
end
//------------------------------------------------------------------------------
cont = 1;
//Desired setpoints or reference trajectory
xd = [0;0;6;0];
xint = []; //stores all the states during integration
uint = []; //stores all the inputs during integration
costQ = [];
costR = [];
x = [0;0;0;0]; //temporary variable for storing states
for k=1:length(t)-1
    //Calculating the input
    u = -Kdisc(:,cont:cont+3) * (x-xd);
    //Calculating the new states
    x = Ad*x + Bd*u;
    //Storing the new states
    xint = [xint x];
    //Storing the new inputs
    uint = [uint u];
    //Stores the cost in this step
    costQ = [costQ (x-xd)'*Qd*(x-xd)];
    //Stores the cost in this step
    costR = [costR u'*Rd*u];
    //Increments the counter to loop through the gain matrix
    cont = cont + 4;
end
//------------------------------------------------------------------------------
t = t(1:$-1);
//------------------------------------------------------------------------------
figure();
plot(xint(1,:),xint(3,:));
ax=gca();
ax.data_bounds=[-6 -6; 6 6];
figure();
plot(t,xint(1,:),'r');
plot(t,xint(2,:),'g');
plot(t,uint(1,:),'k');
figure();
plot(t,xint(3,:),'r');
plot(t,xint(4,:),'g');
plot(t,uint(2,:),'k');
figure();
plot(t,costQ,'r');
plot(t,costR,'b');
//------------------------------------------------------------------------------
