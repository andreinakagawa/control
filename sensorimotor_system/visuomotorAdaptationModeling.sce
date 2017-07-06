//------------------------------------------------------------------------------
// FEDERAL UNIVERSITY OF UBERLANDIA
// Faculty of Electrical Engineering
// Biomedical Engineering Lab
// Uberlandia, Brazil
//------------------------------------------------------------------------------
// Author: Andrei Nakagawa, MSc
// Contact: andrei.ufu@gmail.com
//------------------------------------------------------------------------------
//To do: Study iterative learning control and try to model
//adaptation
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
S0 = diag([500,0,120,0]); //Estimate for the Riccati matrix 
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
xd = [2;0;5;0];
xint = []; //stores all the states during integration
uint = []; //stores all the inputs during integration
costQ = []; //cost of states
costR = []; //cost of control
x = [2;0;0;0]; //temporary variable for storing states
//Perturbation
//Rotation matrix
Ck = [cos((30*%pi)/180) 0 -sin((30*%pi)/180) 0; sin((30*%pi)/180) 0 cos((30*%pi)/180) 0];
yk = [];
for k=1:length(t)-1
    //perturbation
    r = x - [2;0;0;0];
    r = Ck*r;
    r = r + [2;0];
    r = [r(1);r(2)];
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
    yk  = [yk r];
end
//------------------------------------------------------------------------------
t = t(1:$-1);
//------------------------------------------------------------------------------
figure();
plot(xint(1,:),xint(3,:),'r');
plot(yk(1,:),yk(2,:),'b');
ax=gca();
ax.data_bounds=[-15 -15; 15 15];
figure();
plot(t,xint(2,:),'r');
plot(t,xint(4,:),'g');
plot(t,uint(1,:),'k');
plot(t,uint(2,:),'k');
figure();
plot(t,costQ,'r');
plot(t,costR,'b');
//------------------------------------------------------------------------------
