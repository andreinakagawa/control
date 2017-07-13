//------------------------------------------------------------------------------
// FEDERAL UNIVERSITY OF UBERLANDIA
// Faculty of Electrical Engineering
// Biomedical Engineering Lab
// Uberlandia, Brazil
//------------------------------------------------------------------------------
// Author: Andrei Nakagawa, MSc
// Contact: andrei.ufu@gmail.com
//------------------------------------------------------------------------------
clc
clear
mode(-1);
//------------------------------------------------------------------------------
//  Definição  de  funções
function  [f]=fun(x,Ca0,T0,rho,Cp,dH,Qk,k0,Ea_R,Vr,u)
 // variáveis x(1)=Ca, x(2)=Cb, x(3)=T; //x(4)=Tk 
 Ca=x(1); Cb=x(2); T=x(3);
 K=K_fun(T,k0,Ea_R);
 f(1)=u(1)*(Ca0-Ca)-K(1)*Ca-K(3)*Ca^2;
 f(2)=-u(1)*Cb+K(1)*Ca-K(2)*Cb;
 f(3)= -1/(rho*Cp)*(K(1)*Ca*dH(1)+K(2)*Cb*dH(2)+K(3)*Ca^2*dH(3))+u(1)*(u(2)-T)+Qk/rho/Cp/Vr;
endfunction
//------------------------------------------------------------------------------
function K=K_fun(T,k0,Ea_R)
    for i=1:3
       K(i)=k0(i)*exp(-Ea_R(i)/(T+273.15)); 
    end    
endfunction  
//------------------------------------------------------------------------------
function [A,B]=linear(x,Ca0,T0,rho,Cp,dH,Qk,k0,Ea_R,Vr,u)
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
//----------------------------------
function  [dydt]=fun_din(t,x,Ca0,T0,rho,Cp,dH,Qk,k0,Ea_R,Vr,u)
   dydt=fun(x,Ca0,T0,rho,Cp,dH,Qk,k0,Ea_R,Vr,u)
endfunction
//------------------------------------------------------------------------------
function [xplus,pplus] = kalman(F,G,H,Q,R,xk,pk,yk,uk)
    //time update (prediction)
    xminus = F*xk + G*uk;
    pminus = F*pk*F' + G*Q*G';    
    //measurement update (correction)    
    measureError = yk - (H*xminus);
    kalmanGain = pminus*H' * inv(H*pminus*H' + R);
    //state estimate
    xplus = xminus + kalmanGain * measureError;
    //covariance
    pplus = pminus - kalmanGain*H*pminus;
endfunction
//------------------------------------------------------------------------------
//Model parameters
 Ca0=5.1;// mol/L
 Qk=-4496;// kJ/h
 T0=130;// oC
 rho=0.9342;
 Cp=3.01;
 k0=[1.287e12;1.287e12;9.043e9];
 Ea_R=[9758.3;9758.3;8560];// K
 dH=[4.2;-11.0;-41.85];
 Vr=10;//L
 F_Vr=18.83; //  h^(-1)
 u=[F_Vr;T0];
 //-----------------------------------------------------------------------------
 //Time parameters
 t0=0;
 tf=1;
 dt=0.001;
 t = t0:dt:tf;
 //-----------------------------------------------------------------------------
//  Estimativa  para  a  solução
x0=[1;1;130];
// Processamento
[x_ss,valorf,iflag]=fsolve(x0,list(fun,Ca0,T0,rho,Cp,dH,Qk,k0,Ea_R,Vr,u));
//------------------------------------------------------------------------------
//Continuous-time system
[Ac,Bc]=linear(x_ss,Ca0,T0,rho,Cp,dH,Qk,k0,Ea_R,Vr,u);
Cc = [1 0 0; 0 1 0; 0 0 1];
contSys = syslin('c',Ac,Bc,Cc);
//------------------------------------------------------------------------------
//Discrete-time system
discSys = dscr(contSys,dt);
//------------------------------------------------------------------------------
//Simulating the plant dynamics
ynoise = [];
x0 = [0;0;0];
xk = [];
u1 = 1*sin(2*%pi*5*t);
u = [1;0];
for k=1:length(t)
    x = discSys.A*x0 + discSys.B*u;
    xk = [xk x];
    yn1 = x(1) + (rand(1,'normal')/500);
    yn2 = x(2) + (rand(1,'normal')/500);
    yn3 = x(3) + (rand(1,'normal')/50);
    yn = [yn1;yn2;yn3];
    ynoise = [ynoise yn];
    x0 = x;
end
//------------------------------------------------------------------------------
//State estimation
//Kalman filter - using custom-made algorithm
q = [10 0; 0 10];
r = [1 0 0; 0 1 0; 0 0 10];
p0 = [1 0 0; 0 1 0; 0 0 1];
x0 = [0;0;0];
x1=x0;
p1=p0;
xk1 = [x0];
//kalman
for k=1:length(t)
    [x1,p1] = kalman(discSys.A,discSys.B,discSys.C,q,r,x0,p0,ynoise(:,k),u);
    xk1 = [xk1 x1];
    x0 = x1;
    p0 = p1;
end
xk1 = xk1(:,1:$-1);
pk1 = p0;
//------------------------------------------------------------------------------
figure();
plot(t,ynoise(1,:),'g');
plot(t,xk(1,:),'b');
plot(t,xk1(1,:),'k');
title('x1');
legend({'y', 'x', 'xest'},-1);
xs2jpg(gcf(), 'simulation_kalman_x1.jpg'); // Export to a JPG file
figure();
plot(t,ynoise(2,:),'g');
plot(t,xk(2,:),'b');
plot(t,xk1(2,:),'k');
legend({'y', 'x', 'xest'},-1);
title('x2');
xs2jpg(gcf(), 'simulation_kalman_x2.jpg'); // Export to a JPG file
figure();
plot(t,ynoise(3,:),'g');
plot(t,xk(3,:),'b');
plot(t,xk1(3,:),'k');
legend({'y', 'x', 'xest'},-1);
title('x3');
xs2jpg(gcf(), 'simulation_kalman_x3.jpg'); // Export to a JPG file
//------------------------------------------------------------------------------
