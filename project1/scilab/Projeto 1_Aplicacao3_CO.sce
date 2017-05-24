clc
clear
mode(-1)
//Projeto 1: Discrete Time Linear Quadratic Optimal Control for an Electrical Servo Drive System
//Grupo: Andrei 
//       Henrique Carlos Oyama
// Maio 2017
// Definição  de  funções

//----------------------------------
//function  [dydt]=fun_din(t,x,Ca0,T0,rho,Cp,dH,Qk,k0,Ea_R,Vr,Kcont,uss,xss)
   //dydt=fun(x,Ca0,T0,rho,Cp,dH,Qk,k0,Ea_R,Vr,Kcont,uss,xss)
//endfunction   
//  ---------------------------------------------------------------
    function dxdt=planta(t,x,xd,Ac,Bc,Kcont)
        u=-Kcont*(x-xd);
        dxdt=Ac*x+Bc*u;
    endfunction    
//----------------------------------------------------------------------------------------
function [Ac,Bc]=linear(r,pho,J,Cm,Ce,L,Re)
    Ac(1,1)=0;
    Ac(1,2)=r;
    Ac(1,3)=0;
   // Ac(1,4)=0
    Ac(2,1)=0;
    Ac(2,2)=-pho/J;
    Ac(2,3)=Cm/J;
    //Ac(2,4)=-1/J;
    Ac(3,1)=0;                                           
    Ac(3,2)=-Ce/L;
    Ac(3,3)=-Re/L;
    //Ac(3,4)=0;
   // Ac(4,1)=0
   // Ac(4,2)=0
   // Ac(4,3)=0
    //Ac(4,4)=0
    Bc=[0; 0; 1/L]//; 0];
endfunction

//  Entrada de Dados
Re=3.1// Ohm
L=0.16//H
Ce=0.58//Vs/rad
Cm=0.58//Nm/A
J=0.028//Nms^2/rad
r=1/15
tau=0.003//s
tf=0.9//s
thetad=2//rad
pho=0.001008// VER!
m=0.78//Nm  VER!
xd=[4; 0; 0];
 //xss=[1.2348546; 0.9000188;134.13978];
 //uss=[18.83; 130];
[Ac,Bc]=linear(r,pho,J,Cm,Ce,L,Re)
Cc=[1 0 0; 0 1 0; 0 0 1];

 S=diag([600,150,0]);
 Q=diag([10,0,3.1]);
 R=0.6;

 //P1=ricc(Ac,Bc*inv(R)*Bc',Q,'cont');
 P=riccati(Ac,Bc*inv(R)*Bc',Q,'c','eigen');
 //disp(norm(P1-P2,1));
 Kcont=inv(R)*Bc'*P;
 Ac_cl=Ac-Bc*Kcont;
 Lamb_cl=spec(Ac_cl)
 x0=[0 ;0; 0];
 t0=0;
 tfim=50;
 delta_t=0.01;
 t=t0:delta_t:tfim;
 lista=list(planta,xd,Ac,Bc,Kcont);
 x=ode(x0,t0,t,lista);
 U=[];
 for i=1:length(x(1,:))
   u=-Kcont*(x(:,i)-xd);
   U=[U u];
 end 
scf(0)
 clf()
 plot(t,x(1,:),'r')
 plot(t,x(2,:),'b')
 plot(t,x(3,:),'m')
scf(1)
clf()
 plot(t,U(1,:),'r')
disp(' ** FIM ** ') 
