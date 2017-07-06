clc
clear
mode(-1)
//CPII UFU 2017
//define macro which traces an ellipse
function  []=ellipse(m1,m2,s1,s2,s12)
 t=0:.1:.1+%pi*2;
 c=2*cos(t);...
 s=2*sin(t);...
 rho=s12/sqrt(s1*s2);...
 cr=sqrt(s1)*c+m1*ones(c);...
 sr=sqrt(s2)*(rho*c+sqrt(1-rho*rho)*s)+m2*ones(s);...
 plot(cr',sr')
endfunction 
//generate test process to be sent to kalman filter
//initialize state statistics (mean and err. variance)
m0=[10 10]';p0=[2 0;0 2];
//create system
f=[1.1 0.1;0 0.8];g=[1 0;0 1];h=[1 0;0 1];
//noise statistics
q=[.03 0.01;.01 0.03];r=2*eye(2,2);
//initialize system process
rand('seed',2);rand('normal');
p0c=chol(p0);
x0=m0+p0c'*rand(ones(m0));
yt=[];
//initialize kalman filter
xke0=m0;pk0=p0;
//initialize plotted variables
x=x0;xke=m0;
ell=[pk0(1,1) pk0(2,2) pk0(1,2)]';
//loop
n=20;
for k=1:n,
 //generate the state and observation at time k (i.e. x(k+1) and y(k))
 [x1,y]=system(x0,f,g,h,q,r);
 x=[x x1];
 yt=[yt y];
 x0=x1;
 //track the state with the standard kalman filter
 [xke1,pk1,xd,pd]=kalm(y,xke0,pk0,f,g,h,q,r);
 xke=[xke xke1];
 ell=[ell [pk1(1,1) pk1(2,2) pk1(1,2)]'];
 xke0=xke1;
 pk0=pk1;
//end loop
end,
//plot result
  a=min([x(1,:)-2*sqrt(ell(1,:)),xke(1,:)]);a=-0.1*abs(a)+a;
  b=max([x(1,:)+2*sqrt(ell(1,:)),xke(1,:)]);b=.1*abs(b)+b;
  c=min([x(2,:)-2*sqrt(ell(2,:)),xke(2,:)]);c=-0.1*abs(c)+c;
  d=max([x(2,:)+2*sqrt(ell(2,:)),xke(2,:)]);d=.1*abs(d)+d;
//plot frame, real state (x), and estimate (xke)
 scf(0)
 clf()
 plot([a a b],[d c c]),
 plot2d(x(1,:)',x(2,:)',[2],"000"),
 plot2d(xke(1,:)',xke(2,:)',[1],"000"),
//plot ellipses of constant likelihood (2 standard dev’s)
for k=1:n+1,
  ellipse(x(1,k),x(2,k),ell(1,k),ell(2,k),ell(3,k)),
  end,
//mark data points (* for real data, o for estimates)
plot2d(x(1,:)',x(2,:)',[-2],"000"),
plot2d(xke(1,:)',xke(2,:)',[-3],"000")
disp(' END')
