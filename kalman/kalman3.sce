clc
clear
mode(-1)
// CPII - UFU - 2017
//test of the steady-state kalman filter
rand('seed',5);rand('normal');
q=[.03 0.01;.01 0.03];
n=20;
u=rand(2,n);
f=[1.1 0.1;0 0.8];
g=(chol(q))';
m0=[10 10]';
p0=[2 0;0 2];
x0=m0+(chol(p0))'*rand(2,1);
x=ltitr(f,g,u,x0);
r=[2 0;0 2];
v=(chol(r))'*rand(2,n);
y=x+v;
h=eye(2,2);
[xe]=sskf(y,f,h,q,r,m0);
//plot result
a=min([x(1,:),xe(1,:)]);
a=-0.1*abs(a)+a;
b=max([x(1,:),xe(1,:)]);b=.1*abs(b)+b;
c=min([x(2,:),xe(2,:)]);c=-0.1*abs(c)+c;
d=max([x(2,:),xe(2,:)]);d=.1*abs(d)+d;
//plot frame, real state (x), and estimate (xke)
scf(0)
clf()
plot([a a b],[d c c]),
plot2d(x(1,:)',x(2,:)',[1],'000',' ')
plot2d(xe(1,:)',xe(2,:)',[2],'000',' '),
plot2d(xe(1,:)',xe(2,:)',[-3],'000',' '),
scf(1)
clf()
plot(x(1,:)',x(2,:)','b'),
plot(xe(1,:)',xe(2,:)','r-o'),

disp(' END')
