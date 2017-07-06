clc
clear
mode(-1)
rand("seed",5);
rand("normal");
q=[.03 0.01;.01 0.03];
u=rand(2,11);
f=[1.1 0.1;0 0.8];
g=(chol(q))';
m0=[10 10]';
p0=[2 0;0 2];
x0=m0+(chol(p0))'*rand(2,1);
x=ltitr(f,g,u,x0);
r=[2 0;0 2];
v=(chol(r))'*rand(2,11);
y=x+v;
h=eye(2,2);
[xe pe]=sskf(y,f,h,q,r,m0)
scf(0)
clf()
plot(xe(1,:),xe(2,:),'r-o')
plot(x(1,:),x(2,:),'b')
disp(pe)
disp('FIM')
