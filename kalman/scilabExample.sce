clear;
clc;
// Construction of the sinusoid
w=%pi/4; // angular frequency
T=0.1; // period
t=0:T:100;
signal=cos(w*t);
// Sinusoid with noise
v=rand(t,"normal");
y=signal+v;
// Plot the sinusoid with noise
//subplot(2,1,1);
//plot(t,y);
//xtitle("sinusoid with noise","t");
// System
n=2; // system order
f=[cos(w*T) -sin(w*T); sin(w*T) cos(w*T)];
g=0;
h=[1 0];
p0=[1000 0; 0 0];
R=1;
Q=0;
x0=zeros(n,1);
// Initialize for loop
x1=x0;
p1=p0;

xp = x0;
xk = x0;
pp = p0;

// Kalman filter
for i=1:length(t)-1
    [x1(:,i+1),p1,x,p]=kalm(y(i),x1(:,i),p1,f,g,h,Q,R);
end

Xk = []
Xp = []
//Kalman filter according to scilab algorithm
for i=1:length(t)
    //measurement update (correction)
    y = signal(i);    
    E = y - h*xp; //innovation
    s= h*pp*h' + R;
    kk = pp*h'*inv(s);    
    xk = xk + kk*E;
    pp = pp - kk*h*pp;    
    //time update (prediction)
    xp = f*xk;
    pp = f*pp*f' + Q;
    Xk = [Xk xk];
    Xp = [Xp xp];
end

//Xk = []
////Kalman filter according to algorithm from kalman filter for dummies
//for i=1:length(t)-1
//    //signal output        
//    zk = signal(i);
//            
//    //Measurement update
//    //Kalman gain
//    kk = pp*h'*inv(h*pp*h' + R);
//    //update estimate
//    xk = xk + kk*(zk - h*xp);
//    //update error covariance
//    pp = pp - pp*kk*h;
//    Xk = [Xk xk];
//    
//    //Time update
//    //predicted state
//    xp = f*xk;
//    ///error covariance
//    pp = f*pp*f' + Q;
//    
//    
//end


//// Plot the results (in red) to compare with the sinusoid (in green)
//figure(1);
//subplot(3,1,1);
//plot(t,signal,"color","blue");
//subplot(3,1,2);
//plot(t,y,"color","blue");
//subplot(3,1,3);
//plot(t,signal,"color","blue");
//plot(t,x1(1,:),"color","red");
//xtitle("Comparison between sinusoid (green) and extraction with Kalman filter (red)","t");
