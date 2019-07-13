% BME 252 DT frequency response and filter demo
% Nima Maftoom July-032019

close all

%% Example in the course note (July 3)
w=-pi:0.1:pi;
TFZ= @(z) (1)./(z.^3);
HZ1=TFZ(exp(1j*w));
figure
subplot(2,1,1)
% loglog(w,abs(HZ1));
plot(w,abs(HZ1));
subplot(2,1,2)
% semilogx(w,unwrap(angle(HZ1)));
plot(w,unwrap(angle(HZ1)));

%% Example in the course note (July 8)
w=1*(-pi:0.001:pi);
TFZ= @(z) (z.^2+1)./(z.^2);
HZ1=TFZ(exp(1j*w));
figure
subplot(2,1,1)
% loglog(w,abs(HZ1));
plot(w,abs(HZ1));
subplot(2,1,2)
% semilogx(w,unwrap(angle(HZ1)));
plot(w,unwrap(angle(HZ1)));
%% Example A in the course note (July 3)
% solution of the difference equation
% y[n]=0.5y[n-1]-0.1y[n-2]+0.1x[n]-0.5x[n-1]+x[n-2]

% part (a)
%response to x[n]=[1,2,-1],n=0,1,2
 N=10;
 a=[1 -0.5 0.1];
 b=[0.1 -0.5 1];
 x=[1 2 -1 zeros(1,7)];
 y=filter(b,a,x);
 stem(0:N-1,y);
 legend('y[n]');
 
% part (b) 
%impulse response
figure
N=10;
a=[1 -0.5 0.1];
b=[0.1 -0.5 1];
d=[1 zeros(1,N-1)];
h=filter(b,a,d);
stem(0:N-1,h);
xlim([-.5 10]);
legend('h[n]');

% part (c)
w=1*(-pi:0.001:pi);
TFZ= @(z) (1-0.5*z+0.1*z.^2)./(0.1-0.5*z+z.^2);
HZ1=TFZ(exp(1j*w));
figure
subplot(2,1,1)
% loglog(w,abs(HZ1));
plot(w,abs(HZ1));
subplot(2,1,2)
% semilogx(w,unwrap(angle(HZ1)));
plot(w,unwrap(angle(HZ1)));

%% Example B in the course note (July 3)
% solution of the difference equation
% y[n]-1.1y[n-1]+0.9y[n-2]=x[n]

% part (a)
%response to x[n]=[1,2,-1],n=0,1,2
a=[1 -1.1 0.9];
b=1;
n=0:100;
x1=[1 2 -1];
n2=3:100;
x2=zeros(size(n2));
x=[x1 x2];
y=filter(b,a,x);
stem(n,y);
axis([ -2 102 -3 3.5]);
legend('y[n]');

% part (b)
%impulse response
a=[1 -1.1 0.9];
b=1;
n=0:100;
x=[1,zeros(1,100)]
h=filter(b,a,x);
stem(n,h);
axis([ -2 102 -1 1.2]);
legend('h[n]');
