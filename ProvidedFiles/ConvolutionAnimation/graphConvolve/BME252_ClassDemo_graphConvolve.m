%BME 252 convolution examples
%% Example in class:
a0=[1,2,1,0,0,0,0,0,0,0,0,0,0];
b0=[0,2,1,0];
pausespacing=0.3;
[conv,Frames] = graphConvolve(a0,b0,pausespacing);

%% Applications
%input signal
a1=ones(1,30);
a1(5:8)=0;
a1(20:24)=-1;
x=0:29;
pausespacing=0.15;
figure (100)
subplot(2,1,1)
stem(x,a1,'LineWidth',3);
subplot(2,1,2)
b1=4;
stem(x(1),b1,'LineWidth',3,'Color',[1,0,0]);
xlim([0,30])
%% multiply
b1=4;

[conv,Frames] = graphConvolve(a1,b1,pausespacing);
%
%% Delay 
b2=[0,0,0,1];
figure(100)
subplot(2,1,2)
stem(x(1:4),b2,'LineWidth',3,'Color',[1,0,0]);
xlim([0,30])
%%
[conv,Frames] = graphConvolve(a1,b2,pausespacing);
%% Delay a portion of an ECG signal by 10 samples
load ecgsig.mat
ecgseg=ecgsig(300:350);
b3=zeros(1,10);
b3(10)=1;
pausespacing=0.01;
[conv,Frames] = graphConvolve(ecgseg',b3,pausespacing);
%% Blur / smooth
b4=0.5*[.5,1,.5];
figure(100)
subplot(2,1,2)
stem(x(1:3),b4,'LineWidth',3,'Color',[1,0,0]);
xlim([0,30])
%%
pausespacing=0.15;
[conv,Frames] = graphConvolve(a1,b4,pausespacing);
%% Edge detection
b5=[1,-1];
figure(100)
subplot(2,1,2)
stem(x(1:2),b5,'LineWidth',3,'Color',[1,0,0]);
xlim([0,30])
%%
[conv,Frames] = graphConvolve(a1,b5,pausespacing);



