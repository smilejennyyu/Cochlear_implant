% BME 252 CT frequency response and filter demo
% Nima Maftoom June-062018

close all
%% A low-pass filter with transfer function H1(S)=1/(S+1)
w=0:.01:20;
TF=@(s) (1)./(s+1); % Transfer function
% H=1./(1j.*w+1);
H11=TF(1j*w);
% plot the frequency response linear scale
figure%1
subplot(2,1,1) %plot in 2x1 grid at position 1
plot(w,abs(H11));

subplot(2,1,2) %plot in 2x1 grid at position 2
plot(w,angle(H11));
% plot the frequency response log scale
figure %2
subplot(2,1,1)
loglog(w,abs(H11)); 

subplot(2,1,2)
semilogx(w,angle(H11));

%% Same function method #2
numerator=1; %b
denominator= [1,1]; %a
H12=freqs(numerator,denominator,w); %laplace transform (b,a,w)
figure%3
subplot(2,1,1)
plot(w,abs(H12));

subplot(2,1,2)
plot(w,angle(H12));
%% Same filter method #3
H13=tf(numerator,denominator)
figure %4
% plot the frequency response 
bode (H13); %beautifully labeled and correct axis!

%% Another low-pass filter with transfer function H2(S)=2/(4S+2000)
% Look at where The numerator and denominator coefficients go. 
% Look at the matlab prompt also: 
% The transfer function is reported there as well:
H2=tf(2,[4,2000])
figure %5
bode (H2);

figure%7
impulseplot(H2)


%%
% A low-pass filter with transfer function H1(S)=wc/(S+wc)
wc=2*pi*200; % cutoff-frequency 200 Hz
H3=tf(wc,[1,wc])

bode(H3)

figure
impulseplot(H3)
%%
% Defining transfer function with pole, zero and gain
% The system has these poles: lambda1=-0.2+1j and its conjugate lambda2=-0.2-1j
% no zeros and gain is 1
% H4(s)=1/{(s-lambda1)(s-lambda2)}
lambda1=-0.2+1j;
lambda2=conj(lambda1);
sys1=zpk([],[lambda1,lambda2],1);
% converting to rational polynomial form (look at matlab prompt)
H4=tf(sys1)
figure
bode(sys1);

figure
impulseplot(H4)
%%
% the above filter could be defined like the following also
sys1check=tf(1,[1,0.4,1.04])
sys1checkZP=zpk(sys1check);
figure
bode(sys1check);

%% Butterworth Lowpass Filter
N=5;
wc=500;
[num,den]=butter(N,wc,'s');
LPF=tf(num,den)
figure
bode(LPF)

%% Butterworth bandpass Filter
N=5;
[num,den]=butter(N,[300,700],'s');
BPF=tf(num,den)
figure
bode(BPF)


