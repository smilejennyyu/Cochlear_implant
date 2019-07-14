% 3.1 Create program to read the file
filename = 'Khan_Girl.wav';
[raw_data,sample_rate] = audioread(filename);

% 3.2 Check if input sound is stereo
[m,n] = size(raw_data);
if n > 1 % If stereo:
    % Combine two channels into one and take average
    raw_data_mono = sum(raw_data,2) / size(raw_data,2);
    %disp(raw_data_mono);
else
    raw_data_mono = raw_data;
end

% 3.3 Play the sound
%sound(raw_data_mono, sample_rate);

% 3.4 Write sound to a new file 
output_file = strcat('new_', filename);
audiowrite(output_file, raw_data_mono, sample_rate);

% 3.5 Plot sound wave as a function of sample number
stop_time = size(raw_data_mono)/sample_rate;
time_step = 1/sample_rate;
t = time_step:time_step:stop_time;
sample_number = 1:1:size(raw_data_mono);
%t = 1:1:size(raw_data_mono);

%disp(strcat('raw_data_mono size: ', string(size(raw_data_mono))));
%disp(strcat('t size: ', string(size(t))));
%figure('Name', 'Raw Data Mono vs. Sample Number');
%plot(sample_number, raw_data_mono,'g');
%title('3.5 Mono Sound Waveform(Amplitude vs. Sample Number)');

% 3.6 If sampling rate is greater than 16k, downsample it
data_16k = resample(raw_data_mono, 16000, sample_rate);
%disp(strcat('data_16k size: ', string(size(data_16k))));
rate_16k = 16000;
stop_time_16k = size(data_16k)/rate_16k;
t_16 = 1/rate_16k:1/rate_16k:stop_time_16k;

% figure('Name', 'Raw Data Resampled');
% plot(t_16(1:1000), data_16k(1:1000),'r');
% sound(data_16k, rate_16k);
% figure('Name', 'Data Comparison');
% plot(t_16(1:1000), data_16k(1:1000),'r',t(1:2750), raw_data_mono(1:2750),'g');




%------Phase 2 start here-------%
%------------------------------------------------------------------------------------%
% declaring filter parameters
bandwidth=659;
pass_1=100;
pass_2= pass_1+659;

% length = length(filtered_data16k);

% figure('Name', 'Raw Data 16k');
% plot(t_16(1:1000), data_16k(1:1000),'r');
% filtered_data16k = filter(Equiripple1(pass_1, pass_2), data_16k);
% figure('Name', 'Bandpass output');
% plot(t_16(1:1000), filtered_data16k(1:1000),'r');
% 
% rec_filtered_data16k = abs(filtered_data16k);
% figure('Name', 'Rectified bandpass output');
% plot(t_16(1:1000), rec_filtered_data16k(1:1000),'r');
% hold on;
% % lp1 = dsp.FIRFilter('Numerator', firpm(20,[0 0.03 0.1 1],[1 1 0 0]));
% % env_rec_filtered_data16k = lp1(rec_filtered_data16k);
% env_rec_filtered_data16k  = filter(lp30, rec_filtered_data16k);
% plot(t_16(1:1000), env_rec_filtered_data16k(1:1000),'b');
% fvtool(lp30);
% sound(data_16k, rate_16k);


% figure('Name', 'last Raw Data 16k');
% plot(t_16(length-1000:length), data_16k(length-1000:length),'r');
% filtered_data16k = filter(Equiripple1(pass_1, pass_2), data_16k);
% figure('Name', 'last Bandpass output');
% plot(t_16(length-1000:length), filtered_data16k(length-1000:length),'r');
% 
% rec_filtered_data16k = abs(filtered_data16k);
% figure('Name', 'last Rectified bandpass output');
% plot(t_16(length-1000:length), rec_filtered_data16k(length-1000:length),'r');
% hold on;
% env_rec_filtered_data16k  = filter(lowbutter, rec_filtered_data16k);
% plot(t_16(length-1000:length), env_rec_filtered_data16k(length-1000:length),'b');




freq = 150;
cos_1kHz = cos(freq*2*pi*t_16); % lesson learned, use the same sample rate for t!
cycle_plot = 100;
time_plot = cycle_plot/freq;
figure('Name', 'Waveform of Cosine Function');
plot(t_16(1:time_plot*rate_16k),(cos_1kHz(1:time_plot*rate_16k)),'g' );
figure('Name', 'Rectified Waveforms of Cosine Function');
plot(t_16(1:time_plot*rate_16k),(abs(cos_1kHz(1:time_plot*rate_16k))),'r' );
filtered_cos1k = filter(Equiripple1(pass_1, pass_2), cos_1kHz);
figure('Name', 'Bandpass cos 1k 1 in time');
plot(t_16(1:time_plot*rate_16k),(filtered_cos1k(1:time_plot*rate_16k)) );
figure('Name', 'Rectified Bandpass Filter');
filtered_Rect_cos1k = abs(filtered_cos1k);
plot(t_16(1:time_plot*rate_16k),(filtered_Rect_cos1k(1:time_plot*rate_16k)),'r' );
hold on;
Enved_filtered_Rect_cos1k = filter(lp30,filtered_Rect_cos1k);
% figure('Name', 'Enved Rectified Bandpass Filter');
plot(t_16(1:time_plot*rate_16k),(Enved_filtered_Rect_cos1k(1:time_plot*rate_16k)),'b' );
%fvtool(GenEquiripple(pass_1,pass_2));
% fvtool(lowbutter);


function Hd30 = lp30
%LP30 Returns a discrete-time filter object.

% FIR least-squares Lowpass filter designed using the FIRLS function.

% All frequency values are in Hz.
Fs = 16000;  % Sampling Frequency

N     = 30;   % Order
Fpass = 400;  % Passband Frequency
Fstop = 450;  % Stopband Frequency
Wpass = 1;    % Passband Weight
Wstop = 1;    % Stopband Weight

% Calculate the coefficients using the FIRLS function.
b  = firls(N, [0 Fpass Fstop Fs/2]/(Fs/2), [1 1 0 0], [Wpass Wstop]);
Hd30 = dfilt.dffir(b);
end

function Hdleast = lpleastsquare
%LPLEASTSQUARE Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.4 and DSP System Toolbox 9.6.
% Generated on: 14-Jul-2019 15:40:21

% FIR least-squares Lowpass filter designed using the FIRLS function.

% All frequency values are in Hz.
Fs = 16000;  % Sampling Frequency

N     = 10;   % Order
Fpass = 400;  % Passband Frequency
Fstop = 450;  % Stopband Frequency
Wpass = 1;    % Passband Weight
Wstop = 1;    % Stopband Weight

% Calculate the coefficients using the FIRLS function.
b  = firls(N, [0 Fpass Fstop Fs/2]/(Fs/2), [1 1 0 0], [Wpass Wstop]);
Hdleast = dfilt.dffir(b);
end


function Hd = lowbutter
%LOWBUTTER Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.4 and Signal Processing Toolbox 8.0.
% Generated on: 14-Jul-2019 15:00:01

% Butterworth Lowpass filter designed using FDESIGN.LOWPASS.

% All frequency values are in Hz.
Fs = 16000;  % Sampling Frequency

Fpass = 400;         % Passband Frequency
Fstop = 450;         % Stopband Frequency
Apass = 0.1;         % Passband Ripple (dB)
Astop = 50;          % Stopband Attenuation (dB)
match = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, Fs);
Hd = design(h, 'butter', 'MatchExactly', match);
end

% [EOF]


function Hdl3 = lowpass3
%LOWPASS3 Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.4 and Signal Processing Toolbox 8.0.
% Generated on: 14-Jul-2019 14:48:09

% Equiripple Lowpass filter designed using the FIRPM function.

% All frequency values are in Hz.
Fs = 16000;  % Sampling Frequency

Fpass = 400;              % Passband Frequency
Fstop = 450;              % Stopband Frequency
Dpass = 0.0057563991496;  % Passband Ripple
Dstop = 0.0031622776602;  % Stopband Attenuation
dens  = 20;               % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});
Hdl3 = dfilt.dffir(b);
end

function Hdl2 = lowpassEqui2
%LOWPASSEQUI2 Returns a discrete-time filter object.


% Equiripple Lowpass filter designed using the FIRPM function.

% All frequency values are in Hz.
Fs = 16000;  % Sampling Frequency

Fpass = 400;               % Passband Frequency
Fstop = 450;               % Stopband Frequency
Dpass = 0.0057563991496;   % Passband Ripple
Dstop = 3.1622776602e-06;  % Stopband Attenuation
dens  = 20;                % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});
Hdl2 = dfilt.dffir(b);
end

function HdLow2 = lowpassCheb1
% bad example, use for justification
%LOWPASSCHEB1 Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.4 and DSP System Toolbox 9.6.
% Generated on: 13-Jul-2019 18:24:40

% Chebyshev Type I Lowpass filter designed using FDESIGN.LOWPASS.

% All frequency values are in Hz.
Fs = 16000;  % Sampling Frequency

N     = 100;  % Order
Fpass = 400;  % Passband Frequency
Apass = 1;    % Passband Ripple (dB)

% Construct an FDESIGN object and call its CHEBY1 method.
h  = fdesign.lowpass('N,Fp,Ap', N, Fpass, Apass, Fs);
HdLow2 = design(h, 'cheby1');
end

function HdLow = lowpassEqui
%LOWPASSEQUI Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.4 and DSP System Toolbox 9.6.
% Generated on: 13-Jul-2019 18:19:36

% Equiripple Lowpass filter designed using the FIRPM function.

% All frequency values are in Hz.
Fs = 16000;  % Sampling Frequency

Fpass = 400;             % Passband Frequency
Fstop = 450;             % Stopband Frequency
Dpass = 0.028774368332;  % Passband Ripple
Dstop = 1e-05;           % Stopband Attenuation
dens  = 20;              % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});
HdLow = dfilt.dffir(b);
end


function Hd5 = Cheb2
%CHEB2 Returns a discrete-time filter object.

% Chebyshev Type II Bandpass filter designed using FDESIGN.BANDPASS.

% All frequency values are in Hz.
Fs = 16000;  % Sampling Frequency

Fstop1 = 50;          % First Stopband Frequency
Fpass1 = 100;         % First Passband Frequency
Fpass2 = 759;         % Second Passband Frequency
Fstop2 = 809;         % Second Stopband Frequency
Astop1 = 60;          % First Stopband Attenuation (dB)
Apass  = 1;           % Passband Ripple (dB)
Astop2 = 80;          % Second Stopband Attenuation (dB)
match  = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY2 method.
h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, ...
                      Astop2, Fs);
Hd5 = design(h, 'cheby2', 'MatchExactly', match);
end

function Hd4 = Cheb1
%CHEB1 Returns a discrete-time filter object.

% Chebyshev Type I Bandpass filter designed using FDESIGN.BANDPASS.

% All frequency values are in Hz.
Fs = 16000;  % Sampling Frequency

Fstop1 = 50;          % First Stopband Frequency
Fpass1 = 100;         % First Passband Frequency
Fpass2 = 759;         % Second Passband Frequency
Fstop2 = 809;         % Second Stopband Frequency
Astop1 = 60;          % First Stopband Attenuation (dB)
Apass  = 1;           % Passband Ripple (dB)
Astop2 = 80;          % Second Stopband Attenuation (dB)
match  = 'passband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY1 method.
h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, ...
                      Astop2, Fs);
Hd4 = design(h, 'cheby1', 'MatchExactly', match);

% [EOF]
end

function Hd3 = butterworthIIR
%BUTTERWORTHIIR Returns a discrete-time filter object.

% Butterworth Bandpass filter designed using FDESIGN.BANDPASS.

% All frequency values are in Hz.
Fs = 16000;  % Sampling Frequency

Fstop1 = 50;          % First Stopband Frequency
Fpass1 = 100;         % First Passband Frequency
Fpass2 = 759;         % Second Passband Frequency
Fstop2 = 809;         % Second Stopband Frequency
Astop1 = 60;          % First Stopband Attenuation (dB)
Apass  = 1;           % Passband Ripple (dB)
Astop2 = 80;          % Second Stopband Attenuation (dB)
match  = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, ...
                      Astop2, Fs);
Hd3 = design(h, 'butter', 'MatchExactly', match);

% [EOF]
end

function Hd2 = GenEquiripple(pass1, pass2)
%GENEQUIRIPPLE Returns a discrete-time filter object.

% Generalized REMEZ FIR Bandpass filter designed using the FIRGR function.

% All frequency values are in Hz.
Fs = 16000;  % Sampling Frequency

Fstop1 = pass1-50;              % First Stopband Frequency
Fpass1 = pass1;             % First Passband Frequency
Fpass2 = pass2;             % Second Passband Frequency
Fstop2 = pass2+50;             % Second Stopband Frequency
Dstop1 = 0.001;           % First Stopband Attenuation
Dpass  = 0.028774368332;  % Passband Ripple
Dstop2 = 0.001;           % Second Stopband Attenuation
dens   = 20;              % Density Factor

% Calculate the coefficients using the FIRGR function.
b  = firgr('minorder', [0 Fstop1 Fpass1 Fpass2 Fstop2 Fs/2]/(Fs/2), [0 ...
           0 1 1 0 0], [Dstop1 Dpass Dstop2], {dens});
Hd2 = dfilt.dffir(b);
end

function Hd1 = Equiripple1(pass1, pass2)
%EQUIRIPPLE1 Returns a discrete-time filter object.

% Equiripple Bandpass filter designed using the FIRPM function.

% All frequency values are in Hz.
Fs = 16000;  % Sampling Frequency

Fstop1 = pass1-50;              % First Stopband Frequency
Fpass1 = pass1;             % First Passband Frequency
Fpass2 = pass2;             % Second Passband Frequency
Fstop2 = pass2+50;             % Second Stopband Frequency

Dstop1 = 0.001;           % First Stopband Attenuation
Dpass  = 0.028774368332;  % Passband Ripple
Dstop2 = 0.001;           % Second Stopband Attenuation
dens   = 20;              % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 1 ...
                          0], [Dstop1 Dpass Dstop2]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});
Hd1 = dfilt.dffir(b);
end


