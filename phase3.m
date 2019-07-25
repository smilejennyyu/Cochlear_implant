%------Phase 1 start here-------%
%------------------------------------------------------------------------------------%
% 3.1 Create program to read the file
filename = 'emmaWatson.mp4';
[raw_data,sample_rate] = audioread(filename);
plot_upper = 1;
plot_lower = 5999;

% 3.2 Check if input sound is stereo
[m,n] = size(raw_data);
if n > 1 % If stereo:
    % Combine two channels into one and take average
    %       raw_data_mono = sum(raw_data,2) / size(raw_data,2);
    % TA said we should not take the average
    raw_data_mono = sum(raw_data,2);
    %disp(raw_data_mono);
else
    raw_data_mono = raw_data;
end

% 3.3 Play the sound
% sound(raw_data_mono, sample_rate);

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
% figure('Name', 'Raw Data Mono vs. Sample Number');
% plot(sample_number, raw_data_mono,'g');
% title('3.5 Mono Sound Waveform(Amplitude vs. Sample Number)');

% 3.6 If sampling rate is greater than 16k, downsample it
data_16k = resample(raw_data_mono, 16000, sample_rate);
%disp(strcat('data_16k size: ', string(size(data_16k))));
rate_16k = 16000;
stop_time_16k = size(data_16k)/rate_16k;
t_16 = 1/rate_16k:1/rate_16k:stop_time_16k;

figure('Name', 'Raw Data Resampled against time');
plot(t_16(plot_upper:plot_lower), data_16k(plot_upper:plot_lower),'r');
% sound(data_16k, rate_16k);
% figure('Name', 'Data Comparison');
% plot(t_16(1:1000), data_16k(1:1000),'r',t(1:2750), raw_data_mono(1:2750),'g');




%------Phase 2 start here-------%
%------------------------------------------------------------------------------------%
% declaring parameters
passband_num=12;
bandwidth=659;
pass_1_array=100:bandwidth:passband_num*bandwidth;
pass_2_array=100+bandwidth:bandwidth:(passband_num+1)*bandwidth;


% ------ Plot acutal sound signal by calling function ---------%
figure('Name', 'Raw Data 16k against sample number');
plot(data_16k(plot_upper:plot_lower), 'r');
% plot(t_16(1000:1999), data_16k(1000:1999),'r');
% plotSoundAllPassband(pass_1_array, pass_2_array, data_16k, t_16, passband_num);
% plotSoundTwoPassband(pass_1_array, pass_2_array, data_16k, t_16);



% ------ Plot cosine signal to check correctness ---------%
% freq =1800;
% cos_1kHz = cos(freq*2*pi*t_16); % lesson learned, use the same sample rate for t!
% cycle_plot = 70;
% time_plot = cycle_plot/freq;
% figure('Name', 'Two Waveforms of Cosine Function');
% plot(t(1:time_plot*rate_16k),(cos_1kHz(1:time_plot*rate_16k)) );
% plotCosTwoPassband(pass_1_array, pass_2_array, cos_1kHz, time_plot, rate_16k, t_16);


% ------------------------------------------%
% ------------------------------------------%
% --------- Phase 3 starts here ------------%

freq =180;
cos_center = cos(freq*2*pi*t_16);
cycle_plot = 4;
time_plot = cycle_plot/freq;
% figure('Name', 'Two Waveforms of Cosine Function vesus time');
% plot(t(1:time_plot*rate_16k),(cos_center(1:time_plot*rate_16k)) );
% figure('Name', 'Two Waveforms of Cosine Function vesus sample number');
% plot((cos_center(1:time_plot*rate_16k)) );

% plotCosAllPassband(pass_1_array, pass_2_array, time_plot, rate_16k, t_16, passband_num);
plotModSound(pass_1_array, pass_2_array, data_16k, t_16, passband_num, plot_upper, plot_lower)

function cos_center = genCosCenter (lower, upper, t_16)
    freq =sqrt(lower*upper);
    cos_center = cos(freq*2*pi*t_16);
%     cycle_plot = 70;
%     time_plot = cycle_plot/freq;
%     figure('Name', 'Two Waveforms of Cosine Function vesus time');
%     plot(t(1:time_plot*rate_16k),(cos_center(1:time_plot*rate_16k)) );
%     figure('Name', 'Two Waveforms of Cosine Function vesus sample number');
%     plot((cos_center(1:time_plot*rate_16k)) );
end


% Phase 3 Task 11 %
function plotModSound(pass_1_array, pass_2_array, data_16k, t_16, passband_num, plot_upper, plot_lower)
    for p = 11:passband_num
        if p==12
            % For last band, the end parameter need a bit of adjustment
            filtered_data16k = filter(Cheb2(pass_1_array(p),7950), data_16k);
            filtered_cos1k = filter(Cheb2(pass_1_array(p),7950), genCosCenter(pass_1_array(p), 7950, t_16));
        else
            % Use Cheb2 as the band pass filter function
            filtered_data16k = filter(Cheb2(pass_1_array(p),pass_2_array(p)), data_16k);
            filtered_cos1k = filter(Cheb2(pass_1_array(p),pass_2_array(p)),  genCosCenter(pass_1_array(p), pass_2_array(p), t_16));
        end
        rec_filtered_data16k = abs(filtered_data16k);
        env_rec_filtered_data16k = filter(lp30, rec_filtered_data16k);
        modulated_data16k = times(env_rec_filtered_data16k(plot_upper:plot_lower),filtered_cos1k(plot_upper:plot_lower));
%          modulated_data16k = env_rec_filtered_data16k .* filtered_cos1k;
%         title=strcat('Env & Rec Sound in frequency',num2str(pass_1_array(p)),'-',num2str(pass_2_array(p)));
%         figure('Name', title);
% %         Only plotted a tiny cycle, but actual signal is at size data_16k
%         plot(rec_filtered_data16k(plot_upper:plot_lower), 'r');
%         hold on;
%         plot(env_rec_filtered_data16k(plot_upper:plot_lower), 'b');
%         xlabel('Sample Number') 
%         ylabel('Amplitude')
        % second graph of amplitude modulated signal
        title=strcat('Amplitude Modulated Sound in frequency',num2str(pass_1_array(p)),'-',num2str(pass_2_array(p)));
        figure('Name', title);
        plot(modulated_data16k, 'g');
        xlabel('Sample Number') 
        ylabel('Amplitude')
    end
end


% For Task 4 and Task 8 in Phase 2 %
function plotSoundTwoPassband(pass_1_array, pass_2_array, data_16k, t_16)
    plot_upper = 1000;
    plot_lower = 1999;
    % Declaring parameters for passband parameters for first and last channels
    low_start=pass_1_array(1);
    low_end=pass_2_array(1);
    high_start=pass_1_array(length(pass_1_array));
    
    % Plot lowest channel 
    filtered_data16k_low = filter(Cheb2(low_start, low_end), data_16k);
    figure('Name', 'Bandpass output of lowest channel');
%     plot(filtered_data16k_low(plot_upper:plot_lower), 'r');
    plot(t_16((plot_upper:plot_lower)), filtered_data16k_low((plot_upper:plot_lower)),'r');
    xlabel('Time (s)'); 
    ylabel('Amplitude');
    % Plot rectified lowest channel and the envelope using low pass filter
    rec_filtered_data16k_low = abs(filtered_data16k_low);
    env_rec_filtered_data16k_low  = filter(lp30, rec_filtered_data16k_low);
    figure('Name', 'Envelope of the rectified lowest channel');
%     plot(rec_filtered_data16k_low(plot_upper:plot_lower), 'r');
    plot(t_16(plot_upper:plot_lower), rec_filtered_data16k_low(plot_upper:plot_lower),'r'); 
    hold on;
    plot(t_16(plot_upper:plot_lower), env_rec_filtered_data16k_low(plot_upper:plot_lower),'b'); 
%     plot(env_rec_filtered_data16k_low(plot_upper:plot_lower), 'b');
    xlabel('Time (s)') 
    ylabel('Amplitude')
    
    % Plot highest channel 
    filtered_data16k_high = filter(Cheb2(high_start, 7950), data_16k);
    figure('Name', 'Bandpass output of highest channel');
%     plot(filtered_data16k_high(plot_upper:plot_lower), 'r');
    plot(t_16(plot_upper:plot_lower), filtered_data16k_high(plot_upper:plot_lower),'r');
    xlabel('Time (s)') 
    ylabel('Amplitude')
    % Plot rectified highest channel and the envelope using low pass filter
    rec_filtered_data16k_high = abs(filtered_data16k_high);
    env_rec_filtered_data16k_high  = filter(lp30, rec_filtered_data16k_high);
    figure('Name', 'Envelope of the rectified highest channel');
%     plot(rec_filtered_data16k_high(plot_upper:plot_lower), 'r');
    plot(t_16(plot_upper:plot_lower), rec_filtered_data16k_high(plot_upper:plot_lower),'r');
    hold on;
%     plot(env_rec_filtered_data16k_high(plot_upper:plot_lower), 'b');
    plot(t_16(plot_upper:plot_lower), env_rec_filtered_data16k_high(plot_upper:plot_lower),'b');
    xlabel('Time (s)') 
    ylabel('Amplitude')
end

% For loop function to filter and plot each channels %
function plotSoundAllPassband(pass_1_array, pass_2_array, data_16k, t_16, passband_num)
    % Loop through each band
    for p = 1:passband_num
        if p==12
            % For last band, the end parameter need a bit of adjustment
            filtered_data16k = filter(Cheb2(pass_1_array(p),7950), data_16k);
        else
            % Use Cheb2 as the band pass filter function
            filtered_data16k = filter(Cheb2(pass_1_array(p),pass_2_array(p)), data_16k);
        end
        title=strcat('Bandpass sound file in frequency range ',num2str(pass_1_array(p)),'-',num2str(pass_2_array(p)));
        figure('Name', title);
        plot(t_16(1000:1999), filtered_data16k(1000:1999),'r');
    end
end

% Checker function using cosine %
function plotCosTwoPassband(pass_1_array, pass_2_array, cos_1kHz, time_plot, rate_16k, t_16)
    low_start=pass_1_array(1);
    low_end=pass_2_array(1);
    high_start=pass_1_array(length(pass_1_array));
    
    filtered_cos1k_low = filter(Cheb2(low_start,low_end), cos_1kHz);
    figure('Name', 'Bandpass cosine in frequency (low)');
    plot(t_16(1:time_plot*rate_16k),(filtered_cos1k_low(1:time_plot*rate_16k)) );
    filtered_cos1k_high = filter(Cheb2(high_start,7950), cos_1kHz);
    figure('Name', 'Bandpass cosine in frequency (high)');
    plot(t_16(1:time_plot*rate_16k),(filtered_cos1k_high(1:time_plot*rate_16k)));
end

% Checker function using cosine %
function plotCosAllPassband(pass_1_array, pass_2_array, time_plot, rate_16k, t_16, passband_num)
    for p = 1:passband_num
        if p==12
            
            filtered_cos1k = filter(Cheb2(pass_1_array(p),7950), genCosCenter(pass_1_array(p), 7950, t_16));
        else
            
            filtered_cos1k = filter(Cheb2(pass_1_array(p),pass_2_array(p)),  genCosCenter(pass_1_array(p), pass_2_array(p), t_16));
        end
        title=strcat('Bandpass cosine function in frequency',num2str(pass_1_array(p)),'-',num2str(pass_2_array(p)));
        figure('Name', title);
        % Only plotted a tiny cycle, but actual signal is at size data_16k
        plot(t_16(1:time_plot*rate_16k),(filtered_cos1k(1:time_plot*rate_16k)) );
        
    end
end

%----- low pass filter bank ----%
% lp30 is the chosen low pass filter
function Hd30 = lp30
%LP30 Returns a discrete-time filter object.

% FIR least-squares Lowpass filter designed using the FIRLS function.

% All frequency values are in Hz.
Fs = 16000;  % Sampling Frequency

N     = 50;   % Order
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


%----- band pass filter bank ----%
% Cheb2 is the chosen band pass filter
function HdC = Cheb2(pass1,pass2)
%CHEB2ASTOP100 Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.4 and DSP System Toolbox 9.6.
% Generated on: 14-Jul-2019 19:37:31

% Chebyshev Type II Bandpass filter designed using FDESIGN.BANDPASS.

% All frequency values are in Hz.
Fs = 16000;  % Sampling Frequency
Fstop1 = pass1-50;              % First Stopband Frequency
Fpass1 = pass1;             % First Passband Frequency
Fpass2 = pass2;             % Second Passband Frequency
Fstop2 = pass2+50;             % Second Stopband Frequency
Astop1 = 100;         % First Stopband Attenuation (dB)
Apass  = 1;           % Passband Ripple (dB)
Astop2 = 100;         % Second Stopband Attenuation (dB)
match  = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY2 method.
h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, ...
                      Astop2, Fs);
HdC = design(h, 'cheby2', 'MatchExactly', match);
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
display(pass1);
display(pass2);
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

%------Phase 2 end here-------%


