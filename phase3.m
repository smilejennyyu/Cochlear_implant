%%------Phase 1 start here-------------------------------------------------------%%
% 3.1 Create program to read the file
filename = 'Bass_Sample.wav';
[raw_data,sample_rate] = audioread(strcat('soundFiles/',filename));
plot_upper = 1;
plot_lower = 5999;

% 3.2 Check if input sound is stereo
[m,n] = size(raw_data);
if n > 1 % If stereo:
    % Combine two channels into one and take average, TA said we should not take the average
    raw_data_mono = sum(raw_data,2);
else
    raw_data_mono = raw_data;
end

% 3.3 Play the sound
% sound(raw_data_mono, sample_rate);

% 3.4 Write sound to a new file 
output_file = strcat('new_', filename);
audiowrite(strcat('outputSound/',output_file), raw_data_mono, sample_rate);

% 3.5 Plot sound wave as a function of sample number
% figure('Name', 'Raw Data Mono against sample number');
% plot(raw_data_mono, 'r');

% 3.6 If sampling rate is greater than 16k, downsample it
data_16k = resample(raw_data_mono, 16000, sample_rate);
rate_16k = 16000;
stop_time_16k = size(data_16k)/rate_16k;
t_16 = 1/rate_16k:1/rate_16k:stop_time_16k;
% figure('Name', 'Raw Data Resampled against time');
% plot(t_16(plot_upper:plot_lower), data_16k(plot_upper:plot_lower),'r');
% Uncomment line below to hear downsampled sound file
% sound(data_16k, rate_16k);


%% ------Phase 2 start here------- %%
%------------------------------------------------------------------------------------%
% declaring parameters
passband_num=18;
bandwidth=329;
% pass_1_array=100:bandwidth:passband_num*bandwidth;
% pass_2_array=100+bandwidth:bandwidth:(passband_num+1)*bandwidth;

% 18 exponetial custom tight (good)
% pass_1_array=[100 166.4 246.2 342.1 457.3 595.8 762.2 969.2 1202.5 1491.3 1838.4 2255.6 2757 3359.6 4083.8 4954.1 6000.1 7256.9]
% pass_2_array=[166.4 246.2 342.1 457.3 595.8 762.2 969.2 1202.5 1491.3 1838.4 2255.6 2757 3359.6 4083.8 4954.1 6000.1 7256.9 8000]

% 18 exp 10 gap
% pass_1_array=[100 176.4 266.2 372.1 497.3 645.8 822.2 1032.2 1282.5 1581.3 1938.4 2365.6 2877 3489.6 4223.8 5104.1 6160.1 7426.9]
% pass_2_array=[166.4 256.2 362.1 487.3 635.8 812.2 1022.2 1272.5 1571.3 1928.4 2355.6 2867 3479.6 4213.8 5094.1 6150.1 7416.9 8000]

% 18 exp 15 gap
pass_1_array = [100 181.4 276.2 387.1 517.3 670.8 852.2 1067.2 1322.5 1626.3 1988.4 2420.6 2937 3554.6 4293.8 5179.1 6240.1 7511.899]
pass_2_array = [166.4 261.2 372.1 502.3 655.8 837.2 1052.2 1307.5 1611.3 1973.4 2405.6 2922 3539.6 4278.8 5164.1 6225.1 7496.899 8000]

% ------ Plot acutal sound signal by calling function ---------%
figure('Name', 'Raw Data 16k against sample number');
plot(data_16k, 'r');
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



%% --------- Phase 3 starts here ------------%%

% Parameters %
sum_signal = zeros(size(data_16k));

% TODO: Change file!!!!! because of faint sound
freq =180;
cos_center = cos(freq*2*pi*t_16);
cycle_plot = 4;
time_plot = cycle_plot/freq;


% plotCosAllPassband(pass_1_array, pass_2_array, time_plot, rate_16k, t_16, passband_num);
% Task 11 
mod_sound = ModulateSound(pass_1_array, pass_2_array, data_16k, t_16, passband_num, plot_upper, plot_lower, sum_signal);
% Below to hear re-constructed sound file
% sound(data_16k, rate_16k);
sound(mod_sound, rate_16k);
% Task 13 Write sound to a new file 
output_file = strcat('new_', filename);
mod_output_file = strcat('new_output_', filename); 
audiowrite(strcat('finalOutput/Final_24_eq_',mod_output_file), mod_sound, rate_16k); % Does not work right now

%Take standard deviation of data16 k
std_data_16k = std(data_16k);
std_mod_sound = std(mod_sound);
sum_diff = mean(abs(mod_sound-data_16k));
square_diff = sqrt(mean((mod_sound-data_16k).^2));
average_diff = abs(mean(mod_sound)-mean(data_16k));

%% ------- Functions to call and simplify are below ------- %%
% Function to generate cos function at central freq of the given band
function cos_center = genCosCenter (lower, upper, t_16)
    freq = (lower+upper)/2;
%     freq =sqrt(lower*upper);
    cos_center = cos(freq*2*pi*t_16);
    % Uncomment below to see the plots of cos func
%     cycle_plot = 70;
%     time_plot = cycle_plot/freq;
%     figure('Name', 'Two Waveforms of Cosine Function vesus time');
%     plot(t(1:time_plot*rate_16k),(cos_center(1:time_plot*rate_16k)) );
%     figure('Name', 'Two Waveforms of Cosine Function vesus sample number');
%     plot((cos_center(1:time_plot*rate_16k)) );
end

% Phase 3 Task 11 %
function sum_signal = ModulateSound(pass_1_array, pass_2_array, data_16k, t_16, passband_num, plot_upper, plot_lower, sum_signal, evn_average_diff)
    for p = 1:passband_num
        if p==passband_num
            % For last band, the end parameter need a bit of adjustment
            filtered_data16k = filter(Cheb2(pass_1_array(p),7990), data_16k);
            filtered_cos1k = filter(Cheb2(pass_1_array(p),7990), genCosCenter(pass_1_array(p), 7990, t_16));
        else
            % Use Cheb2 as the band pass filter function
            filtered_data16k = filter(Cheb2(pass_1_array(p),pass_2_array(p)), data_16k);
            filtered_cos1k = filter(Cheb2(pass_1_array(p),pass_2_array(p)),  genCosCenter(pass_1_array(p), pass_2_array(p), t_16));
        end
        rec_filtered_data16k = abs(filtered_data16k);
        env_rec_filtered_data16k = filter(lp30, rec_filtered_data16k);
        modulated_data16k = (env_rec_filtered_data16k).*transpose(filtered_cos1k);
        env_diff = abs(mean(rec_filtered_data16k)-mean(env_rec_filtered_data16k));
        evn_average_diff(p) = env_diff;
%         title=strcat('Env & Rec Sound in frequency',num2str(pass_1_array(p)),'-',num2str(pass_2_array(p)));
%         figure('Name', title);
% % %         Only plotted a tiny cycle, but actual signal is at size data_16k
%         plot(rec_filtered_data16k(plot_upper:plot_lower), 'r');
%         hold on;
%         plot(env_rec_filtered_data16k(plot_upper:plot_lower), 'b');
%         xlabel('Sample Number') 
%         ylabel('Amplitude')
        % second graph of amplitude modulated signal
% %         title=strcat('Amplitude Modulated Sound in frequency',num2str(pass_1_array(p)),'-',num2str(pass_2_array(p)));
% %         figure('Name', title);
% %         plot(modulated_data16k, 'g');
% %         xlabel('Sample Number');
% %         ylabel('Amplitude');
        sum_signal = sum_signal + modulated_data16k;
    end
%     display(mean(evn_average_diff));
%     figure('Name', 'Amplitude Modulated Sound Sum');
%     plot(sum_signal, 'g');
%     xlabel('Sample Number');
%     ylabel('Amplitude');
        % normalized signal
        sum_signal = normalize(sum_signal,'norm',Inf);
% 
%     figure('Name', 'Normalized Amplitude Modulated Sound Sum');
%     plot(sum_signal, 'b');
%     xlabel('Normalized Sample Number');
%     ylabel('Amplitude');
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

%% ------- Filter banks below ------- %%
%----- low pass filter bank ----%
% lp30 is the chosen low pass filter
function Hd30 = lp30
%LP30 Returns a discrete-time filter object.

% FIR least-squares Lowpass filter designed using the FIRLS function.

% All frequency values are in Hz.
Fs = 16000;  % Sampling Frequency

N     = 50;   % Order
Fpass = 470;  % Passband Frequency
Fstop = 500;  % Stopband Frequency
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
Fstop1 = pass1-10;              % First Stopband Frequency
Fpass1 = pass1;             % First Passband Frequency
Fpass2 = pass2;             % Second Passband Frequency
Fstop2 = pass2+10;             % Second Stopband Frequency
Astop1 = 80;         % First Stopband Attenuation (dB)
Apass  = 1;           % Passband Ripple (dB)
Astop2 = 80;         % Second Stopband Attenuation (dB)
match  = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY2 method.
h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, ...
                      Astop2, Fs);
HdC = design(h, 'cheby2', 'MatchExactly', match);
end
