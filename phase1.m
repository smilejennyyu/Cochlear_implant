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
figure('Name', 'Raw Data Mono vs. Sample Number');
plot(sample_number, raw_data_mono,'g');
title('3.5 Mono Sound Waveform(Amplitude vs. Sample Number)');

% 3.6 If sampling rate is greater than 16k, downsample it
data_16k = resample(raw_data_mono, 16000, sample_rate);
%disp(strcat('data_16k size: ', string(size(data_16k))));
rate_16k = 16000;
stop_time_16k = size(data_16k)/rate_16k;
t_16 = 1/rate_16k:1/rate_16k:stop_time_16k;

%figure('Name', 'Raw Data Resampled');
%plot(t_16(1:1000), data_16k(1:1000),'r');
%sound(data_16k, rate_16k);
%figure('Name', 'Data Comparison');
%plot(t_16(1:1000), data_16k(1:1000),'r',t(1:2750), raw_data_mono(1:2750),'g');


% 3.7 Generate a signal using 1kHz cosine function that has the same time
% Play the sound and plot two cycles of its waveform versus time
freq = 1000;
cos_1kHz = cos(2000*pi*t);
sound(cos_1kHz, freq);
%disp(t(size(t)));
cycle_plot = 2;
time_plot = cycle_plot/freq;
figure('Name', 'Two Waveforms of Cosine Function');
plot(t(1:time_plot*sample_rate),(cos_1kHz(1:time_plot*sample_rate)) );
title('3.7 Two Waveforms of Cosine(Amplitude vs. Time)');

