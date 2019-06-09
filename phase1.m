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
sound(raw_data_mono, sample_rate);

% 3.4 Write sound to a new file 
output_file = strcat('new_', filename);
audiowrite(output_file, raw_data_mono, sample_rate);

% 3.5 Plot sound wave as a function of sample number
stop_time = size(raw_data_mono)/sample_rate;
time_step = 1/sample_rate;
%sound(raw_data_mono, sample_rate);
t = 1:time_step:stop_time;
%t = 1:1:size(raw_data_mono);

disp(size(raw_data_mono));
disp(size(t));
figure('Name', 'Raw Data Mono');
plot(t, raw_data_mono);

% 3.6 If sampling rate is greater than 16k, downsample it
data_16k = resample(raw_data_mono, 16000, sample_rate);
rate_16k = 16000;
t_16 = 1:1/rate_16k:size(data_16k)/rate_16k;
%figure('Name', 'Raw Data resampled');
plot(t_16, data_16k,'r',t, raw_data_mono,'g');



%generate cos wave signal
figure('Name', 'cos');

%x =1:0.0001:1000;
%sound(cos(20*pi*x),sample_rate);
%disp(t)
%disp(cos(2000*t))
plot(t(1:3),cos(2000*pi*t(1:3)));
%plot(x(1:1000),cos(2*x(1:1000)));
%plot(x, cos(2*pi*x/0.001), 'b');