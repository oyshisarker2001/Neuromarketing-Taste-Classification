clearvars; close all;

% Load your .mat file
myfile = 'Reilly_S1_05.mat';
data = load(myfile);  % loads struct

% Assume variable 'y' contains both fNIRS and EEG channels
y = data.y;

% Sampling frequency and parameters
Fs = 500;                 % Sampling frequency in Hz
startTime = 26;           % Trial start time (seconds)
endTime = 36;             % Trial end time (seconds)
window_sec = 2;           % Window length for segmentation (seconds)
step_sec = 1;             % Step size for segmentation (seconds)

% Ensure 'y' is oriented as (time_points x channels)
if size(y,1) < size(y,2)
    y = y';  % transpose if channels are rows instead of columns
end


% Convert start and end times to sample indices
startPoint = floor(startTime * Fs) + 1;
endPoint = floor(endTime * Fs) + 1;

% Clip endPoint if it exceeds available data length
maxTimePoints = size(y, 1);
if endPoint > maxTimePoints
    warning('Adjusting endPoint from %d to %d to fit data length', endPoint, maxTimePoints);
    endPoint = maxTimePoints;
end

% Extract fNIRS channels (channels 2 to 17)
% Ensure that channels exist; adapt if fewer channels in y
numChannelsTotal = size(y,2);
if numChannelsTotal < 17
    error('Data has only %d channels, cannot extract channels 2 to 17', numChannelsTotal);
end
fNIRS_data = y(startPoint:endPoint, 2:17);

% Optionally extract EEG channels (channels 18 to 29) if they exist
EEG_data = [];
if numChannelsTotal >= 29
    EEG_data = y(startPoint:endPoint, 18:29);
end

% Check data size after extraction
num_samples = size(fNIRS_data, 1);
num_channels = size(fNIRS_data, 2);
fprintf('Extracted fNIRS data size: %d samples x %d channels\n', num_samples, num_channels);

% ----- Preprocessing -----
% 1. Bandpass filter design (Butterworth)
order = 3;                 % Filter order
low_cutoff = 0.05;         % Low cutoff frequency in Hz
high_cutoff = 0.2;         % High cutoff frequency in Hz
Wn = [low_cutoff high_cutoff] / (Fs/2);   % Normalized cutoff frequencies

[b, a] = butter(order, Wn, 'bandpass');

% 2. Apply zero-phase filtering
fNIRS_filtered = filtfilt(b, a, fNIRS_data);

% 3. Baseline correction (zero mean per channel)
fNIRS_detrended = fNIRS_filtered - mean(fNIRS_filtered, 1);

% ----- Segmenting into windows -----
window_size = window_sec * Fs;  % samples per window
step_size = step_sec * Fs;      % step size in samples

% Compute number of segments
segment_count = floor((num_samples - window_size) / step_size) + 1;

if segment_count < 1
    error(['Not enough data samples (%d) for one window of size %d. ', ...
           'Adjust start/end times, window_sec, or collect more data.'], num_samples, window_size);
end

% Pre-allocate segments array: (window_size x num_channels x segment_count)
segments = zeros(window_size, num_channels, segment_count);

for i = 1:segment_count
    start_idx = (i-1)*step_size + 1;
    end_idx = start_idx + window_size - 1;
    segments(:, :, i) = fNIRS_detrended(start_idx:end_idx, :);
end

% Label all segments with the taste label (change as appropriate)
taste_label = 'Umami';    % example label
labels = repmat({taste_label}, segment_count, 1);

% Permute segments to (segments x time_steps x channels)
segments_permuted = permute(segments, [3, 1, 2]);

% Save processed segments, labels, and Fs to .mat file
save('Processed_fNIRS_Segments.mat', 'segments_permuted', 'labels', 'Fs');

% ----- Optional: Power spectral density (PSD) calculation for first segment/channel -----
signal = squeeze(segments(1, :, 1));  % first segment, all time points, first channel
window = [];
nfft = length(signal);
[P, F] = periodogram(signal, window, nfft, Fs);

% Calculate band power in common EEG bands (example)
theta_power = bandpower(P, F, [4 8], 'psd');
alpha_power = bandpower(P, F, [8 13], 'psd');
beta_power  = bandpower(P, F, [13 30], 'psd');


fprintf('Theta power in first segment: %.4f\n', theta_power);
fprintf('Alpha power in first segment: %.4f\n', alpha_power);
fprintf('Beta power in first segment: %.4f\n', beta_power);
