function frame_is_unvoiced = ...
  DecideVUV(x, fs, probability_structure, median_f0, vuv_setting)
% voice activation detector
%
% frame_is_unvoiced = ...
%   DecideVUV(x, fs, probability_structure, best_channel_out)
%
% frame_is_unvoiced = ...
%   DecideVUV(x, fs, probability_structure, best_channel_out, vuv_setting)
%
%  Input argumrnt
%    x : speech sample data, vector
%    fs : sampling frequency (Hz)
%    probability_structure : structure variable, following fields are used
%             temporal_positions: analysis frame center position (s)
%                    frame_shift: frame shiift (s)
%                  inst_freq_map: instantaneous frequency map (Hz)
%                   residual_map: residual map of front end (dB)
%    observation_probability_map: probability map
%                  amplitude_map: amplitude map (rms)
%             channels_in_octave: detector channels per octave
%          center_frequency_list: list of filter center frequencies (Hz)
%    median_f0 : median F0
%    vuv_setting : structure variable, see GenerateVUVSetting.m for detail
%
%  Return value
%    output : structure variable
%      unvoiced_frame_refined : unvoiced mask: 1: unviced 0: voiced
%      unvoiced_frame_w : internal use. amplitude based mask
%      vuv_setting : parameter set, used for decision

% Copyright 2016 Google Inc. All Rights Reserved
% Author: hidekik@google.com (Hideki Kawahara)

%TODO(hidekik) This function is too simple. Huge room for refinement

narginchk(4, 5);
if nargin == 4
  vuv_setting = GenerateVUVSetting;
end;
n_harmonics = vuv_setting.n_harmonics;
smooth_length = vuv_setting.smooth_length; % 30 ms Hanning
slient_threshold_vuv = vuv_setting.slient_threshold_vuv; % dB
lower_level_threshold = vuv_setting.lower_level_threshold; % dB
low_level_bias = vuv_setting.low_level_bias;
unvoiced_threshold_vuv = vuv_setting.unvoiced_threshold_vuv; % dB
%%
total_chennels = size(probability_structure.observation_probability_map,1);
frame_time = probability_structure.temporal_positions;
n_frames = length(frame_time);
[~, median_index] = ...
  min(abs(probability_structure.center_frequency_list - median_f0));
selected_channels = max(1, min(total_chennels, median_index + ...
  (-probability_structure.channels_in_octave:log2(n_harmonics) * ...
  probability_structure.channels_in_octave)));

%%
low_power_w = sum(probability_structure.amplitude_map(selected_channels, :) .^ 2);
half_length = max(1,round(smooth_length / probability_structure.frame_shift /2));
w = hanning(2 * half_length + 1);
w = w / sum(w);
smoothed_low_power = fftfilt(w, [low_power_w(:);zeros(length(w) * 2,1)]);
smoothed_low_power_dB = ...
  ConvertPowerToDecibel(smoothed_low_power(half_length + (1:n_frames)));
unvoiced_frame_w = ...
  double(smoothed_low_power_dB < max(smoothed_low_power_dB) + ...
         slient_threshold_vuv);
%% simple spectrogram
half_wh_length = round(smooth_length * fs / 2);
wh = hanning(half_wh_length * 2 + 1);
wh = wh / sum(wh);
fftl = 2 ^ ceil(log2(length(wh)));
fx = (0:fftl - 1) / fftl * fs;
sgram = zeros(fftl, n_frames);
data_length = length(x);
window_index = -half_wh_length:half_wh_length;
for ii = 1:n_frames
  index_current = round(frame_time(ii) * fs) + 1;
  segment_index = max(1, min(data_length, index_current + window_index));
  sgram(:, ii) = (abs(fft(x(segment_index) .* wh, fftl))) .^ 2;
end;
power_upper2k = sum(sgram(fx > 2000 & fx < fs / 2, :), 1);
power_lower2k = sum(sgram(fx < 2000 & fx > median_f0 / 2, :), 1);
power_upper2k_dB = ConvertPowerToDecibel(power_upper2k(:));
power_lower2k_dB = ConvertPowerToDecibel(power_lower2k(:));
max_power_low = ConvertPowerToDecibel(max(power_lower2k(:)));
%% unvoiced frame calculation
unvoiced_frame_w_plus = ...
  (power_upper2k_dB > power_lower2k_dB + low_level_bias & ...
  smoothed_low_power_dB < max(smoothed_low_power_dB) + ...
  unvoiced_threshold_vuv) | ...
  power_lower2k_dB < max_power_low + lower_level_threshold;
frame_is_unvoiced = ...
  double(unvoiced_frame_w(:) + unvoiced_frame_w_plus(:) > 0.5);
end
