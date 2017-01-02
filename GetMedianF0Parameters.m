function output = GetMedianF0Parameters(probability_structure)
%% produce the most likely median F0
%
% output = GetMedianF0Parameters(probability_structure)
%
% Input argument
%   probability_structure : probability structure
%
% Return value
%   output: median F0 candidate

% Copyright 2016 Google Inc. All Rights Reserved
% Author: hidekik@google.com (Hideki Kawahara)

prob_map = probability_structure.observation_probability_map;
amp_map = probability_structure.amplitude_map;
fc_list = probability_structure.center_frequency_list;
channels_in_octave = probability_structure.channels_in_octave;

[~, ac_chanel] = min(abs(fc_list - 60));
rms_dB = 20 * log10(sum(amp_map(ac_chanel:end, :)));
amount = 0.6;
boundary_level = GetDefaultVuvThreshold(rms_dB, amount);
reliable_cum_prob = cumsum(mean(prob_map(:, rms_dB > boundary_level), 2));
global_jump = reliable_cum_prob(channels_in_octave:end) - ...
  reliable_cum_prob(1:end - channels_in_octave + 1);
[~, max_index] = max(global_jump);
median_index = round(channels_in_octave / 2 + max_index);
median_f0 = fc_list(median_index);
output = median_f0;
end

