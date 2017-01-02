function output = ...
  EstimateInitialF0AndDeviation(probability_structure, trajectory, option)
% estimate F0 initial value from instantaneous frequency map, residual map
% and given trajectory
%
% output = ...
%   EstimateInitialF0AndDeviation(probability_structure, trajectory)
%
% output = ...
%   EstimateInitialF0AndDeviation(probability_structure, trajectory, option)
%
% Arguments
%
%   probability_structure: structure variable, following fields are used
%     inst_freq_map: matrix, time-frequency map of instantaneous frequency
%     residual_map: matrix, time-frequency map of residuals
%     center_frequency_list: vector, list of filter center frequency
%     channels_in_octave: scalar, number of channels in one octave
%  trajectory: vector: sequence of F0 trajectory in channel ID number
%  option: structure variable with the following fields
%    residual_margin: scalar, safeguard for residual selection
%    frequency_margin: scalar, safeguard for limiting instantaneous
%      frequency
%
% Return value
%
%  output: structure variable with the following fields
%    f0_initial: vector, initial estimate of F0 sequence (Hz)
%    sd_initial: vector, initial estimate of relative deviation
%    temporal_positions: vector, frame center location (s)
%    residual_margin: scalar, actual value used
%    frequency_margin: scalar, actual value used
%    sd_on_trajectory: vector, residual value on give trajectory
%    elapsed_time: scalar, time needed to process this function

% Copyright 2016 Google Inc. All Rights Reserved
% Author: hidekik@google.com (Hideki Kawahara)

start_tic = tic;
residual_margin = 1.5; % default value
frequency_margin = 1.15; % default value
narginchk(2, 3);
if nargin == 3
  if isfield(option, 'residual_margin')
    residual_margin = option.residual_margin;
  end;
  if isfield(option, 'frequency_margin')
    frequency_margin = option.frequency_margin;
  end;
end;
inst_freq_map = probability_structure.inst_freq_map;
residual_map = probability_structure.residual_map;
fc_list = probability_structure.center_frequency_list;
n_channels = size(inst_freq_map,1);
n_frames = size(inst_freq_map, 2);
f0_initial = zeros(n_frames, 1);
sd_initial = zeros(n_frames, 1);
sample_f0_floor = ...
  probability_structure.center_frequency_list(min(trajectory));
sample_f0_ceiling = ...
  probability_structure.center_frequency_list(max(trajectory));
sd_on_trajectory = zeros(n_frames, 1); % for debug
n_in_octave = probability_structure.channels_in_octave;
scan_base = -ceil(n_in_octave / 2 + 2):ceil(n_in_octave / 4);
for ii = 1:n_frames
  base_index = trajectory(ii) + scan_base;
  base_index = base_index(base_index > 0 & base_index <= n_channels);
  reference_residual = residual_map(trajectory(ii), ii);
  sd_on_trajectory(ii) = sqrt(reference_residual);
  base_index = base_index(residual_map(base_index, ii) < ...
    reference_residual * residual_margin);
  reference_if = fc_list(trajectory(ii));
  tmp_residual = residual_map(base_index, ii);
  tmp_if = inst_freq_map(base_index, ii);
  tmp_if = max(reference_if / frequency_margin, ...
    min(reference_if * frequency_margin, tmp_if));
  if length(base_index) > 1
    best_weight = GetBestMixingWeight(sqrt(tmp_residual));
  else
    best_weight = 1;
  end;
  f0_initial(ii) = sum(best_weight .* tmp_if);
  sd_initial(ii) = sqrt(sum(best_weight.^2 .* tmp_residual));
end;
output = ...
  struct('f0_initial', max(sample_f0_floor, ...
                           min(sample_f0_ceiling, f0_initial)), ...
         'sd_initial', sd_initial, ...
         'temporal_positions', probability_structure.temporal_positions, ...
         'sample_f0_floor', sample_f0_floor, ...
         'sample_f0_ceiling', sample_f0_ceiling, ...
         'residual_margin', residual_margin, ...
         'frequency_margin', frequency_margin, ...
         'sd_on_trajectory', sd_on_trajectory, ...
         'elapsed_time', toc(start_tic));
end
