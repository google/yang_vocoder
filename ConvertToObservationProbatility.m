function output = ...
  ConvertToObservationProbatility(f0_structure, frame_shift)
% Convert wavelet filter bank output to ovserbation probability map
%
%   output = ConvertToObservationProbatility(f0_structure, frame_shift)
%
% Arguments
%
% f0_structure: structure with the following fields.
%     Note that only used fieds are listed below.
%   temporal_positions: vector: center location of the analysis window.
%     represented in terms of second
%   sampling_frequency: scalar: sampling frequency (Hz)
%   center_frequency_list: vector: center frequency of each wavelet filter
%   (Hz)
%   inst_freq_map: matrix: instantaneous frequency of each filter output
%   (Hz)
%   residual_map: matrix: relative level of random component (power ratio)
%   amplitude_map: matrix: amplitude of filter output (absolute value)
% frame_shift: optional variable, scalar, real value in second
%   Interval between adjacent frame center positions
%
% Return value
%
% output: structure with the following fields.
%   temporal_positions: vector: center location of the analysis window.
%   frame_shift: scalar: interval between frame locations
%   inst_freq_map: matrix: instantaneous frequency of each channel (Hz)
%   residual_map: relative resitual of each channel
%           (relative power, calibrated for Nuttall-11 window.)
%   observation_probability_map: matrix: probability time series of each
%           channel
%   amplitude_map: matrix: amplitude of filter output (absolute value)
%   channels_in_octave: scalar: number of wavelet channels in one octave
%   center_frequency_list: vector: list of center frequencies of wavelet
%           filters (Hz)

% Copyright 2016 Google Inc. All Rights Reserved
% Author: hidekik@google.com (Hideki Kawahara)

if nargin == 1, frame_shift = 0.005; end;

startTic = tic;
calibration_constant = 64; % This is for NuttallWindows(~, 11)
time_axis = f0_structure.temporal_positions;
fs = f0_structure.sampling_frequency;
temporal_positions = (0:frame_shift:time_axis(end))';
frame_index = ...
  min(f0_structure.sample_length,max(1, round(fs * temporal_positions)));
number_of_channels = length(f0_structure.center_frequency_list);
inst_freq_map = f0_structure.inst_freq_map(:, frame_index);
residual_map = f0_structure.residual_map(:, frame_index);
amplitude_map = f0_structure.amplitude_map(:, frame_index);
tmp_variance_map = ...
  ConvertPowerToDecibel(residual_map) + calibration_constant;
snr_map = 10.0 .^ (tmp_variance_map / 20);
best_weight = snr_map * 0;
for ii = 1 : length(temporal_positions)
  best_weight(:, ii) = GetBestMixingWeight(snr_map(:, ii));
end;
observation_probability_map = tmp_variance_map * 0;
center_frequency_list = f0_structure.center_frequency_list;
tmp_half_bw = sqrt(center_frequency_list(2) / center_frequency_list(1));
fixed_best_weight = best_weight .* snr_map;
fc_in_cent_U = 1200 * log2(center_frequency_list * tmp_half_bw);
fc_in_cent_L = 1200 * log2(center_frequency_list / tmp_half_bw);
for ii = 1:length(temporal_positions)
  for kk = 1 : number_of_channels
    fqi_cent = 1200 * log2(max(20, inst_freq_map(kk, ii)));
    tmpU = erf((fc_in_cent_U - fqi_cent) / snr_map(kk, ii) / sqrt(2));
    tmpL = erf((fc_in_cent_L - fqi_cent) / snr_map(kk, ii) / sqrt(2));
    tmpProb = (tmpU - tmpL) / sum(tmpU - tmpL);
    observation_probability_map(:, ii) = ...
      observation_probability_map(:, ii) ...
      + fixed_best_weight(kk, ii) * tmpProb(:);
  end;
  observation_probability_map(:, ii) = observation_probability_map(:, ii) ...
    / sum(observation_probability_map(:, ii));
end;
output = struct('temporal_positions', temporal_positions, ...
                'frame_shift', frame_shift', ...
                'inst_freq_map', inst_freq_map, ...
                'residual_map', residual_map, ...
                'observation_probability_map', observation_probability_map, ...
                'amplitude_map', amplitude_map, ...
                'channels_in_octave', f0_structure.channels_in_octave', ...
                'center_frequency_list', center_frequency_list);
end
