function output = ...
  CalculateAperiodicitySgram(source)
% Calculate aperiodicity spectrogram from harmonic information
%
%  output = CalculateAperiodicitySgram(source)
%
%  Input argumrnt
%    source : structure, result with time warp, following fields are used
%      sampling_frequency : sampling frequency (Hz)
%      refined_f0 : F0 trajectory, vector (Hz)
%      frame_time : analysis frame center position (s)
%      aperiodicity_matrix : aperiodicity of each harmonic component (dB)
%
%  Return value
%    output : frequency smoothed aperiodicity spectrogram
%         (default smoothing width is one octave. Hard coded.)
%

% Copyright 2016 Google Inc. All Rights Reserved
% Author: hidekik@google.com (Hideki Kawahara)

narginchk(1, 1);
mag_factor = 1.2;
fs = source.sampling_frequency;
f0 = source.refined_f0;
f0_floor = min(f0);
tx = source.frame_time;
n_frames = length(tx);
smoothing_bw = 2 ^ (1/2); % width is 1 octave
%
fftl = 2 ^ ceil(log2(fs / f0_floor * 4 * 2 * mag_factor + 2));
fx = (0:fftl-1)'/fftl * fs;
fx(fx > fs/2) = fx(fx > fs /2 )-fs;
frequency_axis = fx(1:fftl / 2 + 1);
frequency_U = frequency_axis * smoothing_bw;
frequency_L = frequency_axis / smoothing_bw;
raw_aperiodicity = zeros(fftl / 2 + 1, n_frames);
smoothed_aperiodicity = zeros(fftl / 2 + 1, n_frames);
bin_picker = 1:size(source.aperiodicity_matrix,1);
for ii = 1:n_frames
  interpolated_aperiodicity = ...
    interp1(bin_picker * f0(ii), ...
            source.aperiodicity_matrix(:, ii), ...
            frequency_axis, 'linear', 'extrap');
  raw_aperiodicity(:, ii) = interpolated_aperiodicity(:);
  smoothed_aperiodicity(:, ii) = ...
    SmoothProportionalToFrequency(frequency_axis, frequency_U, frequency_L, ...
                                raw_aperiodicity(:, ii));
end;
%output = struct(...
%  'raw_aperiodicity', raw_aperiodicity,...
%  'smoothed_aperiodicity', smoothed_aperiodicity,...
%  'temporal_positions', tx,...
%  'elapsed_time', toc(start_tic));
output = smoothed_aperiodicity;
end

function output = ...
  SmoothProportionalToFrequency(frequency_axis, frequency_U, frequency_L, ...
                                raw_aperiodicity_slice)
% Triangular smoothing as a result of recursive rectangular smoothing
% This is for better sounding random component
power_slice = ConvertDecibelToPower(raw_aperiodicity_slice);
tmp = ...
  SmoothFrequencyProportionalRect(frequency_axis, frequency_U, frequency_L, ...
                                  power_slice);
tmp = ...
  SmoothFrequencyProportionalRect(frequency_axis, frequency_U, frequency_L, ...
                                  tmp);
output = ConvertPowerToDecibel(tmp);
output(1) = raw_aperiodicity_slice(1);
end

function output = ...
  SmoothFrequencyProportionalRect(frequency_axis, frequency_U, frequency_L, ...
                                power_slice)
% Rectangular smoothing
tmp_cumulative_power = cumsum(power_slice * frequency_axis(2));
power_U = interp1(frequency_axis, tmp_cumulative_power, frequency_U, ...
                  'linear','extrap');
power_L = interp1(frequency_axis, tmp_cumulative_power, frequency_L, ...
                  'linear','extrap');
output = tmp_cumulative_power;
output(2:end) = (power_U(2:end) - power_L(2:end)) ./ ...
  (frequency_U(2:end) - frequency_L(2:end));
end
