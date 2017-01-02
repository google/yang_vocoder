function output = AnalyzeSpeechSpectra(x, fs, f0, frame_time, opt)
%  AnalyzeSpeechSpectra -- Calculate frame-based parameters
%
%  output = AnalyzeSpeechSpectra(x, fs, f0, frame_time)
%
%  output = AnalyzeSpeechSpectra(x, fs, f0, frame_time, opt)
%
% Argument
%
%   x: vector: speech data
%   fs: scalar: sampling frequency (Hz)
%   f0: vector: fundamental frequency (Hz)
%   frame_time: vector: time of each frame (s)
%   opt: structure: optional, fields are as follows
%     mag_factor: scalar: time window stretching factor, default: 1.2
%
% Return value
%
% output: structur with following fields
%   residual_sgram: for debug
%   inst_freq_map: for deug
%   temporal_positions: frame time in second
%   spectral_envelope_dB: spectrum envelope (dB)
%   harmonic_power_dB: power of each harmonic component (dB)
%   frequency_axis: frequency axis of sgram_spline_power (Hz)

% Copyright 2016 Google Inc. All Rights Reserved
% Author: hidekik@google.com (Hideki Kawahara)

narginchk(4, 5);
if nargin == 5
  if isfield(opt, 'mag_factor')
    mag_factor = opt.mag_factor;
  end;
else
  mag_factor = 1.2;
end;
nuttall_coeff = [0.338946 0.481973 0.161054 0.018027]';
f0_floor = min(f0);
fftl = 2 ^ ceil(log2(fs / f0_floor * 4 * 2 * mag_factor + 2));
n_frames = length(f0);
res_gram = zeros(fftl, n_frames);
s_gram = zeros(fftl, n_frames);
fx = (0:fftl-1)'/fftl * fs;
fx(fx > fs/2) = fx(fx > fs /2 )-fs;
freq_axis = fx(1:fftl / 2 + 1);
sgram_spline = zeros(fftl / 2 + 1, n_frames);
inst_freq_map = zeros(fftl / 2 + 1, n_frames);
power_matrix = zeros(ceil(fs / 2 / f0_floor),n_frames);
for ii = 1:n_frames
  index_time = max(1, min(length(x), round(frame_time(ii)*fs)));
  scan_index = (-round(2 * fs / f0(ii) * mag_factor):round(2*fs/f0(ii) * mag_factor))';
  t_scan = scan_index / fs * f0(ii) / 4 / mag_factor;
  scan_index2 = (-round(2 * fs / f0(ii) * mag_factor) * 2: ...
    2 * round(2 * fs / f0(ii) * mag_factor))';
  w = cos(2 * pi * t_scan * [0 1 2 3]) * nuttall_coeff;
  bias1 = round(2 * fs / f0(ii) * mag_factor);
  w2 = conv(w, w);
  bias2 = bias1 * 2;
  [spec1, unit_spec1, segment] = ...
    AnalyzeSpectraAtWindowCenter(x, fs, w, fx, index_time, scan_index, ...
    fftl, bias1);
  [~, unit_spec2, ~] = ...
    AnalyzeSpectraAtWindowCenter(x, fs, w2, fx, index_time, scan_index2, ...
    fftl, bias2);
  res_spec = abs(unit_spec1 - unit_spec2);
  res_gram(:,ii) = res_spec;
  s_gram(:,ii) = spec1;
  [harmonic_values_dB, harmonic_locations, n_harmonics] = ...
    AssignPowerToHarmonics(spec1, fftl, f0(ii), fs, freq_axis);
  power_matrix(:, ii) = harmonic_values_dB(end);
  power_matrix(1:length(harmonic_values_dB), ii) = harmonic_values_dB(:);
  tmp = interp1([-2 * f0(ii); -f0(ii); 0; harmonic_locations; fs / 2], ...
    harmonic_values_dB([2 1 1 1:n_harmonics n_harmonics]), ...
    freq_axis, 'spline', 'extrap');
  %--- calculate dB spectral envelope
  sgram_spline(:, ii) = tmp(:);
  %--- calculate instantaneous frequency
  wd = sin(2 * pi * t_scan * [1 2 3]) * (nuttall_coeff(2:4) .* [1 2 3]');
  wd = + wd * pi / (2 / f0(ii) * mag_factor);
  x1 = spec1;
  xd = fft(wd .* segment, fftl);
  tmp_inst_freq = ...
    ((real(x1 .* imag(xd) - imag(x1) .* real(xd))) ./ abs(x1) .^ 2) ...
    / 2 / pi + fx;
  inst_freq_map(:,ii) = tmp_inst_freq(1:fftl / 2 + 1);
end;
output = struct('residual_sgram', res_gram(1:fftl / 2 + 1, :), ...
  'inst_freq_map', inst_freq_map, ...
  'temporal_positions', frame_time, ...
  'spectral_envelope_dB', sgram_spline, ...
  'harmonic_power_dB', power_matrix, ...
  'frequency_axis', freq_axis);
end

function [spec1, unit_spec1, segment] = ...
  AnalyzeSpectraAtWindowCenter(x, fs, w, fx, index_time, scan_index, ...
  fftl, bias1)
segment = x(max(1, min(length(x), index_time + scan_index)));
phase_bias = exp(1i * 2 * pi * bias1 / fs * fx);
spec1 = fft(segment .* w, fftl);
unit_spec1 = spec1 ./ abs(spec1) .* phase_bias;
end

function [power_dB, harmonic_locations, n_harmonics] = ...
  AssignPowerToHarmonics(spec1, fftl, current_f0, fs, freq_axis)
power_spectrum = abs(spec1(1:fftl / 2 + 1)) .^ 2;
cumulated_spectrum = cumsum(power_spectrum);
harmonic_locations = (current_f0:current_f0:fs / 2 - current_f0 / 2)';
n_harmonics = length(harmonic_locations);
harmonic_values_U = ...
  interp1(freq_axis, cumulated_spectrum, harmonic_locations + ...
  current_f0 / 2, 'linear', 'extrap');
harmonic_values_L = ...
  interp1(freq_axis, cumulated_spectrum, harmonic_locations - ...
  current_f0 / 2, 'linear', 'extrap');
power_dB = ConvertPowerToDecibel((harmonic_values_U - harmonic_values_L) / ...
  (current_f0 / freq_axis(2)));
end
