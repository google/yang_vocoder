function output = CalculateSpectrumEnvelope(spectrum)
% Calculate aperiodicity spectrogram from harmonic information
%
%  output = CalculateSpectrumEnvelope(spectrum)
%
%  Argumrnts
%    spectrum : structure variable, following fields are used
%      sampling_frequency : sampling frequency (Hz)
%      used_f0: F0 used in spectrum analysis, vector (Hz)
%      temporal_positions : analysis frame center position (s)
%      harmonic_power_dB : power level of each harmonic component (dB)
%
%  Return value
%    output : interpolated spectrogram (power)

% Copyright 2016 Google Inc. All Rights Reserved
% Author: hidekik@google.com (Hideki Kawahara)

narginchk(1, 1);
mag_factor = 1.2;
f0 = spectrum.used_f0;
fs = spectrum.sampling_frequency;
f0_floor = min(f0);
tx = spectrum.temporal_positions;
n_frames = length(tx);
fftl = 2 ^ ceil(log2(fs / f0_floor * 4 * 2 * mag_factor + 2));
fx = (0:fftl - 1)' / fftl * fs;
fx(fx > fs / 2) = fx(fx > fs / 2 ) - fs;
frequency_axis = fx(1:fftl / 2 + 1);
sgram_spline = zeros(fftl / 2 + 1, n_frames);
for ii = 1:n_frames
  harmonic_locations = (f0(ii):f0(ii):fs / 2 - f0(ii) / 2)';
  n_harmonics = length(harmonic_locations);
  % TODO(hidekik): These two lines and interpolation needs refactoring together with
  % regeneration of reference files.
  % ****** fragile. refactor with care ******
  harmonic_values_dB = spectrum.harmonic_power_dB(:, ii);
  harmonic_values_dB = harmonic_values_dB(1:n_harmonics);
  % TODO(hidekik) This interpolation trick usig [2 1 1 ...] is not a good idea
  tmp = interp1([-2 * f0(ii); -f0(ii); 0; harmonic_locations; fs / 2], ...
                harmonic_values_dB([2 1 1 1:n_harmonics n_harmonics]), ...
                frequency_axis, 'spline', 'extrap'); % should be rev.
  sgram_spline(:, ii) = tmp(:);
end;
output = ConvertDecibelToPower(sgram_spline);
end
