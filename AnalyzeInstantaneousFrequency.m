function output = AnalyzeInstantaneousFrequency(x, fs, wavelet)
% Analize instantaneous frequency using wavelet filter bank
%   output = AnalyzeInstantaneousFrequency(x, fs, wavelet)
%
% Arguments
%
% x: column vector: input speech signal
% fs: scalar: sampling frequency of input (Hz)
% wavelet: structure with the following fields.
%   center_frequency_list: vector: center frequencies of filters (Hz)
%   sampling_frequency: scalar: sampling frequency (Hz)
%   channels_in_octave: scalar: number of wavelet channels in one octave
%   detector_bank: structure variable: designed wavelet filters with the
%       following fields
%     filter: main filter
%     filter_diff: paired filter
%     bias: offset samples from the center of the filter
%
% Return value
%
% output: structure with the following fields
%   temporal_positions: vector: center location of the analysis window.
%     represented in terms of second
%   sampling_frequency: scalar: sampling frequency (Hz)
%   sample_length: length of the input data (sample)
%   number_of_channels: number of filter channels
%   center_frequency_list: vector: center frequency of each wavelet filter
%        (Hz)
%   inst_freq_map: matrix: instantaneous frequency of each filter output
%        (Hz)
%   residual_map: matrix: relative level of random component (power ratio)
%   amplitude_map: matrix: amplitude of filter output (absolute value)

% Reference:
% H. Kawahara, Y. Agiomyrgiannakis, H. Zen: Using instantaneous frequency
% and aperiodicity detection to estimate F0 for high-quality speech
% synthesis, ISCA SSW9, 2016 [submitted]

% Copyright 2016 Google Inc. All Rights Reserved
% Author: hidekik@google.com (Hideki Kawahara)

narginchk(3, 3);
if fs ~= wavelet.sampling_frequency
  disp('Sampling frequency mismatch!');
  disp(['Sampling frequency for wavelet:' ...
    num2str(wavelet.sampling_frequency) ' Hz']);
  disp(['Sampling frequency for signal:' num2str(fs) ' Hz']);
  output = [];
  return;
end;
sample_length = length(x);
inst_freq_map = zeros(wavelet.number_of_channels, sample_length);
amplitude_map = zeros(wavelet.number_of_channels, sample_length);
residual_map = zeros(wavelet.number_of_channels, sample_length);
center_frequency_list = wavelet.center_frequency_list;
filler = zeros(wavelet.detector_bank(1).bias * 2, 1);
x_extended = [x; filler];
x_extended = x_extended ... % next is a safequard (not comp. with octave)
  + randn(length(x_extended),1) * std(x_extended) / 100000;
n_data = length(x);
% The following loop calcuates front-end aperiodicity detector outputs.
% Refer to SSW9 paper: Eqs. (16) to (23)
for ii = 1:wavelet.number_of_channels
  complex_filter = wavelet.detector_bank(ii).filter;
  complex_filter_diff = wavelet.detector_bank(ii).filter_diff;
  smoother = abs(complex_filter) / sum(abs(complex_filter));
  bias = wavelet.detector_bank(ii).bias;
  x1 = FilterWithoutTimeBias(complex_filter, x_extended,  bias, n_data);
  xd = FilterWithoutTimeBias(complex_filter_diff, x_extended, bias, n_data);
  x1_normal = x1 ./ abs(x1);
  x2 = FilterWithoutTimeShift(complex_filter, x1_normal, filler, bias);
  x2_normal = x2 ./ abs(x2);
  inst_freq_map(ii, :) = ...
    CalculateInstantaneousFrequency(x1, xd, center_frequency_list(ii))';
  amplitude_map(ii, :) = abs(x1)';
  residual_map(ii, :) = ... % Eq. (23)
    FilterWithoutTimeShift(smoother, abs(x1_normal - x2_normal) .^ 2, ...
                           filler, bias)';
end;
inst_freq_map = max(center_frequency_list(1), ...  % safeguard
                    min(center_frequency_list(end), inst_freq_map));
residual_map = max(0.00000001, residual_map); % safeguard
% ----
output = struct('center_frequency_list', center_frequency_list, ...
                'sampling_frequency', fs, ...
                'number_of_channels', length(center_frequency_list), ...
                'channels_in_octave', wavelet.channels_in_octave, ...
                'sample_length', sample_length, ...
                'temporal_positions', (0:sample_length-1)' / fs, ...
                'inst_freq_map', inst_freq_map, ...
                'amplitude_map', amplitude_map, ...
                'residual_map', residual_map);
end

function output = CalculateInstantaneousFrequency(x1r,xd, center_frequency)
% Intantaneous frequency calculation using Eq.(16) of SSW9 paper.
  numrr = abs(x1r) .^ 2;
  denomr = real(x1r) .* imag(xd) - imag(x1r) .* real(xd);
  output = (denomr ./ numrr) / (2 * pi) + center_frequency;
end

function y = FilterWithoutTimeShift(fir_filter, x, filler, bias)
y = fftfilt(fir_filter, [x; filler]);
y = y(bias + (1:length(x)));
end

function y = FilterWithoutTimeBias(fir_filter, x, bias, n)
y = fftfilt(fir_filter, x);
y = y(bias + (1:n));
end
