function wavelet = DesignWaveletFilterBank(option_structure)
% Design wavelet filter bank for periodicity detection
%   wavelet = designWaveletFilterBank(option_structure)

% Arguments
%
% option_structure: structure variable with the following fields.
%   sampling_frequency: scalar: sampling frequency (Hz)
%   f0_floor: scalar: lowest F0 search range (Hz)
%   f0_ceiling: scalar: highest F0 search range (Hz)
%   nInOctave: scalar: number of filters in one octave
%
% Return value
%
% wavelet: structure variable with the following fields.
%   center_frequency_list: vector: center frequencies of filters (Hz)
%   sampling_frequency: scalar: sampling frequency (Hz)
%   channels_in_octave: scalar: number of wavelet channels in one octave
%   detector_bank: structure variable: designed wavelet filters with the
%       following fields
%     filter: main filter
%     filter_diff: paired filter
%     bias: offset samples from the center of the filter
%     number_of_channels: total number of filter channels

% Copyright 2016 Google Inc. All Rights Reserved
% Author: hidekik@google.com (Hideki Kawahara)

fs = option_structure.sampling_frequency;
center_frequency_list = option_structure.f0_floor * 2.0 .^ ...
  (0:1 / option_structure.channels_in_octave: ...
  log2(option_structure.f0_ceiling / option_structure.f0_floor));
%wavelet.center_frequency_list = center_frequency_list;
number_of_channels = length(center_frequency_list);
kMagnificationFactor = 2;
detector_bank = struct;
for ii = 1:number_of_channels
  fc = center_frequency_list(ii);
  half_length = round(kMagnificationFactor * fs / fc);
  tt2 = (-half_length:half_length)' / fs;
  bias = half_length;
  w = NuttallWindows(2 * half_length + 1, 11);
  w = w / sum(w);
  complex_wavelet = w .* exp(2 * pi * 1i * tt2 * fc);
  wd = diff(w) * fs;
  wdi = interp1(tt2(1:end - 1) + 1 / fs / 2, wd, tt2, 'linear', 'extrap')';
  complex_wavelet_diff = exp(1i * 2 * pi * fc * tt2) .* wdi(:);

  detector_bank(ii).filter = complex_wavelet;
  detector_bank(ii).filter_diff = complex_wavelet_diff;
  detector_bank(ii).bias = bias;
end;
wavelet = struct('detector_bank', detector_bank, ...
                 'number_of_channels', number_of_channels, ...
                 'sampling_frequency', fs, ...
                 'center_frequency_list', center_frequency_list, ...
                 'channels_in_octave', option_structure.channels_in_octave);
return;
