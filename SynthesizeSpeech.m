function output = ...
  SynthesizeSpeech(spectrum_envelope, source)
% Synthesize speech from harmonic information
%
% output = ...
%  SynthesizeSpeech(spectrum_envelope, source)
%
% Arguments
% spectrum_envelope: structure with following fields
%   harmonic_power_dB: Harmonics level in dB
%   sampling_frequency: sampling frequency 9Hz)
% source: structure with following fields
%   f0: fundamental frequency (Hz)
%   aperiodicity_matrix: residual level in dB (Has to be calibrated)
%   vuv : voiced/unvoiced indicator of each frame: 1: voiced
%
% Return value
%   output : synthesized signal

% Copyright 2016 Google Inc. All Rights Reserved
% Author: hidekik@google.com (Hideki Kawahara)

% check input and prepare information for synthesis
narginchk(2, 2);
tx = source.frame_time;
f0 = source.f0;
fs = source.sampling_frequency;
harmonic_level = spectrum_envelope.harmonic_power_dB;
harmonic_deviation = source.aperiodicity_matrix;
unvoiced_frame = double(source.vuv < 0.5);
[event_index, event_locations, f0i] = CalculateEventLocations(tx, fs, f0);
unvoiced_sample = ...
  interp1(tx, double(unvoiced_frame), event_locations, 'linear', 'extrap');
harmonics_level_at_event = ...
  interp1(tx, harmonic_level', event_locations, 'linear', 'extrap')';
harmonics_deviation_at_event = ...
  StabilizeHarmonicDeviations(tx, harmonic_deviation, event_locations);
output = ...
  SynthesizeByEventBasedMethod(event_index, event_locations,...
  harmonics_level_at_event, ...
  harmonics_deviation_at_event, fs, f0i, ...
  unvoiced_sample);
end

