function output = ...
  StabilizeHarmonicDeviations(tx, harmonic_deviation, event_locations)
% StabilizeHarmonicDeviations -- remove statistical fluctuation of
% random component by smoothing
%
% output = StabilizeHarmonicDeviations(tx, event_locations)
%
% Argument
%   tx : time axis in audio sampling rate (s)
%   harmonic_deviation : random component level in audio rate (dB)
%   event_locations : excitation event location of periodic component (s)
%
% Return value
%   output : temporally smoothed random component

% Copyright 2016 Google Inc. All Rights Reserved
% Author: hidekik@google.com (Hideki Kawahara)

kSpacing = 0.01; % 10 ms separation for smoothing
harmonics_deviation_at_event = ...
  interp1(tx, harmonic_deviation', event_locations, 'linear', 'extrap')';
harmonics_deviation_at_event_B = ...
  interp1(tx, harmonic_deviation', event_locations - kSpacing, 'linear', ...
          'extrap')';
harmonics_deviation_at_event_F = ...
  interp1(tx, harmonic_deviation', event_locations + kSpacing, 'linear', ...
          'extrap')';
output = 0.5 * harmonics_deviation_at_event ...
  + 0.25 * (harmonics_deviation_at_event_B + ...
  harmonics_deviation_at_event_F); % this seems not very effective
end