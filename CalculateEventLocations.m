function [event_index, event_locations, f0i] = ...
  CalculateEventLocations(tx, fs, f0)
% CalculateEventLocations -- calculate excitation event locations from F0
%
% [event_index, event_locations] = CalculateEventLocations(tx, fs, f0)
%
% Argument
%   tx : time axis in audio sample rate (s)
%   fs : sampling frequency (Hz)
%   f0 : fundamental frequency in frame rate (Hz)
%
% Return value
%   event_index : excitation event location in sampled bin
%   event_locations : excitation event location (s)
%   f0i : fundamental frequency in audio sample rate (Hz)

% Copyright 2016 Google Inc. All Rights Reserved
% Author: hidekik@google.com (Hideki Kawahara)

tmp_sample_time = (0:1 / fs:tx(end))';
f0_floor = min(f0);
f0_ceiling = max(f0);
tmp_n_data = length(tmp_sample_time);
f0i = max(f0_floor, ...
          min(f0_ceiling, ...
              interp1(tx, f0, tmp_sample_time, 'linear', 'extrap')));
phase = cumsum(2 * pi * f0i / fs);
tmp_cos = cos(phase);
tmp_index = (1:tmp_n_data)';
event_index = tmp_index(tmp_cos([1 1:end - 1]) < tmp_cos & ...
  tmp_cos([2:end end]) <= tmp_cos);
event_index = event_index(event_index > 1 & event_index < tmp_n_data);
event_locations = ...
  (-0.5 + abs(tmp_cos(event_index) - tmp_cos(event_index - 1)) ./ ...
  (abs(tmp_cos(event_index) - tmp_cos(event_index - 1)) + ...
  abs(tmp_cos(event_index) - tmp_cos(event_index + 1))) - ...
  1 + event_index) / fs; % palabolic interpolation
end
