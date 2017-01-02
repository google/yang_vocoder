function validation_ok = ...
  ValidateMaxSegmentalSnrDifference(filename1, filename2, acceptable_error)
% compare two wave files and validate if their difference is within
% the given limit
%
% validation_OK = ...
%    ValidateMaxRelativeDifference(filename1, filename2, acceptable_error)
%
% Argument
%   filename1, filename2: full path of tested files
%   acceptable_error: acceptable segmental SNR difference
%
% Return value
%   validation_OK: validation result, true: passed, false: failed

% Copyright 2016 Google Inc. All Rights Reserved
% Author: hidekik@google.com (Hideki Kawahara)

validation_ok = false;
frame_length = 0.010; % 10 ms default
[x1, fs1] = AudioRead(filename1);
[x2, fs2] = AudioRead(filename2);
if fs1 ~= fs2
  return;
end;
n_data = min(length(x1), length(x2));
frame_sample_length = round(frame_length * fs1);
frame_locations = 1:frame_sample_length:n_data - frame_sample_length - 1;
scan_index = 0:frame_sample_length - 1;
max_snr_error = 0;
for ii = 1:length(frame_locations)
  segment_index = frame_locations(ii) + scan_index;
  [relative_snr_error, segment_power1, segment_power2] = ...
    ComputeSNR(x1, x2, segment_index);
  if relative_snr_error > max_snr_error && ...
      segment_power1 + segment_power2 > 0
    max_snr_error = relative_snr_error;
  end;
end;
if max_snr_error < acceptable_error
  validation_ok = true;
end
end

function [relative_snr_error, segment_power1, segment_power2] = ...
  ComputeSNR(x1, x2, segment_index)
segment_power1 = sum(x1(segment_index) .^ 2);
segment_power2 = sum(x2(segment_index) .^ 2);
error_power = sum((x2(segment_index) - x1(segment_index)) .^ 2);
relative_snr_error = ...
  sqrt(2 * error_power / (segment_power1 + segment_power2));
end
