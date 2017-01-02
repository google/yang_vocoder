function validation_OK = ...
  ValidateMaxRelativeDifference(filename1, filename2, acceptable_error)
% compare two wave files and validate if their difference is within
% the given limit
%
% validation_OK = ...
%    ValidateMaxRelativeDifference(filename1, filename2, acceptable_error)
%
% Argument
%   filename1, filename2: full path of tested files
%   acceptable_error: acceptable relative error
%
% Return value
%   validation_OK: validation result, true: passed, false: failed

% Copyright 2016 Google Inc. All Rights Reserved
% Author: hidekik@google.com (Hideki Kawahara)

validation_OK = false;
[x1, fs1] = AudioRead(filename1);
[x2, fs2] = AudioRead(filename2);
if fs1 ~= fs2
  return;
end;
n_data = min(length(x1), length(x2));
max_difference = max(abs(x1(1:n_data) - x2(1:n_data)));
relative_error = max_difference / max(max(abs(x1)), max(abs(x2)));
if relative_error < acceptable_error
  validation_OK = true;
end
end
