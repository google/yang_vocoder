function best_weights = GetBestMixingWeightForHarmonics(standard_deviation_list)
% Best mixing of instantaneous frequencies with common mean
% and built-in constraint (safeguard version designed for harmonics information mixing)
%
% best_weights = GetBestMixingWeightForHarmonics(sdListOrg)
%
% Input argument
%   standard_deviation_list    : standard deviation of each variable
%
% Output argument
%   best_weights   : best set of weights for mixing random variables
%
% Note
%   This mixing consists of heuristics. Re-design is needed.

% Copyright 2016 Google Inc. All Rights Reserved.
% Author: hidekik@google.com (Hideki Kawahara)

original_list_length = length(standard_deviation_list);
safe_list = (1:original_list_length)';
safe_list = safe_list(standard_deviation_list > 0 & ...
  standard_deviation_list < 100 * min(abs(standard_deviation_list)));
sdLsafe_SD_list = standard_deviation_list(safe_list);
n_of_safe_SD = length(sdLsafe_SD_list);
if n_of_safe_SD > 0
  H = ones(n_of_safe_SD-1, n_of_safe_SD-1) * sdLsafe_SD_list(n_of_safe_SD) ^ 2;
  H = H + diag(sdLsafe_SD_list(1:n_of_safe_SD - 1) .^ 2);
  v = ones(n_of_safe_SD - 1, 1) * sdLsafe_SD_list(n_of_safe_SD) ^ 2;
  a = H \ v; % revised 23/May/2016 from inv(H) * v
  if sum(a) > 1
    a = a / sum(a);
  end;
  w = [a; 1 - sum(a)];
  best_weights = zeros(original_list_length, 1);
  best_weights(safe_list) = w;
else
  best_weights = zeros(original_list_length, 1);
  best_weights(1) = 1;
end;
end
