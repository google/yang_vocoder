function output = GenerateVUVSetting
% Default condition generation for V/UV decision
%
% output = GenerateVUVSetting
%

% Copyright 2016 Google Inc. All Rights Reserved
% Author: hidekik@google.com (Hideki Kawahara)

output = struct('n_harmonics', 2, ...
                'smooth_length', 0.02, ... % 20 ms Hanning
                'slient_threshold_vuv', -30, ... % dB
                'lower_level_threshold', -25, ... % dB
                'low_level_bias', -5, ... % dB
                'unvoiced_threshold_vuv', -20); % dB
end
