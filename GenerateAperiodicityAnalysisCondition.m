function output = GenerateAperiodicityAnalysisCondition
% Default analysis condition generation for aperiodicity analysis
%
% output = GenerateAperiodicityAnalysisCondition
%
% This defines operator H_m and T_m in SSW9 paper

% Copyright 2016 Google Inc. All Rights Reserved
% Author: hidekik@google.com (Hideki Kawahara)


conditions_for_T = struct('n_harmonics', {4, 4, 4}, ...
                         'window_magnification', {1.5, 1.5, 1.5}, ...
                         'enable_downsampling', {1, 1, 0});
output = struct('n_harmonics_for_H', 3, ...
                'window_magnification_for_H', 1.5, ...
                'enable_downsampling_for_H', 1, ...
                'conditions_for_T', conditions_for_T, ...
                'generate_aperiodicity_sgram', 1);
end
