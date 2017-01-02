function output = GenerateOptionForSourceAnalysisNV
% Default parameter definition for wavelet based IF map calculation
%  output = GenerateOptionForSourceAnalysis
%
% Return value
%   output: structure variable: field names and meanings are below
%     f0_floor: scalar: lowest F0 search range (Hz)
%     f0_ceiling: scalar: highest F0 search range (Hz)
%     nInOctave: scalar: number of filters in one octave
%     enable_downsampling: downsampling of signals switch. 1: enabled
%     frame_shift: analysis frame shift (s). default 0.005 s
%     enable_f0_refinement: F0 refinement and aperiodicity analysis switch
%                           1: enabled (default)
%     refinement_analysis_condition: setting for refinment and aperiodicity
%                           See GenerateAperiodicityAnalysisCondition.m

% Copyright 2016 Google Inc. All Rights Reserved.
% Author: hidekik@google.com (Hideki Kawahara)

output = struct('f0_floor', 40, ...
                'f0_ceiling', 1000, ...
                'channels_in_octave', 12, ...
                'enable_downsampling', 1, ...
                'frame_shift', 0.005, ...
                'enable_f0_refinement', 1, ...
                'refinement_analysis_condition', ...
                GenerateAperiodicityAnalysisCondition, ...
                'tracking_condition', DefaultTrackingOptions, ...
                'vuv_setting', GenerateVUVSetting);
end
