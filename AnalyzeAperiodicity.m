function output = ...
  AnalyzeAperiodicity(x, fs, f0, frame_time, analysis_condition)
% F0 refinement and aperiodicity analysis
%
% output = ...
%   AnalyzeAperiodicity(x, fs, f0, frame_time, analysis_condition)
%
%  Input argument
%    x: speech signal, one dimensional column vector
%    fs: sampling frequency (Hz)
%    f0: initial F0 estimate (Hz)
%    frame_time: time location of the center of the analysis window
%    analysis_condition: structure variable with following fields
%      n_harmonics_for_H : number of harmonic components to use in
%        harmonics based F0 refinement
%      window_magnification_for_H : time window stretching facotor from
%        the permissible shortest length
%      enable_downsampling_for_H : 1: enables down sampling before analysis
%      conditions_for_T : structure variable array for defining condition
%        for time-warping based F0 refinement with the following fields
%        n_harmonics : number of harmonics to be uses
%        window_magnification : time window stretching facotor from
%          the permissible shortest length
%        enable_downsampling : 1: enables down sampling before analysis
%
%  Return value
%    output : structure variable with the following fields
%      aperiodicity_matrix : aperiodicity of each harmonic component
%      refined_f0 : refined F0 (Hz)
%      frame_time : analysis window center location at each frame (s)

% Copyright 2016 Google Inc. All Rights Reserved
% Author: hidekik@google.com (Hideki Kawahara)

narginchk(4, 5);
if nargin < 5, analysis_condition = GenerateAperiodicityAnalysisCondition; end;
start_tic = tic;
n_harmonics_for_H = analysis_condition.n_harmonics_for_H;
window_magnification_for_H = analysis_condition.window_magnification_for_H;
enable_downsampling_for_H = analysis_condition.enable_downsampling_for_H;
% Since at the very beginning F0 information is not very reliable,
% it is necessary to refine it without using F0 adaptive time warping
% Effects of time warping propagates to all the following process.
% The fundamental component is usually damaged by microphone low-cut for
% preventing proximity effect and background noise is stronger in lower
% frequency region. In addition, time warping procedure is fragile.
% It should be used carefully.
[refined_f0, ~] = RefineF0byHarmonics(x, fs, f0, frame_time, ...
  n_harmonics_for_H, window_magnification_for_H, enable_downsampling_for_H);
% Recursively apply time-warping based F0 refinement
conditions_for_T = analysis_condition.conditions_for_T;
n_steps = size(conditions_for_T, 2);
f0_last = refined_f0;
for ii = 1:n_steps
  [refined_f0, aperiodicity_matrix] = ...
    RefineF0byTimeWarping(x, fs, f0_last, frame_time, ...
                          conditions_for_T(ii).n_harmonics , ...
                          conditions_for_T(ii).window_magnification, ...
                          conditions_for_T(ii).enable_downsampling);
  f0_last = refined_f0;
end;
output = struct('aperiodicity_matrix', aperiodicity_matrix, ...
                'refined_f0', f0_last, 'sampling_frequency', fs, ...
                'frame_time', frame_time, 'elapsed_time', toc(start_tic));
end
