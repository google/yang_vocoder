function output = ...
  AnalyzeSpeechSource(x_original, fs_original, conditions)
% Execute fundamental frequency analysis and voiced / unvoiced detection
%
%  output = ...
%    AnalyzeSpeechSource(x_original, fs_original)
%
%  output = ...
%    AnalyzeSpeechSource(x_original, fs_original, conditions)
%
%  Input argumrnt
%    x_original : speech data, vertical vector, one dimensional
%    fs_original : sampling frequency (Hz)
%    conditions : structure variable following fields
%      f0_ceiling : higher F0 search range (Hz)
%      enable_downsampling : default 1, 1: downsampling OK
%      frame_shift : frame shift interval (s)
%
%  Return value
%    output : structure variable
%      f0 : fundamental frequency (Hz), only from fundamental
%           when refined, refined F0 updates this value
%      vuv : voicing detection result 1: voiced
%      f0_initial: initial estimate of F0 (Hz)
%      used_sampling_frequenency: sampling frequency used for initial
%                      estimate (Hz)
%      sampling_frequency: sampling frequency of the signal (Hz)
%
%      frame_time : frame center pisition (s)
%      probability_structure : probability map and other data
%        (see ConvertToObservationProbatility.m )
%      best_channel_out : selected best channel and other data
%        (see trackPractically.m)
%
%      instantaneous_frequency_structure: analysed data in audio sample
%                   rate, for debug

% Copyright 2016 Google Inc. All Rights Reserved
% Author: hidekik@google.com (Hideki Kawahara)

%% Preparation,
% For seed up processing, downsample input signal if it is allowed.
% Read analysis conditions if it is assigned.
narginchk(2, 3);
if nargin == 2
  conditions = GenerateOptionForSourceAnalysisNV;
end;
highest_frequency = conditions.f0_ceiling * 4;
original_Nyquist = fs_original / 2;
n_downsampling = ceil(original_Nyquist / highest_frequency);
if n_downsampling > 1 && conditions.enable_downsampling == 1
  conditions.sampling_frequency = ...
    fs_original / n_downsampling;
  x = DecimateUsingFFT(x_original, n_downsampling);
  fs = fs_original / n_downsampling;
else
  conditions.sampling_frequency = fs_original;
  x = x_original;
  fs = fs_original;
end;
frame_shift = conditions.frame_shift;
%% Front-end analysis using wavelet and produc probability map
wavelet = DesignWaveletFilterBank(conditions);
instantaneous_frequency_structure = ...
  AnalyzeInstantaneousFrequency(x, fs, wavelet);
probability_structure = ...
  ConvertToObservationProbatility(instantaneous_frequency_structure, ...
                                  frame_shift);
%% select best channel trajectory consisting of F0 component
median_f0 = GetMedianF0Parameters(probability_structure);
best_channel_structure = ...
  FindBestTrajectory(probability_structure, median_f0, ...
                     conditions.tracking_condition);
%% estimate initial F0 value using the best channel trajectory
initial_f0_structure = ...
  EstimateInitialF0AndDeviation(probability_structure, ...
                                best_channel_structure.trajectory);
%% VUV analysis
median_f0 = probability_structure.center_frequency_list(...
  best_channel_structure.median_index);
frame_is_unvoiced = ...
  DecideVUV(x_original, fs_original, probability_structure, ...
            median_f0, conditions.vuv_setting);
%% F0 refinement and aperiodicity analysis
if conditions.enable_f0_refinement
  aperiodicity = ...
    AnalyzeAperiodicity(x_original, fs_original, ...
                        initial_f0_structure.f0_initial, ...
                        initial_f0_structure.temporal_positions, ...
                        conditions.refinement_analysis_condition);
end;
%% Complete return value fields
output = struct('sampling_frequency', fs_original, ...
                'frame_time', probability_structure.temporal_positions, ...
                'probability_structure', probability_structure, ...
                'best_channel_out', best_channel_structure, ...
                'vuv', double(frame_is_unvoiced < 0.5));
if conditions.enable_f0_refinement
  output.f0 = aperiodicity.refined_f0;
  output.aperiodicity_matrix = aperiodicity.aperiodicity_matrix;
end;
end
