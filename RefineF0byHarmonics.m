function [refined_f0, spectra] = ...
  RefineF0byHarmonics(x_original, fs_original, f0_initial, frame_time, ...
                      n_harmonics, mag_factor, downsample_enable)
% refine F0 using harmonic components
%
% output = ...
%   RefineF0byHarmonics(x_original, fs_original, f0_initial, frame_time, ...
%                       n_harmonics, mag_factor, downsample_enable)
%
% Input argument
%   x_original : speech data
%   f0_initial : sampling frequency (Hz)
%   f0_initial : intial estimate of F0 sequence (Hz)
%   frame_time : center location of analysis window (s)
%   n_harmonics : number of harmonic components to be used for refinement
%   mag_factor : time window stretching factor
%   downsample_enable : flag for downsampling, 1: down sampling if
%                       applicable
%
% Return value
%   refined_f0 : refined F0 estimate (Hz)
%   spectra : structure variable. output of function
%      "AnalyzeSpeechSpectra" called inside this routine

% Copyright 2016 Google Inc. All Rights Reserved
% Author: hidekik@google.com (Hideki Kawahara)

narginchk(7, 7);
opt.mag_factor = mag_factor;
top_margin = 1.1;
threshold_relevant = 0.1;
f0_ceiling = max(f0_initial);
f0_floor = min(f0_initial);
if f0_ceiling == f0_floor
  f0_ceiling = f0_ceiling * 1.4;
  f0_floor = f0_floor / 1.4;
end;
n_downsampling = ceil(fs_original / 2 / f0_ceiling / n_harmonics / 4);
if n_downsampling > 1 && downsample_enable
  x = DecimateUsingFFT(x_original, n_downsampling);
  fs = fs_original / n_downsampling;
else
  x = x_original;
  fs = fs_original;
end;
spectra = ...
  AnalyzeSpeechSpectra(x, fs, f0_initial, frame_time, opt);
inst_freq_map = spectra.inst_freq_map;
residual_map = spectra.residual_sgram;
n_frames = length(frame_time);
frequency_axis = spectra.frequency_axis;
refined_f0 = f0_initial * 0;
total_SD = f0_initial * 0;
index_base = (1:n_harmonics)';
for ii = 1:n_frames
  upper_limit = f0_initial(ii) * n_harmonics * top_margin;
  tmp_harmonics = ...
    SelectValueAtHarmonics(frequency_axis, upper_limit, inst_freq_map, ...
                           ii, n_harmonics, f0_initial);
  tmp_residual = ...
    SelectValueAtHarmonics(frequency_axis, upper_limit, residual_map, ...
                           ii, n_harmonics, f0_initial);
  tmp_residual = max(0.005, tmp_residual);
  relevant_channel = tmp_residual < threshold_relevant;
  best_weights = ...
    GetBestMixingWeightForHarmonics(tmp_residual(relevant_channel) ./ ...
                                    index_base(relevant_channel));
  refined_f0(ii) = sum((tmp_harmonics(relevant_channel) ./ ...
                       index_base(relevant_channel)) .* best_weights);
  total_SD(ii) = sqrt(sum(tmp_residual(relevant_channel) .^2 ./ ...
                      index_base(relevant_channel) .* best_weights .^ 2));
  if abs((1 - refined_f0(ii) / f0_initial(ii))) > 0.12
    refined_f0(ii) = f0_initial(ii);
    total_SD(ii) = 1;
  end;
end;
refined_f0 = max(f0_floor, min(f0_ceiling, refined_f0));
end

function output = ...
  SelectValueAtHarmonics(frequency_axis, upper_limit, map, ii, ...
                         n_harmonics, f0)
output = ...
  interp1(frequency_axis(frequency_axis < upper_limit), ...
          map(frequency_axis < upper_limit, ii), ...
          (1:n_harmonics) * f0(ii), 'linear', 'extrap');
output = output(:);
end
