function output = ...
  SaveAnalysisResult(base_path, base_name, wave_file_name, ...
                     source, spectrum_out)
% Save analysis results to C-compatible binary files
%
% output = ...
%   SaveAnalysisResult(base_path, base_name, source, ...
%                      spectrum_out)
%
% Arguments
%   base_path : path to the directory for saving results
%   base_name : base neme of result files
%   wave_file_name : analysed file name
%   source : output structure variable of source analysis
%   spectrum_out : output structure variable of spectrum_out analysis
%     (Note that source and spectrum_out
%      are enough to resynthesize)
%
% Return value
%   output : structure variable consisting of result file names

% Copyright 2016 Google Inc. All Rights Reserved
% Author: hidekik@google.com (Hideki Kawahara)

output = struct(...
  'f0_name', [base_name wave_file_name(1:end-4) '_f0.bin'], ...
  'vuv_name', [base_name wave_file_name(1:end-4) '_vuv.bin'], ...
  'sp_name', [base_name wave_file_name(1:end-4) '_sp.bin'], ...
  'ap_name', [base_name wave_file_name(1:end-4) '_ap.bin']);
SaveMatrix(source.f0, [base_path output.f0_name]);
SaveMatrix(source.vuv, [base_path output.vuv_name]);
SaveMatrix(source.aperiodicity_matrix, [base_path output.ap_name]);
SaveMatrix(spectrum_out.harmonic_power_dB, [base_path output.sp_name]);
end
