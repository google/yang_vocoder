function output = GetSigmoidNoiseShaper(frequency_axis, f_corner, f_slope)
%   output = GetSigmoidNoiseShaper(frequency_axis, f_corner, f_slope)
%   Power spectrum energy allocation to pulse and random component
%   using sinusoid. Sum of components has constant gain (0dB).
%   input
%       frequency_axis  : frequency axis vector (Hz)
%       f_corner        : transition frequendcy (Hz) -3dB point
%       f_slope         : slope of transition (cent)

% Copyright 2016 Google Inc. All Rights Reserved
% Author: hidekik@google.com (Hideki Kawahara)

frequency_axis(1) = frequency_axis(2);
cent_axis = 1200 * log2(abs(frequency_axis));
output = ...
  1.0 ./ (1 + exp(-(cent_axis - 1200 * log2(f_corner)) / f_slope)) + ...
  0.0000000001;
end
