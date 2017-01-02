function output = DefaultTrackingOptions
% Default tracking condition generation
%
% output = DefaultTrackingOptions
%
% Return value
%  output : structur variable with following fields
%    peak_spread : Gaussion SD for transition
%    deviation_spread : Gaussian SD for weighting around median F0
%    focus_spread : Gaussian SD for vicinity peak selection

% Copyright 2016 Google Inc. All Rights Reserved
% Author: hidekik@google.com (Hideki Kawahara)

output = struct('peak_spread', 0.5, ...
                'deviation_spread', 0.75, ...
                'focus_spread', 0.2);
end
