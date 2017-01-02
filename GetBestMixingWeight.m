function best_weights = GetBestMixingWeight(stand_dev_of_variables)
% Best mixing of random variables with common mean
% and built-in constraint
%   best_weights = GetBestMixingWeight(stand_dev_of_variables)
%
% Input argument
%   stand_dev_of_variables: standard deviation of each variable
%
% Output argument
%   best_weights: best set of weights for mixing random variables

% Copyright 2016 Google Inc. All Rights Reserved.
% Author: hidekik@google.com (Hideki Kawahara)

narginchk(1, 1);
n_variables = length(stand_dev_of_variables);
H = ones(n_variables - 1, n_variables - 1) * ...
  stand_dev_of_variables(n_variables) ^ 2;
H = H + diag(stand_dev_of_variables(1 : n_variables - 1) .^ 2);
v = ones(n_variables - 1, 1) * stand_dev_of_variables(n_variables) ^ 2;
a = H \ v;
best_weights = [a; 1 - sum(a)];
end
