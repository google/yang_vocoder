function boundary_level = GetDefaultVuvThreshold(rms_dB, amount)
% set obvious VUV threshold and amount boundary based on Otsu method.
%
% boundary_level = GetDefaultVuvThreshold(rms_dB, amount)
%
% Input argument
%   rms_dB : signal rms level (dB)
%   amount : amount of selected voiced segment : suggestion 0.6
%
% Return value
%   boundary_level : threshold level to meet the condition (dB)
%
% Reference
% Otsu, N. A thresholding selection method from gray-level histogram.
% IEEE Transactions on Systems, Man and Cybernetics 9 (1979): 62-66.
%

n_class = 30;
[count, bin] = hist(rms_dB, n_class);
p = count / sum(count);
%%

omega_0 = zeros(n_class, 1);
omega_1 = zeros(n_class, 1);
mu_0 = zeros(n_class, 1);
mu_1 = zeros(n_class, 1);
sigma_B2 = zeros(n_class, 1);
mu_T = sum(bin .* p);
for k = 1:n_class - 1
   omega_0(k) = sum(p(1:k));
   mu_0(k) = sum(p(1:k) .* bin(1:k)) / omega_0(k);
   omega_1(k) = sum(p(k+1:n_class));
   mu_1(k) = sum(p(k+1:n_class) .* bin(k+1:n_class)) / omega_1(k);
   sigma_B2(k) = (mu_T * omega_0(k) - sum(p(1:k) .* bin(1:k))) ^2 / ...
       (omega_0(k) * (1 - omega_0(k)));
end;
[~, max_bin_index] = max(sigma_B2);
sorted_high = sort(rms_dB(rms_dB > bin(max_bin_index)));
n_high = length(sorted_high);
boundary_level = sorted_high(round(n_high * (1 - amount)));
end
