function y = DecimateUsingFFT(x, r)
%   This version is minimaly reduced implementation of decimate
%
%   y = DecimateUsingFFT(x, r)
%
% Input argument
%   x : data to be resampled, vector
%   r : decimation ratio, integer

% TODO(hidekik) To be refactored to fftl = 128, because total response length
% is less than 6 ms and smearing of temporal resolution is small.
% Gain at resampled Nyquist frequency is -100dB, meaning aliasing is negligible.
% But this may also affects regression tests and modified carefully

fftl = 512;
fx = (0:fftl - 1) / fftl;
fx(fx > 1/2) = fx(fx > 1 / 2) - 1;
fc = 0.45 / r;
g = double(abs(fx) < fc)';
fir = fftshift(real(ifft(g))) .* NuttallWindows(fftl, 12);
y = fftfilt(fir,[x(:); zeros(fftl, 1)]);
y = y(fftl / 2 + (1:length(x)));
y = y(1:r:end);
end
