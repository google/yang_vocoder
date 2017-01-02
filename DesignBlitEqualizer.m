function output = DesignBlitEqualizer(raw_BLIT_response)
% DesignBlitEqualizer -- Equalize raw BLIT (Band Limitted Impulse Train)
% to have flat gain using FIR equalizer.
% This equalizer assumes Nuttall-based BLIT with zero at 0.55 times
% sampling frequency.
% This is a supporting function for  fractional pitch allocation of
% excitation pulse
%
% output = DesignBlitEqualizer(raw_BLIT_response)
%
% Input argument
%   raw_BLIT_response: impulse response of antialiasing BLIT filter
%
% Output
%   output : FIR equalizer to flatten frequency response

% Copyright 2016 Google Inc. All Rights Reserved
% Author: hidekik@google.com (Hideki Kawahara)

fftl = 8192; % Finally the central 121 points are used.
% The truncation length 121 is decided for the response to decrease
% relatively small enough.
amp = abs(fft(raw_BLIT_response / sum(raw_BLIT_response), fftl));
eqAmp = 1.0 ./ (amp + 10.0 ^ (-45 / 20));
output = real(fftshift(ifft(eqAmp)));
output = output(fftl / 2 + 1 + (-60:60));
end
