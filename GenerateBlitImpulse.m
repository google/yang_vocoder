function blit_impulse = ...
  GenerateBlitImpulse(fractional_index, half_width_of_BLIT)
%  Generate a band limited pulse (BLIT) using Nuttall window
%    blit_impulse = ...
%      GenerateBlitImpulse(fractional_index, half_width_of_BLIT)
%
% Input argument
%   fractional_index : non-integer sample index around pulse position
%   half_width_of_BLIT : normalized width of the anti-aliased impulse

%   Reference:
%   [1] Nuttall, A.; ``Some windows with very good sidelobe behavior,''
%   Acoustics, Speech and Signal Processing, IEEE Transactions on ,
%   vol.29, no.1, pp. 84--91, Feb 1981.
%   doi: 10.1109/TASSP.1981.1163506

narginchk(2, 2);
aa = [0.355768 0.487396 0.144232 0.012604];
theta = pi * fractional_index / half_width_of_BLIT;
blit_impulse = aa * cos((0:3)' * theta) ...
  .* (abs(fractional_index) <= half_width_of_BLIT);
end
