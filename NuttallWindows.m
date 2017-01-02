function w = NuttallWindows(window_length, window_index)
%   Nuttall window with continuous first order derivative
%   No.11 and No.12 items in [1] are selectable
%
%   w = NuttallWindows(N, window_index)
%
% Argument
%   window_length : length of window (samples)
%   window_index : item number in the table of the reference
%
% Return value
%   w : window coefficients

%   Reference:
%   [1] Nuttall, A.; ``Some windows with very good sidelobe behavior,''
%   Acoustics, Speech and Signal Processing, IEEE Transactions on ,
%   vol.29, no.1, pp. 84--91, Feb 1981.
%   doi: 10.1109/TASSP.1981.1163506

w = [];
x = (0:window_length - 1)' * 2 * pi / (window_length - 1);
switch window_index
  case 11
    aa = [0.338946 0.481973 0.161054 0.018027];
  case 12
    aa = [0.355768 0.487396 0.144232 0.012604];
end;
bb = [aa(1); -aa(2); aa(3); -aa(4)];
w = cos(x * [0 1 2 3]) * bb;
end
