function is_octave = isOctave
% check if the environment is Octave
% is_octave = isOctave
%
% Return value
%   is_octave: result, true: Octave, false: Matlab

if isempty(strfind(version,'R20'))
  is_octave = true;
else
  is_octave = false;
end;
end
