function AudioWrite(fname, y, fs)
% wrapper function of audiowrite for Octave compatibiliyty
%
% AudioWrite(fname, y, fs)
%

narginchk(3, 3);
if ~isOctave
  audiowrite(fname, y, fs);
else
  wavwrite(y, fs, fname);
end;
end