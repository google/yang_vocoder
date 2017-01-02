function [x, fs] = AudioRead(fname)
% wrapper function of audioread for Octave compatibiliyty
%
% [x, fs] = AudioRead(fname)
%

narginchk(1, 1);
if ~isOctave
  [x, fs] = audioread(fname);
else
  [x, fs] = wavread(fname);
end;
end
