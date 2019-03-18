function ret = isOctave()
%ISOCTAVE Summary of this function goes here
%   Detailed explanation goes here

    ret = true;
    try
        OCTAVE_VERSION;
    catch
        ret = false;
    end
end
