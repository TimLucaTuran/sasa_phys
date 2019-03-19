function [ leftIndex, rightIndex ] = embedExtractRegionCalculation( len, size )
%EMBEDEXTRACTREGIONCALCULATION Summary of this function goes here
%   Detailed explanation goes here

    if CHECK_INPUT_PARAMTERS
        if ((dims(len) ~= 0) || (dims(size) ~= 0))
            error('len & size must be scalar!');
        end

        if (len < 0 || size < 0)
            error('len & size must be positive scalars!');
        end
    end

    if (size <= len)
        leftIndex = fix((len - size) / 2);        
        if ((rem(len,2) == 0) && (rem(size,2) ~= 0))
            % if len is even and size is odd ...
            leftIndex = leftIndex + 1;
        end
        
        rightIndex = leftIndex + size;
        leftIndex = leftIndex+1;
    else
        leftIndex = fix((size - len) / 2);
        if ((rem(len,2) ~= 0) && (rem(size,2) == 0))
            % if len is even and size is odd ...
            leftIndex = leftIndex + 1;
        end

        rightIndex = leftIndex + len;
        leftIndex = leftIndex+1;
    end
    
% endfunction