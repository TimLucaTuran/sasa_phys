function fOut = EmbedExtractCentral2D(field, len, fillValue)
%UNTITLED1 Summary of this function goes here
%  Detailed explanation goes here

    if (nargin < 3)
        fillValue = 0; 
    end

    if CHECK_INPUT_PARAMTERS
        if (ndims(field) ~= length(len))
            error('field size and length dimension must match!');
        end

        if (ndims(field) ~= 2)
            error('field dimension must be two!');
        end
    end
    
    [leftIndex, rightIndex] = embedExtractRegionCalculation( len(1), size(field,1) );
    if (size(field,1) <= len(1))
        fOut = fillValue * ones(len(1), size(field,2));
        fOut(leftIndex : rightIndex, :) = field;
    else
        fOut = field(leftIndex : rightIndex, :);
    end
    
    field = fOut;
    
    [leftIndex, rightIndex] = embedExtractRegionCalculation( len(2), size(field,2) );
    if (size(field,2) <= len(2))
        fOut = fillValue * ones(size(field,1), len(2));
        fOut(:, leftIndex : rightIndex) = field;
    else
        fOut = field(:, leftIndex : rightIndex);
    end    
    
%endfunction