function indexOut = calcFieldIndexOfCentralField( fieldLength, sp )
%transforms the sampling coordinate (pixel coordinate) of a centraled field
%in the numerical field coordinate --> by shifting of fieldLength/2.

    if CHECK_INPUT_PARAMTERS && (dims(fieldLength) > 1)
        error('fieldLength dimension cannot exceed 1!');
    end
    
    spWasNoCell = false;
    if (~iscell(sp))
        spWasNoCell = true;
        if (dims(fieldLength) > 0)
            error('if field dimension is bigger than 1, sp must be a cell array!');
        else
            sp = {sp};
        end
    end

    if CHECK_INPUT_PARAMTERS
        if (dims(sp) > 1)
            error('sp must be a 1D cell array!');
        end

        if (length(sp) ~= length(fieldLength))
            error('length of sp and fieldLength must match!');
        end
    end
    
    index = {};
    for n = 1:length(sp)
        index{n} = sp{n} + (ones(size(sp{n})) * (-fix(-fieldLength(n)/2) + 1));
    end
       
    if (spWasNoCell)
        indexOut = index{1};
    else
        indexOut = index;
    end   

%endfunction