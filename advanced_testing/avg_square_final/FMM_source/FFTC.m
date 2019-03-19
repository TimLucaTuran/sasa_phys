function [ ftOut, dkx ] = FFTC( field, dx, dim )
%Performs a fft under the assumption that the origin of the field lies on 
%the central pixel coordinate.   

    %% the matlab integrated dft-function named fft performs an fft whereby
    %% the "real" field-zero-coordinate and the numerical zero-field coordinate
    %% coinside. So the procedure is the following: 1) shift the "real"
    %% field coordinate that it matches the numerical field coordinate
    %% (done by: fftshift). 2) perform a fft. 3) undo the shift, so that
    %% the real field origin lies again in the middle of the numerical
    %% field.
    
    global ft;
    ft = field;
    s = size(ft);
    
    if (nargin < 3)
    % if no specific dimension is given to transform to only ...                

        dim = 1:ndims(ft);
        if CHECK_INPUT_PARAMTERS && (dims(ft) ~= length(dx))
            error('Dimensions of Field and SamplingDistances must match');
        end        

        if (dims(ft) == 1)
            fftshift(dim);
            ft = fft(ft);
            invfftshift(dim);
            ft = ft ./ sqrt(prod(s));
            dkx = 2*pi / (max(s)*dx);
        elseif (dims(ft) == 2)
            fftshift(dim);
            ft = fft2(ft);
            invfftshift(dim);
            ft = ft ./ sqrt(prod(s));
            dkx = 2*pi ./ (s.*dx);
        else
            fftshift(dim);
            ft = fftn(ft);
            invfftshift(dim);
            ft = ft ./ sqrt(prod(s));
            dkx = 2*pi ./ (s.*dx);        
        end        
    else
        if (max(dim) > length(s))
            s = [s, ones(1, max(dim) - length(s))];
        end
        
        dkx = zeros(size(dim));
        fftshift(dim);
        for d = 1:length(dim)
            ft = fft(ft, [], dim(d));
            ft = ft ./ sqrt(s(dim(d)));
            dkx(d) = 2*pi / (s(dim(d)) * dx(d));
        end
        invfftshift(dim);
    end
    
    ftOut = ft;
    clear global ft;    

%endfunction

function fftshift(dim)

    global ft;

    if (ndims(ft) <= 2)
        shift(1) = calcFieldIndexOfCentralField(size(ft,1), 0) - 1;
        shift(2) = calcFieldIndexOfCentralField(size(ft,2), 0) - 1;
    else
        shift = zeros(1, ndims(ft));
        for d = 1:ndims(ft)
            shift(d) = calcFieldIndexOfCentralField(size(ft, d), 0) - 1;
        end
    end
    
    dummy = 1:ndims(ft);
    dummy(dim(dim <= ndims(ft))) = [];
        
    shift(dummy) = 0;

    ft = circshift(ft, -shift);

%endfunction

function invfftshift(dim)

    global ft;

    if (ndims(ft) <= 2)
        shift(1) = calcFieldIndexOfCentralField(size(ft,1), 0) - 1;
        shift(2) = calcFieldIndexOfCentralField(size(ft,2), 0) - 1;
    else
        shift = zeros(1, ndims(ft));
        for d = 1:ndims(ft)
            shift(d) = calcFieldIndexOfCentralField(size(ft, d), 0) - 1;
        end
    end
    
    dummy = 1:ndims(ft);
    dummy(dim(dim <= ndims(ft))) = [];
    shift(dummy) = 0;

    ft = circshift(ft, shift);

%endfunction