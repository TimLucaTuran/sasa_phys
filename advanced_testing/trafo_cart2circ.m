function CIRC_OUT = trafo_cart2circ(CART_IN)
% coordinate trafo from Cartesian to circular basis

% create cell of index-cell-vectors
% the cell-vectors denote s-matrix slices of the entire parameter space
DIM_IN = ndims(CART_IN);
II = cell(4,4);
for loop1 = 1:4
    for loop2 = 1:4
        II{loop1,loop2} = smat_phys_ind_arb_sz(DIM_IN,loop1,loop2);
    end
end

% front transmission
TFXX = CART_IN( II{1,1}{:} );
TFXY = CART_IN( II{1,2}{:} );
TFYX = CART_IN( II{2,1}{:} );
TFYY = CART_IN( II{2,2}{:} );

TFRR = TFXX + TFYY + 1j*(TFXY - TFYX);
TFRL = TFXX - TFYY - 1j*(TFXY + TFYX);
TFLR = TFXX - TFYY + 1j*(TFXY + TFYX);
TFLL = TFXX + TFYY - 1j*(TFXY - TFYX);

% back transmissiom
TBXX = CART_IN( II{3,3}{:} );
TBXY = CART_IN( II{3,4}{:} );
TBYX = CART_IN( II{4,3}{:} );
TBYY = CART_IN( II{4,4}{:} );

TBRR = TBXX + TBYY + 1j*(TBXY - TBYX);
TBRL = TBXX - TBYY - 1j*(TBXY + TBYX);
TBLR = TBXX - TBYY + 1j*(TBXY + TBYX);
TBLL = TBXX + TBYY - 1j*(TBXY - TBYX);

% front reflection
RFXX = CART_IN( II{3,1}{:} );
RFXY = CART_IN( II{3,2}{:} );
RFYX = CART_IN( II{4,1}{:} );
RFYY = CART_IN( II{4,2}{:} );

RFRR = RFXX + RFYY + 1j*(RFXY - RFYX);
RFRL = RFXX - RFYY - 1j*(RFXY + RFYX);
RFLR = RFXX - RFYY + 1j*(RFXY + RFYX);
RFLL = RFXX + RFYY - 1j*(RFXY - RFYX);

% back reflection
RBXX = CART_IN( II{1,3}{:} );
RBXY = CART_IN( II{1,4}{:} );
RBYX = CART_IN( II{2,3}{:} );
RBYY = CART_IN( II{2,4}{:} );

RBRR = RBXX + RBYY + 1j*(RBXY - RBYX);
RBRL = RBXX - RBYY - 1j*(RBXY + RBYX);
RBLR = RBXX - RBYY + 1j*(RBXY + RBYX);
RBLL = RBXX + RBYY - 1j*(RBXY - RBYX);

% put new matrix 2gether
% T-front block
CIRC_OUT( II{1,1}{:} ) = TFRR;
CIRC_OUT( II{1,2}{:} ) = TFRL;
CIRC_OUT( II{2,1}{:} ) = TFLR;
CIRC_OUT( II{2,2}{:} ) = TFLL;

% T-back block
CIRC_OUT( II{3,3}{:} ) = TBRR;
CIRC_OUT( II{3,4}{:} ) = TBRL;
CIRC_OUT( II{4,3}{:} ) = TBLR;
CIRC_OUT( II{4,4}{:} ) = TBLL;

% R-back block
CIRC_OUT( II{1,3}{:} ) = RBRR;
CIRC_OUT( II{1,4}{:} ) = RBRL;
CIRC_OUT( II{2,3}{:} ) = RBLR;
CIRC_OUT( II{2,4}{:} ) = RBLL;

% R-front block
CIRC_OUT( II{3,1}{:} ) = RFRR;
CIRC_OUT( II{3,2}{:} ) = RFRL;
CIRC_OUT( II{4,1}{:} ) = RFLR;
CIRC_OUT( II{4,2}{:} ) = RFLL;

% multiply by norm from lambda-matrices (sqrt(2)^-1).^2
CIRC_OUT = 0.5*CIRC_OUT;