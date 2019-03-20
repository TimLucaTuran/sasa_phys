import numpy as np
import cython
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.



cdef void star_product_core(complex[:,:,:,:] s_in1,
                            complex[:,:,:,:] s_in2,
                            complex[:,:,:,:] s_out,
                            int H, int L):

    cdef complex TF1XX, TF1XY, RB1XX, RB1XY, TF1YX, TF1YY, RB1YX, RB1YY, RF1XX
    cdef complex RF1XY, TB1XX, TB1XY, RF1YX, RF1YY, TB1YX, TB1YY
    cdef complex TF2XX, TF2XY, RB2XX, RB2XY, TF2YX, TF2YY, RB2YX, RB2YY, RF2XX
    cdef complex RF2XY, TB2XX, TB2XY, RF2YX, RF2YY, TB2YX, TB2YY

    cdef Py_ssize_t h, l

    for h in range(H):
        for l in range(L):

            # S-matrix 1
            TF1XX = s_in1[h,l,0,0]
            TF1XY = s_in1[h,l,0,1]
            RB1XX = s_in1[h,l,0,2]
            RB1XY = s_in1[h,l,0,3]
            TF1YX = s_in1[h,l,1,0]
            TF1YY = s_in1[h,l,1,1]
            RB1YX = s_in1[h,l,1,2]
            RB1YY = s_in1[h,l,1,3]
            RF1XX = s_in1[h,l,2,0]
            RF1XY = s_in1[h,l,2,1]
            TB1XX = s_in1[h,l,2,2]
            TB1XY = s_in1[h,l,2,3]
            RF1YX = s_in1[h,l,3,0]
            RF1YY = s_in1[h,l,3,1]
            TB1YX = s_in1[h,l,3,2]
            TB1YY = s_in1[h,l,3,3]

            # S-matrix 2
            TF2XX = s_in2[h,l,0,0]
            TF2XY = s_in2[h,l,0,1]
            RB2XX = s_in2[h,l,0,2]
            RB2XY = s_in2[h,l,0,3]
            TF2YX = s_in2[h,l,1,0]
            TF2YY = s_in2[h,l,1,1]
            RB2YX = s_in2[h,l,1,2]
            RB2YY = s_in2[h,l,1,3]
            RF2XX = s_in2[h,l,2,0]
            RF2XY = s_in2[h,l,2,1]
            TB2XX = s_in2[h,l,2,2]
            TB2XY = s_in2[h,l,2,3]
            RF2YX = s_in2[h,l,3,0]
            RF2YY = s_in2[h,l,3,1]
            TB2YX = s_in2[h,l,3,2]
            TB2YY = s_in2[h,l,3,3]



            # Plain analytic form of the staproduct
            TFXX = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
                RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
                RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(TF1YX*(RB1XX*RF2XY*TF2XX+
                RB1XY*RF2YY*TF2XX+TF2XY+(-1)*RB1XX*RF2XX*TF2XY+(-1)*RB1XY*
                RF2YX*TF2XY)+TF1XX*(TF2XX+(-1)*RB1YX*RF2XY*TF2XX+(-1)*
                RB1YY*RF2YY*TF2XX+RB1YX*RF2XX*TF2XY+RB1YY*RF2YX*TF2XY))
            # -------------------------------------------------------------------------
            TFXY = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
                RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
                RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(TF1YY*(RB1XX*RF2XY*TF2XX+
                RB1XY*RF2YY*TF2XX+TF2XY+(-1)*RB1XX*RF2XX*TF2XY+(-1)*RB1XY*
                RF2YX*TF2XY)+TF1XY*(TF2XX+(-1)*RB1YX*RF2XY*TF2XX+(-1)*
                RB1YY*RF2YY*TF2XX+RB1YX*RF2XX*TF2XY+RB1YY*RF2YX*TF2XY))
            # -------------------------------------------------------------------------
            TFYX = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
                RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
                RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(TF1YX*(RB1XX*RF2XY*TF2YX+
                RB1XY*RF2YY*TF2YX+TF2YY+(-1)*RB1XX*RF2XX*TF2YY+(-1)*RB1XY*
                RF2YX*TF2YY)+TF1XX*(TF2YX+(-1)*RB1YX*RF2XY*TF2YX+(-1)*
                RB1YY*RF2YY*TF2YX+RB1YX*RF2XX*TF2YY+RB1YY*RF2YX*TF2YY))
            # -------------------------------------------------------------------------
            TFYY = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
                RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
                RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(TF1YY*(RB1XX*RF2XY*TF2YX+
                RB1XY*RF2YY*TF2YX+TF2YY+(-1)*RB1XX*RF2XX*TF2YY+(-1)*RB1XY*
                RF2YX*TF2YY)+TF1XY*(TF2YX+(-1)*RB1YX*RF2XY*TF2YX+(-1)*
                RB1YY*RF2YY*TF2YX+RB1YX*RF2XX*TF2YY+RB1YY*RF2YX*TF2YY))
            # -------------------------------------------------------------------------
            RBXX = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
                RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
                RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(RB2XX*(((-1)+RB1YX*RF2XY)*
                ((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+RB1YY*RF2XY*RF2YX)+
                RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+RB1XY*RB1YX*RF2XX)*
                RF2YY)+RB1XY*RB1YX*RF2YY*TB2XX*TF2XX+RB1XY*TB2YX*TF2XX+(-1)
                *RB1XY*RB1YX*RF2XY*TB2YX*TF2XX+RB1YX*TB2XX*TF2XY+(-1)*
                RB1XY*RB1YX*RF2YX*TB2XX*TF2XY+RB1YY*TB2YX*TF2XY+RB1XY*
                RB1YX*RF2XX*TB2YX*TF2XY+RB1XX*(RB1YY*TB2YX*(RF2XY*TF2XX+(
                -1)*RF2XX*TF2XY)+TB2XX*(TF2XX+(-1)*RB1YY*RF2YY*TF2XX+RB1YY*
                RF2YX*TF2XY)))
            # -------------------------------------------------------------------------
            RBXY = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
                RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
                RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(RB2XY*(((-1)+RB1YX*RF2XY)*
                ((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+RB1YY*RF2XY*RF2YX)+
                RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+RB1XY*RB1YX*RF2XX)*
                RF2YY)+RB1XY*RB1YX*RF2YY*TB2XY*TF2XX+RB1XY*TB2YY*TF2XX+(-1)
                *RB1XY*RB1YX*RF2XY*TB2YY*TF2XX+RB1YX*TB2XY*TF2XY+(-1)*
                RB1XY*RB1YX*RF2YX*TB2XY*TF2XY+RB1YY*TB2YY*TF2XY+RB1XY*
                RB1YX*RF2XX*TB2YY*TF2XY+RB1XX*(RB1YY*TB2YY*(RF2XY*TF2XX+(
                -1)*RF2XX*TF2XY)+TB2XY*(TF2XX+(-1)*RB1YY*RF2YY*TF2XX+RB1YY*
                RF2YX*TF2XY)))
            # -------------------------------------------------------------------------
            RBYX = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
                RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
                RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(RB2YX*(((-1)+RB1YX*RF2XY)*
                ((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+RB1YY*RF2XY*RF2YX)+
                RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+RB1XY*RB1YX*RF2XX)*
                RF2YY)+RB1XY*RB1YX*RF2YY*TB2XX*TF2YX+RB1XY*TB2YX*TF2YX+(-1)
                *RB1XY*RB1YX*RF2XY*TB2YX*TF2YX+RB1YX*TB2XX*TF2YY+(-1)*
                RB1XY*RB1YX*RF2YX*TB2XX*TF2YY+RB1YY*TB2YX*TF2YY+RB1XY*
                RB1YX*RF2XX*TB2YX*TF2YY+RB1XX*(RB1YY*TB2YX*(RF2XY*TF2YX+(
                -1)*RF2XX*TF2YY)+TB2XX*(TF2YX+(-1)*RB1YY*RF2YY*TF2YX+RB1YY*
                RF2YX*TF2YY)))
            # -------------------------------------------------------------------------
            RBYY = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
                RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
                RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(RB2YY*(((-1)+RB1YX*RF2XY)*
                ((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+RB1YY*RF2XY*RF2YX)+
                RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+RB1XY*RB1YX*RF2XX)*
                RF2YY)+RB1XY*RB1YX*RF2YY*TB2XY*TF2YX+RB1XY*TB2YY*TF2YX+(-1)
                *RB1XY*RB1YX*RF2XY*TB2YY*TF2YX+RB1YX*TB2XY*TF2YY+(-1)*
                RB1XY*RB1YX*RF2YX*TB2XY*TF2YY+RB1YY*TB2YY*TF2YY+RB1XY*
                RB1YX*RF2XX*TB2YY*TF2YY+RB1XX*(RB1YY*TB2YY*(RF2XY*TF2YX+(
                -1)*RF2XX*TF2YY)+TB2XY*(TF2YX+(-1)*RB1YY*RF2YY*TF2YX+RB1YY*
                RF2YX*TF2YY)))
            # -------------------------------------------------------------------------
            RFXX = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
                RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
                RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(RF1XX*(((-1)+RB1YX*RF2XY)*
                ((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+RB1YY*RF2XY*RF2YX)+
                RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+RB1XY*RB1YX*RF2XX)*
                RF2YY)+RB1YY*RF2XY*RF2YX*TB1XX*TF1XX+RF2YX*TB1XY*TF1XX+(-1)
                *RB1YX*RF2XY*RF2YX*TB1XY*TF1XX+RF2XY*TB1XX*TF1YX+(-1)*
                RB1XY*RF2XY*RF2YX*TB1XX*TF1YX+RB1XX*RF2XY*RF2YX*TB1XY*
                TF1YX+RF2YY*TB1XY*TF1YX+RF2XX*(RF2YY*TB1XY*(RB1YX*TF1XX+(-1)
                *RB1XX*TF1YX)+TB1XX*(TF1XX+(-1)*RB1YY*RF2YY*TF1XX+RB1XY*
                RF2YY*TF1YX)))
            # -------------------------------------------------------------------------
            RFXY = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
                RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
                RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(RF1XY*(((-1)+RB1YX*RF2XY)*
                ((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+RB1YY*RF2XY*RF2YX)+
                RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+RB1XY*RB1YX*RF2XX)*
                RF2YY)+RB1YY*RF2XY*RF2YX*TB1XX*TF1XY+RF2YX*TB1XY*TF1XY+(-1)
                *RB1YX*RF2XY*RF2YX*TB1XY*TF1XY+RF2XY*TB1XX*TF1YY+(-1)*
                RB1XY*RF2XY*RF2YX*TB1XX*TF1YY+RB1XX*RF2XY*RF2YX*TB1XY*
                TF1YY+RF2YY*TB1XY*TF1YY+RF2XX*(RF2YY*TB1XY*(RB1YX*TF1XY+(-1)
                *RB1XX*TF1YY)+TB1XX*(TF1XY+(-1)*RB1YY*RF2YY*TF1XY+RB1XY*
                RF2YY*TF1YY)))
            # -------------------------------------------------------------------------
            RFYX = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
                RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
                RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(RF1YX*(((-1)+RB1YX*RF2XY)*
                ((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+RB1YY*RF2XY*RF2YX)+
                RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+RB1XY*RB1YX*RF2XX)*
                RF2YY)+RB1YY*RF2XY*RF2YX*TB1YX*TF1XX+RF2YX*TB1YY*TF1XX+(-1)
                *RB1YX*RF2XY*RF2YX*TB1YY*TF1XX+RF2XY*TB1YX*TF1YX+(-1)*
                RB1XY*RF2XY*RF2YX*TB1YX*TF1YX+RB1XX*RF2XY*RF2YX*TB1YY*
                TF1YX+RF2YY*TB1YY*TF1YX+RF2XX*(RF2YY*TB1YY*(RB1YX*TF1XX+(-1)
                *RB1XX*TF1YX)+TB1YX*(TF1XX+(-1)*RB1YY*RF2YY*TF1XX+RB1XY*
                RF2YY*TF1YX)))
            # -------------------------------------------------------------------------
            RFYY = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
                RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
                RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(RF1YY*(((-1)+RB1YX*RF2XY)*
                ((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+RB1YY*RF2XY*RF2YX)+
                RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+RB1XY*RB1YX*RF2XX)*
                RF2YY)+RB1YY*RF2XY*RF2YX*TB1YX*TF1XY+RF2YX*TB1YY*TF1XY+(-1)
                *RB1YX*RF2XY*RF2YX*TB1YY*TF1XY+RF2XY*TB1YX*TF1YY+(-1)*
                RB1XY*RF2XY*RF2YX*TB1YX*TF1YY+RB1XX*RF2XY*RF2YX*TB1YY*
                TF1YY+RF2YY*TB1YY*TF1YY+RF2XX*(RF2YY*TB1YY*(RB1YX*TF1XY+(-1)
                *RB1XX*TF1YY)+TB1YX*(TF1XY+(-1)*RB1YY*RF2YY*TF1XY+RB1XY*
                RF2YY*TF1YY)))
            # -------------------------------------------------------------------------
            TBXX = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
                RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
                RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(TB1XY*(RB1XX*RF2YX*TB2XX+
                RB1YX*RF2YY*TB2XX+TB2YX+(-1)*RB1XX*RF2XX*TB2YX+(-1)*RB1YX*
                RF2XY*TB2YX)+TB1XX*(TB2XX+(-1)*RB1XY*RF2YX*TB2XX+(-1)*
                RB1YY*RF2YY*TB2XX+RB1XY*RF2XX*TB2YX+RB1YY*RF2XY*TB2YX))
            # -------------------------------------------------------------------------
            TBXY = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
                RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
                RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(TB1XY*(RB1XX*RF2YX*TB2XY+
                RB1YX*RF2YY*TB2XY+TB2YY+(-1)*RB1XX*RF2XX*TB2YY+(-1)*RB1YX*
                RF2XY*TB2YY)+TB1XX*(TB2XY+(-1)*RB1XY*RF2YX*TB2XY+(-1)*
                RB1YY*RF2YY*TB2XY+RB1XY*RF2XX*TB2YY+RB1YY*RF2XY*TB2YY))
            # -------------------------------------------------------------------------
            TBYX = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
                RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
                RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(TB1YY*(RB1XX*RF2YX*TB2XX+
                RB1YX*RF2YY*TB2XX+TB2YX+(-1)*RB1XX*RF2XX*TB2YX+(-1)*RB1YX*
                RF2XY*TB2YX)+TB1YX*(TB2XX+(-1)*RB1XY*RF2YX*TB2XX+(-1)*
                RB1YY*RF2YY*TB2XX+RB1XY*RF2XX*TB2YX+RB1YY*RF2XY*TB2YX))
            # -------------------------------------------------------------------------
            TBYY = (((-1)+RB1YX*RF2XY)*((-1)+RB1XY*RF2YX)+(-1)*RB1XX*(RF2XX+
                RB1YY*RF2XY*RF2YX)+RB1XX*RB1YY*RF2XX*RF2YY+(-1)*(RB1YY+
                RB1XY*RB1YX*RF2XX)*RF2YY)**(-1)*(TB1YY*(RB1XX*RF2YX*TB2XY+
                RB1YX*RF2YY*TB2XY+TB2YY+(-1)*RB1XX*RF2XX*TB2YY+(-1)*RB1YX*
                RF2XY*TB2YY)+TB1YX*(TB2XY+(-1)*RB1XY*RF2YX*TB2XY+(-1)*
                RB1YY*RF2YY*TB2XY+RB1XY*RF2XX*TB2YY+RB1YY*RF2XY*TB2YY))


            # Assemble the resulting s-matrix using the elements from above


            s_out[h,l,0,0] = TFXX
            s_out[h,l,0,1] = TFXY
            s_out[h,l,0,2] = RBXX
            s_out[h,l,0,3] = RBXY
            s_out[h,l,1,0] = TFYX
            s_out[h,l,1,1] = TFYY
            s_out[h,l,1,2] = RBYX
            s_out[h,l,1,3] = RBYY
            s_out[h,l,2,0] = RFXX
            s_out[h,l,2,1] = RFXY
            s_out[h,l,2,2] = TBXX
            s_out[h,l,2,3] = TBXY
            s_out[h,l,3,0] = RFYX
            s_out[h,l,3,1] = RFYY
            s_out[h,l,3,2] = TBYX
            s_out[h,l,3,3] = TBYY
    return






def star_product(complex[:,:,:,:] s_in1,complex[:,:,:,:] s_in2):
    cdef int H1 = s_in1.shape[0]
    cdef int H2 = s_in2.shape[0]
    cdef int H = max(H1, H2)
    cdef int L = s_in1.shape[1]

    s_out = np.zeros((H, L, 4, 4), dtype=complex)
    #cdef complex[:,:,:,:] s_out_view = s_out
    if H != 1:
        if H1 == 1:
            s_in1 = np.repeat(s_in1, H, axis=0)
        elif H2 == 1:
            s_in2 = np.repeat(s_in2, H, axis=0)
    star_product_core(s_in1, s_in2, s_out, H, L)
    return s_out
