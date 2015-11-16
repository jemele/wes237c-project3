#include "dft.h"
#include "coefficients256.h"
#include <ap_int.h>

void dft_inner(int k, const DTYPE real_sample[SIZE], const DTYPE imag_sample[SIZE], DTYPE &outreal, DTYPE &outimag)
{
#pragma HLS UNROLL
	outreal = 0;
	outimag = 0;
	ap_uint<8> angle = 0;
	for (int t = 0; t < SIZE; ++t, angle += k) {
		const DTYPE cos_angle = cos_coefficients_table[angle];
		const DTYPE sin_angle = -sin_coefficients_table[angle];
		outreal +=  real_sample[t] * cos_angle + imag_sample[t] * sin_angle;
		outimag += -real_sample[t] * sin_angle + imag_sample[t] * cos_angle;
	}
}

void dft(DTYPE real_sample[SIZE], DTYPE imag_sample[SIZE], DTYPE outreal[SIZE], DTYPE outimag[SIZE])
{
#pragma HLS ARRAY_PARTITION variable=imag_sample complete dim=1
#pragma HLS ARRAY_PARTITION variable=real_sample complete dim=1
#pragma HLS DATAFLOW
	const int offset1 = SIZE/4;
	const int offset2 = 2*offset1;
	const int offset3 = 3*offset1;
    for (int k = 0; k < offset1; ++k) {
    	dft_inner(k,real_sample,imag_sample,outreal[k],outimag[k]);
    	dft_inner(k+offset1,real_sample,imag_sample,outreal[k+offset1],outimag[k+offset1]);
    	dft_inner(k+offset2,real_sample,imag_sample,outreal[k+offset2],outimag[k+offset2]);
    	dft_inner(k+offset3,real_sample,imag_sample,outreal[k+offset3],outimag[k+offset3]);

    }
}
