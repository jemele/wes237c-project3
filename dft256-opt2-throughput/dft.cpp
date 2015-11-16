#include "dft.h"
#include "coefficients256.h"
#include <ap_int.h>

void dft_inner(int k, DTYPE real_sample[SIZE], DTYPE imag_sample[SIZE], DTYPE &outreal, DTYPE &outimag)
{
#pragma HLS unroll
	outreal = 0;
	outimag = 0;
	ap_uint<8> angle = 0;
	for (int t = 0; t < SIZE; ++t, angle += k) {
#pragma HLS unroll
		const DTYPE cos_angle = cos_coefficients_table[angle];
		const DTYPE sin_angle = -sin_coefficients_table[angle];
		outreal +=  real_sample[t] * cos_angle + imag_sample[t] * sin_angle;
		outimag += -real_sample[t] * sin_angle + imag_sample[t] * cos_angle;
	}
}

void dft(DTYPE real_sample[SIZE], DTYPE imag_sample[SIZE], DTYPE outreal[SIZE], DTYPE outimag[SIZE])
{
#pragma HLS dataflow
    for (int k = 0; k < SIZE; ++k) {
#pragma HLS pipeline
    	dft_inner(k,real_sample,imag_sample,outreal[k],outimag[k]);
    }
}
