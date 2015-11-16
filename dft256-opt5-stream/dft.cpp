#include "dft.h"
#include "coefficients256.h"
#include <ap_int.h>
#include <hls_stream.h>

static void dft_inner(int k, DTYPE real_sample[SIZE], DTYPE imag_sample[SIZE], hls::stream<DTYPE> &outreal, hls::stream<DTYPE> &outimag)
{
	DTYPE sumreal = 0;
	DTYPE sumimag = 0;
    ap_uint<8> angle = 0;
	for (int t = 0; t < SIZE; ++t, angle += k) {
		const DTYPE cos_angle = cos_coefficients_table[angle];
		const DTYPE sin_angle = -sin_coefficients_table[angle];
		sumreal +=  real_sample[t] * cos_angle + imag_sample[t] * sin_angle;
		sumimag += -real_sample[t] * sin_angle + imag_sample[t] * cos_angle;
	}
    outreal.write(sumreal);
    outreal.write(sumimag);
}

void dft(DTYPE real_sample[SIZE], DTYPE imag_sample[SIZE], hls::stream<DTYPE> &outreal, hls::stream<DTYPE> &outimag)
{
    for (int k = 0; k < SIZE; ++k) {
    	dft_inner(k,real_sample,imag_sample,outreal[k],outimag[k]);
    }
}
