#include "dft.h"
#include "coefficients256.h"

void dft_inner(int k, DTYPE real_sample[SIZE], DTYPE imag_sample[SIZE], DTYPE &outreal, DTYPE &outimag)
{
	outreal = 0;
	outimag = 0;
	unsigned char angle = 0;
	for (int t = 0; t < SIZE; ++t, angle += k) {
		const DTYPE cos_angle = cos_coefficients_table[angle];
		const DTYPE sin_angle = -sin_coefficients_table[angle];
		outreal +=  real_sample[t] * cos_angle + imag_sample[t] * sin_angle;
		outimag += -real_sample[t] * sin_angle + imag_sample[t] * cos_angle;
	}
}

void dft(DTYPE real_sample[SIZE], DTYPE imag_sample[SIZE], DTYPE outreal[SIZE], DTYPE outimag[SIZE])
{
    for (int k = 0; k < SIZE; ++k) {
    	dft_inner(k,real_sample,imag_sample,outreal[k],outimag[k]);
    }
}
