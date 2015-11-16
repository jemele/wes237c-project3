#include "dft.h"
#include "coefficients1024.h"
#include <ap_int.h>

static DTYPE mult(DTYPE a, DTYPE b)
{
    return a*b;
}

static DTYPE add(DTYPE a, DTYPE b)
{
    return a+b;
}

static DTYPE sub(DTYPE a, DTYPE b)
{
    return a-b;
}

static void dft_inner(int k, DTYPE real_sample[SIZE], DTYPE imag_sample[SIZE], DTYPE &outreal, DTYPE &outimag)
{
    outreal = 0;
    outimag = 0;
    ap_uint<10>angle = 0;
    for (int t = 0; t < SIZE; ++t, angle += k) {
        const DTYPE cos_angle = cos_coefficients_table[angle];
        const DTYPE sin_angle = sin_coefficients_table[angle];

        const DTYPE real_part = sub(mult(real_sample[t],cos_angle),mult(imag_sample[t],sin_angle));
        outreal = add(real_part,outreal);
        const DTYPE imag_part = add(mult(real_sample[t],sin_angle),mult(imag_sample[t],cos_angle));
        outimag = add(imag_part,outimag);
    }
}

void dft(DTYPE real_sample[SIZE], DTYPE imag_sample[SIZE], DTYPE outreal[SIZE], DTYPE outimag[SIZE])
{
#pragma HLS DATAFLOW

    for (int k = 0; k < SIZE; ++k) {
        dft_inner(k,real_sample,imag_sample,outreal[k],outimag[k]);
    }
}
