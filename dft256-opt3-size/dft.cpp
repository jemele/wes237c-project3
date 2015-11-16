#include "dft.h"
#include "coefficients256.h"
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
    DTYPE sumreal = 0;
    DTYPE sumimag = 0;
    ap_uint<8>angle = 0;
    for (int t = 0; t < SIZE; ++t, angle += k) {
        const DTYPE cos_angle = cos_coefficients_table[angle];
        const DTYPE sin_angle = sin_coefficients_table[angle];
        const DTYPE real_part = sub(mult(real_sample[t],cos_angle),mult(imag_sample[t],sin_angle));
        sumreal = add(real_part,sumreal);
        const DTYPE imag_part = add(mult(real_sample[t],sin_angle),mult(imag_sample[t],cos_angle));
        sumimag = add(imag_part,sumimag);
    }
    outreal = sumreal;
    outimag = sumimag;
}

void dft(DTYPE real_sample[SIZE], DTYPE imag_sample[SIZE], DTYPE outreal[SIZE], DTYPE outimag[SIZE])
{
    for (int k = 0; k < SIZE; ++k) {
        dft_inner(k,real_sample,imag_sample,outreal[k],outimag[k]);
    }
}
