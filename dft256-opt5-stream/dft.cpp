#include "dft.h"
#include "coefficients256.h"
#include <ap_int.h>
#include <hls_stream.h>

static void dft_inner(int k, DTYPE real_samples[SIZE], DTYPE imag_samples[SIZE], DTYPE &outreal, DTYPE &outimag)
{
    DTYPE sumreal = 0;
    DTYPE sumimag = 0;
    ap_uint<8>angle = 0;
    for (int t = 0; t < SIZE; ++t, angle += k) {
        const DTYPE cos_angle = cos_coefficients_table[angle];
        const DTYPE sin_angle = sin_coefficients_table[angle];
        const DTYPE real_sample = real_samples[t];
        const DTYPE imag_sample = imag_samples[t];
        sumreal +=  real_sample * cos_angle + imag_sample * sin_angle;
        sumimag += -real_sample * sin_angle + imag_sample * cos_angle;
    }
    outreal = sumreal;
    outimag = sumimag;
}

void dft(DTYPE real_samples[SIZE], DTYPE imag_samples[SIZE], stream_t &outreal,stream_t &outimag)
{
#pragma HLS dataflow
#pragma HLS INTERFACE axis port=outreal bundle=OUTPUT_STREAM
#pragma HLS INTERFACE axis port=outimag bundle=OUTPUT_STREAM
    for (int k = 0; k < SIZE; ++k) {
        DTYPE oreal, oimag;
        dft_inner(k,real_samples,imag_samples,oreal,oimag);
        outreal.write(oreal);
        outimag.write(oimag);
    }
}
