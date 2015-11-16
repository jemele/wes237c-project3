#include <ap_fixed.h>
typedef ap_fixed<35,15> DTYPE;
#define SIZE 256
void dft(DTYPE real_sample[SIZE], DTYPE imag_sample[SIZE], DTYPE outreal[SIZE], DTYPE outimag[SIZE]);

