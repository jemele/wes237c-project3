#include <hls_stream.h>
typedef float DTYPE;
#define SIZE 256 		/* SIZE OF DFT */
typedef hls::stream<DTYPE> stream_t;
void dft(DTYPE real_sample[SIZE], DTYPE imag_sample[SIZE], stream_t &outreal, stream_t &outimag);

