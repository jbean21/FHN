# include <fftw3.h>
# include <cmath>

// link with library using -lfftw3f -lm (single double precision) or -lfftw3f -lm (long double precision) on Unix

// -lfftw3_omp -lfftw3 -lm  <------ OpenMP compiler
// -fopenmp

// real (in[i][0]) and imaginary (in[i][1])

// http://www.fftw.org/fftw3.pdf

// INSTALLATION:
// ./configure --enable-openmp && make
// sudo make install

// Reference: Matteo  Frigo  and  Steven  G.  Johnson,  “The  design and implementation of FFTW3,” \textit{Proc. IEEE} \textbf{93} (2), 216–231 (2005)



int main(void) {
  fftw_complex *in, *out;
  fftw_plan p;

  
  int N = 5;

  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  
  
  in[0][0] = 1/sqrt(2*3.141592);
  in[0][1] = 0;
  
  in[1][0] = 0.24197;
  in[1][1] = 0;
  
  in[2][0] = 0.14676;
  in[2][1] = 0;
  
  in[3][0] = 0.08902;
  in[3][1] = 0;
  
  in[4][0] = 0.05399;
  in[4][1] = 0;
  
  fftw_execute(p); /* repeat as needed */
  
  printf("%f + %f i => %f + %f i\n", in[0][0], in[0][1], out[0][0], out[0][1]);
  printf("%f + %f i => %f + %f i\n", in[1][0], in[1][1], out[1][0], out[1][1]);
  printf("%f + %f i => %f + %f i\n", in[2][0], in[2][1], out[2][0], out[2][1]);
  printf("%f + %f i => %f + %f i\n", in[3][0], in[3][1], out[3][0], out[3][1]);
  printf("%f + %f i => %f + %f i\n", in[4][0], in[4][1], out[4][0], out[4][1]);
  
  /* clean up */
  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);

}
