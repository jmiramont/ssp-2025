#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <fftw3.h>
#include <algorithm>      // Include max and min.
#define PI 3.14159265358979323846

// Compile with:
// g++ -shared -o src/stft_wrap/libstft.so -fPIC src/stft_wrap/stft.cpp -lfftw3 -lm


extern "C" {

// Find local minima in the spectrogram magnitude
void find_spectrogram_minima(const std::vector<float>& magnitude_spectrogram, 
                            int rows, 
                            int cols, 
                            int* minima_x, 
                            int* minima_y, 
                            int* minima_count) {
    int count = 0;

    // std::cout << "Hola_3";
    // std::vector<int> x(rows) ;
    // std::vector<int> y(rows) ;

    // Iterate through each element, skipping the first and last rows/cols
    for (int i = 1; i < rows - 1; ++i) {
        for (int j = 1; j < cols - 1; ++j) {
            double val = magnitude_spectrogram[i * cols + j];

            // Check against 8 neighbors
            if (val <= magnitude_spectrogram[(i - 1) * cols + j] &&
                val <= magnitude_spectrogram[(i + 1) * cols + j] &&
                val <= magnitude_spectrogram[i * cols + (j - 1)] &&
                val <= magnitude_spectrogram[i * cols + (j + 1)] &&
                val <= magnitude_spectrogram[(i - 1) * cols + (j - 1)] &&
                val <= magnitude_spectrogram[(i + 1) * cols + (j + 1)] &&
                val <= magnitude_spectrogram[(i - 1) * cols + (j + 1)] &&
                val <= magnitude_spectrogram[(i + 1) * cols + (j - 1)]) {
                
                minima_x[count] = i;
                minima_y[count] = j;
                count++;
            }
        }
    }
    *minima_count = count;

}


void round_gauss_win(double *window, int N, float L){
    float tt = 0;
    int t0 = (N-1)/2+1;

    for (int n = 0; n < N; ++n) {
        tt = n-t0;
        float exp_arg = pow(tt/L,2);
        window[n] = exp(-exp_arg*PI);
    }
}

// Compute STFT
void compute_spectrogram_and_zeros(const double* signal, 
                int signal_length, 
                int nfft, 
                int* zeros_x,
                int* zeros_y,
                int* zeros_count
                ) {
    // std::cout << "Hola \n";
    std::vector<float> spect_out((nfft/2+1)*(signal_length/2), 0.0f);

    float prec = 1e-15;
    double L = sqrt(nfft);
    int l = floor(sqrt(-nfft*log(prec)/PI))+1;
    int N = 2*l+1;
    float t0 = l+1;

    std::vector<double> window(N);
    round_gauss_win(window.data(), N, L);
    
    fftw_complex *in, *out;
    fftw_plan p;

    // Allocate memory for the outputs of the FFT (just one row)
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nfft);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nfft);
    // fftw_plan_with_nthreads(1); // Force single-threaded FFT
    p = fftw_plan_dft_1d(nfft, in, out, FFTW_FORWARD, FFTW_PATIENT);

    for (int i = signal_length/4-t0; i < 3*signal_length/4-t0; i++) {
        for (int j = 0; j < nfft; j++) {
            if (j<N){
                in[j][0] = signal[i+j]*window[j];
            }
            else{
                in[j][0] = 0.0;
            }
            in[j][1] = 0.0; 
        }

        fftw_execute(p); // Compute FFT

        for (int j = 0; j < nfft/2+1; j++) { // Spectrogram rows.
            spect_out[(i-signal_length/4+int(t0))*(nfft/2+1)+j] = float(pow(out[j][0],2) + pow(out[j][1],2));
        }
    }

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    int rows = signal_length/2;
    int cols = nfft/2+1; 
    
    find_spectrogram_minima(spect_out, 
                            rows, 
                            cols, 
                            zeros_x, 
                            zeros_y,
                            zeros_count);
    
}

// Compute ISTFT
// Calcular STFT > Enmascarar > Invertir
void compute_istft(const double* signal, 
                int signal_length, 
                int nfft, 
                bool* mask,
                double* reconstructed_signal) 
                   {

    float prec = 1e-15;
    double L = sqrt(nfft);
    int l = floor(sqrt(-nfft*log(prec)/PI))+1;
    int N = 2*l+1;
    int t0 = l+1;

    std::vector<double> window(N);
    round_gauss_win(window.data(), N, L);
    std::vector<double> window_sum(signal_length, 0.0);

    fftw_complex *in, *out, *mid;
    fftw_plan p;
    fftw_plan pinv;

    // Allocate memory for the outputs of the FFT (just one row)
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nfft);
    mid = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nfft);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nfft);

    p = fftw_plan_dft_1d(nfft, in, mid, FFTW_FORWARD, FFTW_ESTIMATE);
    pinv = fftw_plan_dft_1d(nfft, mid, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    for (int i = signal_length/4-t0; i < 3*signal_length/4+t0; i++){

        for (int j = 0; j < nfft; j++) {
            if (j<N){
                in[j][0] = signal[i+j]*window[j];
            }
            else{
                in[j][0] = 0.0;
            }
            in[j][1] = 0.0; 
        }

        // std::cout << "FFT \n";

        fftw_execute(p); // Compute FFT
         
        // Now for the inversion
        // std::cout << "mid to in2 \n";
        // Here we should modify mid with the mask
        if (i < 3*signal_length/4-t0){
        for (int j = 0; j < nfft/2+1; j++) {
            mid[j][0] = mid[j][0]*mask[(i-signal_length/4+t0)*(nfft/2+1)+j];
            mid[j][1] = mid[j][1]*mask[(i-signal_length/4+t0)*(nfft/2+1)+j];
        }
        for (int j=nfft-1; j>= nfft/2+1; j--) {
            mid[j][0] = mid[j][0]*mask[(i-signal_length/4+t0)*(nfft/2+1)+nfft-1-j];
            mid[j][1] = mid[j][1]*mask[(i-signal_length/4+t0)*(nfft/2+1)+nfft-1-j];
        }
        }

        // std::cout << "Inversion \n";
        fftw_execute(pinv); // Compute IFFT

        // std::cout << "Rearranging Inverse \n";
        for (int j = 0; j < N; j++) {
            double value = out[j][0] / nfft;
            reconstructed_signal[i-signal_length/4+t0 + j] += value;
            window_sum[i-signal_length/4+t0 + j] += window[j];
        }
    }

    // std::cout << "Normalizing \n";
    for (int i = 0; i < signal_length; i++) {
        if (window_sum[i] > 0) {
            reconstructed_signal[i] /= window_sum[i];
        }
    }

    // std::cout << "Chau \n";
    fftw_destroy_plan(p);
    fftw_destroy_plan(pinv);
    fftw_free(in);
    fftw_free(mid);
    fftw_free(out);
    
}


}