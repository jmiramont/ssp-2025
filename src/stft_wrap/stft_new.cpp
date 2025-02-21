#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <fftw3.h>
#include <algorithm>      // Include max and min.
#define PI 3.14159265358979323846

// Find local minima in the spectrogram magnitude
void find_spectrogram_minima(const std::vector<float>& magnitude_spectrogram,
                             int rows, int cols,
                             std::vector<int>& minima_x,
                             std::vector<int>& minima_y) {
    for (int i = 1; i < rows - 1; ++i) {
        for (int j = 1; j < cols - 1; ++j) {
            float val = magnitude_spectrogram[i * cols + j];
            if (val <= magnitude_spectrogram[(i - 1) * cols + j] &&
                val <= magnitude_spectrogram[(i + 1) * cols + j] &&
                val <= magnitude_spectrogram[i * cols + (j - 1)] &&
                val <= magnitude_spectrogram[i * cols + (j + 1)]) {
                minima_x.push_back(i);
                minima_y.push_back(j);
            }
        }
    }
}

extern "C" {

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
void compute_spectrogram_and_zeros(const double* signal, int signal_length, int nfft,
                                   std::vector<int>& zeros_x, std::vector<int>& zeros_y) {
    int rows = signal_length / 2;
    int cols = nfft / 2 + 1;
    std::vector<float> spect_out(rows * cols, 0.0f);
    
    std::vector<double> window(nfft);
    fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nfft);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nfft);
    fftw_plan p = fftw_plan_dft_1d(nfft, in, out, FFTW_FORWARD, FFTW_PATIENT);

    for (int i = signal_length / 4; i < 3 * signal_length / 4; i++) {
        for (int j = 0; j < nfft; j++) {
            in[j][0] = (j < window.size()) ? signal[i + j] * window[j] : 0.0;
            in[j][1] = 0.0;
        }
        fftw_execute(p);
        for (int j = 0; j < cols; j++) {
            spect_out[(i - signal_length / 4) * cols + j] = float(out[j][0] * out[j][0] + out[j][1] * out[j][1]);
        }
    }

    find_spectrogram_minima(spect_out, rows, cols, zeros_x, zeros_y);
    
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
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
        for (int j = nfft/2+1; j < nfft; j++) {
            mid[j][0] = mid[j][0]*mask[(i-signal_length/4+t0)*(nfft/2+1)+j-nfft/2];
            mid[j][1] = mid[j][1]*mask[(i-signal_length/4+t0)*(nfft/2+1)+j-nfft/2];
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