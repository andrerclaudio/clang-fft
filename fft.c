#ifndef _DEFAULT_SOURCE
#define _DEFAULT_SOURCE 1
#endif

#include <time.h>

// FFT Configuration
#define FFT_SIZE 256     // Number of samples for FFT. MUST be a power of 2.
#define SAMPLE_RATE 2048 // Sampling frequency in Hz

#define COOLEY_TUKEY
#define ATSAM
#define GOERTZEL

// Array to store the input signal data for FFT processing.
double inputSignal[FFT_SIZE];
float inputBuffer[FFT_SIZE];

#ifdef COOLEY_TUKEY

#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

// Output buffers: real and imaginary parts of the FFT
float real[FFT_SIZE];
float imag[FFT_SIZE];

// Function Prototypes
void FFT_Init(float *vReal, float *vImag, uint16_t samples, float samplingFrequency);
void FFT_Windowing(uint8_t windowType, uint8_t dir);
void Compute(uint8_t dir);
void ComplexToMagnitude(void);
float MajorPeak();

/**
 * @brief Initializes the FFT buffers.
 *
 * This function initializes the real and imaginary parts of the FFT.
 * The imaginary part is initialized to zero.
 *
 * @param vReal             Pointer to the real part of the FFT buffer. This will contain the input samples and the real part of the FFT after processing.
 * @param vImag             Pointer to the imaginary part of the FFT buffer. This should be zero before the FFT and will contain the imaginary part of the FFT after processing.
 * @param samples           Number of samples (FFT_SIZE). Must be a power of 2.
 * @param samplingFrequency Sampling frequency of the signal. Should be at least twice the maximum signal frequency (Nyquist rate).
 */
void FFT_Init(float *vReal, float *vImag, uint16_t samples, float samplingFrequency)
{
    // Initialize imaginary array
    for (int i = 0; i < samples; i++)
    {
        vImag[i] = 0.0f; // Set imaginary part to 0
    }
}

/**
 * @brief Applies a windowing function to the input data to reduce spectral leakage.
 *
 * This function applies a windowing function to the real part of the input buffer before FFT computation.
 * Common window types include rectangular, Hamming, Hanning, and Blackman windows [1].
 *
 * @param windowType Type of window to apply.
 *                   0: Rectangular window
 *                   1: Hamming window
 *                   Other values: Defaults to Rectangular window
 * @param dir       Direction flag (currently only processes when dir == 0).
 */
void FFT_Windowing(uint8_t windowType, uint8_t dir)
{
    float window[FFT_SIZE]; // Array to store window values

    // Generate window values based on the specified window type
    switch (windowType)
    {
    default: // Default: Rectangular window
    case 0:  // FFT_WIN_TYP_RECTANGLE: Rectangular window
        for (int i = 0; i < FFT_SIZE; i++)
        {
            window[i] = 1.0; // All values are 1
        }
        break;

    case 1: // FFT_WIN_TYP_HAMMING: Hamming window
        for (int i = 0; i < FFT_SIZE; i++)
        {
            window[i] = 0.54 - 0.46 * cos(2 * M_PI * i / (FFT_SIZE - 1)); // Hamming window formula [1]
        }
        break;
    }

    // Apply the window to the real part of the input buffer
    if (dir == 0)
    { // Apply window
        for (int i = 0; i < FFT_SIZE; i++)
        {
            real[i] = real[i] * window[i]; // Multiply each sample by the corresponding window value
        }
    }
    else
    {
        // inverse direction
    }
}

/**
 * @brief Computes the FFT (Fast Fourier Transform) using the Cooley-Tukey Radix-2 DIT algorithm.
 *
 * This function performs an in-place FFT on the real and imaginary arrays.
 *
 * @param dir Direction of the FFT.
 *            0: Forward FFT (time domain to frequency domain)
 *            1: Inverse FFT (frequency domain to time domain)
 */
void Compute(uint8_t dir)
{
    const int n = FFT_SIZE;

    // Bit reversal permutation:  Rearrange data in bit-reversed order
    for (int i = 0; i < n; i++)
    {
        int j = 0;
        int k = i;
        for (int t = 0; t < (int)(log2(n)); t++)
        {
            j = (j << 1) | (k & 1); // Shift j left and add the least significant bit of k
            k >>= 1;                // Shift k right
        }
        if (j > i)
        {
            // Swap real parts
            float tempReal = real[i];
            real[i] = real[j];
            real[j] = tempReal;

            // Swap imaginary parts
            float tempImag = imag[i];
            imag[i] = imag[j];
            imag[j] = tempImag;
        }
    }

    // Cooley-Tukey iterative process
    for (int len = 2; len <= n; len <<= 1) // len is the length of the sub-FFTs
    {
        double angle = 2 * M_PI / len; // Calculate angle increment
        if (dir)
            angle *= -1; // Reverse angle for inverse FFT

        for (int i = 0; i < n; i += len) // i is the starting index of each sub-FFT
        {
            for (int j = 0; j < len / 2; j++) // j iterates within each sub-FFT
            {
                double wReal = cos(angle * j); // Real part of twiddle factor
                double wImag = sin(angle * j); // Imaginary part of twiddle factor

                float uReal = real[i + j];           // Real part of the first element in the butterfly
                float uImag = imag[i + j];           // Imaginary part of the first element in the butterfly
                float tReal = real[i + j + len / 2]; // Real part of the second element in the butterfly
                float tImag = imag[i + j + len / 2]; // Imaginary part of the second element in the butterfly

                float twiddleReal = (tReal * wReal) - (tImag * wImag); // Real part of the twiddle product
                float twiddleImag = (tReal * wImag) + (tImag * wReal); // Imaginary part of the twiddle product

                // Butterfly computation
                real[i + j] = uReal + twiddleReal;
                imag[i + j] = uImag + twiddleImag;
                real[i + j + len / 2] = uReal - twiddleReal;
                imag[i + j + len / 2] = uImag - twiddleImag;
            }
        }
    }

    // Normalize if performing inverse FFT
    if (dir)
    {
        for (int i = 0; i < n; i++)
        {
            real[i] /= n;
            imag[i] /= n;
        }
    }
}

/**
 * @brief Converts the complex FFT output to magnitude (spectrum).
 *
 * This function calculates the magnitude of each complex number in the FFT output.
 * The magnitude is calculated as sqrt(real[i]^2 + imag[i]^2).
 */
void ComplexToMagnitude(void)
{
    for (int i = 0; i < FFT_SIZE; i++)
    {
        real[i] = sqrt(real[i] * real[i] + imag[i] * imag[i]); // Calculate magnitude
    }
}

/**
 * @brief Finds the frequency of the largest magnitude peak in the FFT result.
 *
 * This function searches for the index of the largest magnitude in the real array (which now holds magnitude values).
 * It then converts this index to a frequency using the formula: frequency = index * sampleRate / FFT_SIZE [3].
 *
 * @return The frequency (in Hz) of the largest peak.
 */
float MajorPeak()
{
    int peakIndex = 0;   // Index of the largest peak
    float peakValue = 0; // Magnitude of the largest peak

    // Find the index of the largest magnitude in the first half of the FFT result
    for (int i = 1; i < FFT_SIZE / 2; i++)
    {
        if (real[i] > peakValue)
        {
            peakValue = real[i]; // Update peak value
            peakIndex = i;       // Update peak index
        }
    }

    // Convert the index to frequency
    return (float)peakIndex * SAMPLE_RATE / FFT_SIZE; // frequency = index * sampleRate / FFT_SIZE [3]
}

void cooleyTukey()
{
    // Fill input buffer with sample data (replace with your actual data)
    for (int i = 0; i < FFT_SIZE; i++)
    {
        real[i] = inputBuffer[i]; // Copy input data to real part of FFT input
        imag[i] = 0;              // Initialize imaginary part to 0
    }

    // Initialize FFT
    FFT_Init(real, imag, FFT_SIZE, SAMPLE_RATE);

    // Apply windowing function
    FFT_Windowing(0, 0); // Rectangular window (windowType = 0)

    // Compute FFT
    Compute(0); // Forward FFT (dir = 0)

    // Calculate magnitudes
    ComplexToMagnitude();

    float threshold = 1.0; // Threshold for printing peaks (adjust as needed)

    for (int i = 1; i < FFT_SIZE / 2; i++)
    {
        float frequency = (float)i * SAMPLE_RATE / FFT_SIZE; // Calculate frequency for this bin
        float magnitude = real[i];                           // Get magnitude for this bin

        // Simple peak detection: check if this magnitude is greater than its neighbors AND above the threshold
        if (i > 1 && i < (FFT_SIZE / 2) - 1 && magnitude > real[i - 1] && magnitude > real[i + 1] && magnitude > threshold)
        {
            printf("Magnitude at %.2f Hz: %.2f\n", frequency, magnitude);
        }
    }
}

#endif

#ifdef ATSAM

#include <math.h>
#include <stdio.h>

// Complex number structure for FFT processing and output.
typedef struct
{
    double real; ///< Real part of the complex number.
    double imag; ///< Imaginary part of the complex number.
} Complex;

// Output buffer for FFT results.
Complex fftOutput[FFT_SIZE];

// Buffer for storing magnitude of FFT results (only the first half is needed due to symmetry).
double magnitude[FFT_SIZE / 2];

// Buffer for corresponding frequency values for each FFT bin.
double frequency[FFT_SIZE / 2];

/**
 * @brief Compute the Fast Fourier Transform (FFT) of a complex input array.
 *
 * This function recursively computes the FFT of an array of Complex numbers.
 * It divides the input into even and odd indexed elements, computes their FFT,
 * and then combines the results using the butterfly operation.
 *
 * @param v Pointer to the array of Complex numbers. On input, it contains time-domain samples.
 *          On output, it holds the FFT result.
 * @param n Number of samples in the array. Must be a power of 2.
 * @param tmp Pointer to a temporary buffer array of Complex numbers of size n.
 */
void fft(Complex *v, int n, Complex *tmp)
{
    if (n > 1)
    {
        int k, m;
        Complex z, w, *vo, *ve;
        // Split input array into even and odd indexed elements.
        ve = tmp;         // Even-indexed elements.
        vo = tmp + n / 2; // Odd-indexed elements.
        for (k = 0; k < n / 2; k++)
        {
            ve[k] = v[2 * k];     // Even indices.
            vo[k] = v[2 * k + 1]; // Odd indices.
        }
        // Recursively compute FFT on even and odd parts.
        fft(ve, n / 2, v);
        fft(vo, n / 2, v);
        // Combine the two halves using the butterfly operation.
        for (m = 0; m < n / 2; m++)
        {
            // Compute the twiddle factor: e^(-j*2*pi*m/n)
            w.real = cos(2 * M_PI * m / (double)n);
            w.imag = -sin(2 * M_PI * m / (double)n);
            // Multiply the odd part by the twiddle factor.
            z.real = w.real * vo[m].real - w.imag * vo[m].imag;
            z.imag = w.real * vo[m].imag + w.imag * vo[m].real;
            // Combine even and odd parts to form the FFT output.
            v[m].real = ve[m].real + z.real;
            v[m].imag = ve[m].imag + z.imag;
            v[m + n / 2].real = ve[m].real - z.real;
            v[m + n / 2].imag = ve[m].imag - z.imag;
        }
    }
}

void atsam(void)
{
    // Declare arrays for FFT processing: 'v' for the complex signal and 'tmp' as temporary storage.
    Complex v[FFT_SIZE], tmp[FFT_SIZE];

    // Generate a test signal with three frequency components: 16 Hz, 200 Hz, and 1000 Hz.
    // (Note: The comment mentions 10 Hz, but the code uses 16 Hz for the first component.)
    for (int i = 0; i < FFT_SIZE; i++)
    {
        // Initialize the complex array 'v' using the input signal.
        v[i].real = inputSignal[i]; // Real part from the input signal.
        v[i].imag = 0.0;            // Imaginary part is set to zero.
    }

    // Perform the FFT on the complex signal.
    fft(v, FFT_SIZE, tmp);

    // Copy the FFT results to the global output buffer.
    for (int i = 0; i < FFT_SIZE; i++)
    {
        fftOutput[i] = v[i];
    }

    // Calculate the magnitude for each FFT bin (only the first half is needed due to symmetry).
    for (int i = 0; i < FFT_SIZE / 2; i++)
    {
        magnitude[i] = sqrt(fftOutput[i].real * fftOutput[i].real +
                            fftOutput[i].imag * fftOutput[i].imag);
    }

    // Compute the frequency corresponding to each FFT bin.
    for (int i = 0; i < FFT_SIZE / 2; i++)
    {
        frequency[i] = (SAMPLE_RATE / FFT_SIZE) * i;
    }

    // Detect and print peaks in the magnitude spectrum.
    // A peak is defined as a point that is greater than its immediate neighbors and above a set threshold.
    for (int i = 1; i < FFT_SIZE / 2; i++)
    {
        if (i > 1 && i < (FFT_SIZE / 2) - 1 &&
            magnitude[i] > magnitude[i - 1] &&
            magnitude[i] > magnitude[i + 1] &&
            magnitude[i] > 1.0)
        {
            printf("Magnitude at %.2f Hz: %.2f\n", frequency[i], magnitude[i]);
        }
    }
}

#endif

#ifdef GOERTZEL

#include <math.h>
#include <stdint.h>
#include <stdio.h>

/**
 * @brief Computes the magnitude of a target frequency using the Goertzel algorithm.
 *
 * @param samples         The input signal array.
 * @param numSamples      Number of samples in the input signal.
 * @param targetFrequency The frequency (in Hz) to analyze.
 * @param sampleRate      The sampling frequency of the signal.
 * @return float          Magnitude of the target frequency component.
 */
float compute(const float samples[], int numSamples, float targetFrequency, float sampleRate)
{
    double s_prev = 0.0;
    double s_prev2 = 0.0;
    // Compute the normalized frequency (radians per sample)
    double omega = (2.0 * M_PI * targetFrequency) / sampleRate;
    // Precompute the coefficient used in the recurrence
    double coeff = 2.0 * cos(omega);

    // Process all samples
    for (int i = 0; i < numSamples; i++)
    {
        double s = samples[i] + coeff * s_prev - s_prev2;
        s_prev2 = s_prev;
        s_prev = s;
    }

    // Compute the real and imaginary parts of the result
    double real_part = s_prev - s_prev2 * cos(omega);
    double imag_part = s_prev2 * sin(omega);

    // Return the magnitude of the frequency component
    return (float)sqrt(real_part * real_part + imag_part * imag_part);
}

void Goertzel(void)
{
    // Choose a target frequency to analyze with Goertzel
    float targetFrequency = 16.0;
    float magnitude = compute(inputBuffer, FFT_SIZE, targetFrequency, SAMPLE_RATE);
    printf("Magnitude at %.2f Hz: %.2f\n", targetFrequency, magnitude);

    targetFrequency = 200.0;
    magnitude = compute(inputBuffer, FFT_SIZE, targetFrequency, SAMPLE_RATE);
    printf("Magnitude at %.2f Hz: %.2f\n", targetFrequency, magnitude);

    targetFrequency = 1000.0;
    magnitude = compute(inputBuffer, FFT_SIZE, targetFrequency, SAMPLE_RATE);
    printf("Magnitude at %.2f Hz: %.2f\n", targetFrequency, magnitude);
}

#endif

int main(void)
{
    clock_t start, end;
    double elapsed;

    // Fill input buffer with sample data (replace with your actual data)
    for (int i = 0; i < FFT_SIZE; i++)
    {
        // Create a test signal with three frequency components: 10 Hz, 200 Hz, and 1000 Hz
        inputBuffer[i] = (0.50 * sin(2 * M_PI * 16 * i / SAMPLE_RATE)) +
                         (0.20 * sin(2 * M_PI * 200 * i / SAMPLE_RATE)) +
                         (0.30 * sin(2 * M_PI * 1000 * i / SAMPLE_RATE));

        inputSignal[i] = inputBuffer[i];
    }

#ifdef COOLEY_TUKEY
    start = clock();
    printf("Cooley Tukey Algorithm\n");
    cooleyTukey();
    end = clock();
    elapsed = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Time taken: %f seconds\n\n", elapsed);
#endif

#ifdef ATSAM
    start = clock();
    printf("Microchip ATSAM Algorithm\n");
    atsam();
    end = clock();
    elapsed = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Time taken: %f seconds\n\n", elapsed);
#endif

#ifdef GOERTZEL
    start = clock();
    printf("Goertzel Algorithm\n");
    Goertzel();
    end = clock();
    elapsed = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Time taken: %f seconds\n\n", elapsed);
#endif

    return 0;
}
