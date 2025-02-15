# FFT Example in C

This repository contains an example implementation of the Fast Fourier Transform (FFT) in C using the Cooley-Tukey Radix-2 Decimation-in-Time (DIT) algorithm. The code demonstrates how to initialize FFT buffers, apply windowing functions, compute the FFT, convert the complex output to magnitudes, and detect the major frequency peak in the resulting spectrum.

## Features

- **FFT Initialization:** Sets up the real and imaginary buffers for FFT computation.
- **Windowing Functions:** Applies window functions (Rectangular and Hamming) to reduce spectral leakage.
- **FFT Computation:** Uses the Cooley-Tukey Radix-2 DIT algorithm for efficient FFT computation.
- **Magnitude Calculation:** Converts the complex FFT output into magnitude values.
- **Peak Detection:** Identifies and prints frequency peaks from the FFT result.

## Code Structure

- **Configuration Macros:**
  - `FFT_SIZE`: Number of samples for FFT (must be a power of 2).
  - `SAMPLE_RATE`: Sampling frequency in Hz.
  
- **Buffers:**
  - `inputBuffer`: Array holding the input sample data.
  - `real` and `imag`: Arrays storing the real and imaginary parts of the FFT data.

- **Function Overview:**
  - `FFT_Init(float *vReal, float *vImag, uint16_t samples, float samplingFrequency)`: Initializes the FFT buffers.
  - `FFT_Windowing(uint8_t windowType, uint8_t dir)`: Applies a windowing function to the input data to reduce spectral leakage.
  - `Compute(uint8_t dir)`: Computes the FFT using the Cooley-Tukey Radix-2 DIT algorithm.
  - `ComplexToMagnitude(void)`: Converts the complex FFT output to magnitude (spectrum).
  - `MajorPeak()`: Identifies the frequency of the largest magnitude peak in the FFT result.
  - `main()`: Generates test data, performs the FFT, and prints the frequency spectrum.

## Compilation

To compile the code, use a C compiler such as GCC. For example:

```bash
gcc -o fft_example main.c -lm