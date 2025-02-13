fft: fft.c
	gcc -g fft.c -o fft -lm  # -lm links the math library

clean:
	rm -f *.o *.d $(BIN)