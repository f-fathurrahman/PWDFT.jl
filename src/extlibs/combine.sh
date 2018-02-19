gcc -c -fPIC -O3 fft3d.c

#gcc -shared -o fft3d.so -Wl,--whole-archive fft3d.o libfftw3.a libblas.so.3 -Wl,--no-whole-archive
gcc -shared -o fft3d.so -Wl,--whole-archive fft3d.o libfftw3.so libblas.so.3 -Wl,--no-whole-archive
#gcc -shared -o fft3d_julia_0.4.5.so -Wl,--whole-archive fft3d.o libfftw3.so.3.4.4 -Wl,--no-whole-archive

#clang -c -march=native -fPIC -O3 fft3d.c
#clang -shared -o fft3d_ffr.so -Wl,--whole-archive fft3d.o libfftw3_ffr.a -Wl,--no-whole-archive
