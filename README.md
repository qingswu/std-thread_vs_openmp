# std::thread vs OpenMP
A performance comparison of a threadpool built from C++11's std::thread library and the shared-memory parallel programming API OpenMP.

## Matrix-matrix multiplication test
The code `mmul` (from `mmul.cc` and `dummy.cc`) multiplies two square matrices using both OpenMP and the ThreadPool. The usage is:
```bash
mmul [MATRIX_SIZE] [NUM_THREADS] [NUM_TRIALS]
```
All parameters are optional.

An example from my Macbook Pro:
```bash
$ ./mmul 600 4
                   Min        Mean         Max
    OpenMP: 7.8162e-01  8.0546e-01  8.2313e-01
ThreadPool: 7.7258e-01  7.9581e-01  8.4442e-01
```
