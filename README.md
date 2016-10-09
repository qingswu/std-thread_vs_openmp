My test result with 1 trial and 100 trails:
```bash
$ ./bin/mmul 600 4
                     Min        Mean         Max
      OpenMP: 2.6492e-01  2.7967e-01  2.9745e-01
  ThreadPool: 2.5470e-01  2.7414e-01  2.7992e-01
LFThreadPool: 1.8939e-01  1.9907e-01  2.1628e-01
threadpool11: 1.8301e-01  1.8313e-01  1.8338e-01

$ ./bin/mmul 600 8 100
                     Min        Mean         Max
      OpenMP: 1.5175e-01  1.5457e-01  1.6374e-01
  ThreadPool: 1.5155e-01  1.6049e-01  1.8564e-01
LFThreadPool: 1.4836e-01  1.5173e-01  1.5802e-01
threadpool11: 1.4787e-01  1.5236e-01  1.6303e-01

```
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
