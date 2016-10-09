My test result with 1 trial and 100 trails:
```bash
$ ./bin/mmul 600 8  
                     Min        Mean         Max
      OpenMP: 1.5229e-01  1.5443e-01  1.5800e-01
  ThreadPool: 1.5181e-01  1.5746e-01  1.6226e-01
LFThreadPool: 1.4992e-01  1.5297e-01  1.5830e-01
threadpool11: 1.4851e-01  1.4981e-01  1.5172e-01

$ ./bin/mmul 600 8 10 
                     Min        Mean         Max
      OpenMP: 1.5093e-01  1.5313e-01  1.5747e-01
  ThreadPool: 1.5226e-01  1.5751e-01  1.6271e-01
LFThreadPool: 1.4959e-01  1.5235e-01  1.5608e-01
threadpool11: 1.4789e-01  1.5357e-01  1.5820e-01

$ ./bin/mmul 600 8 100
                     Min        Mean         Max
      OpenMP: 1.4166e-01  1.4565e-01  1.5344e-01
  ThreadPool: 1.4173e-01  1.5490e-01  2.8225e-01
LFThreadPool: 1.3958e-01  1.4553e-01  1.5349e-01
threadpool11: 1.3895e-01  1.5309e-01  1.9100e-01

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
