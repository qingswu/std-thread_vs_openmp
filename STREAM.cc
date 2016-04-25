#include <iostream>

#include "ThreadPool.hh"
#include "Timer.hh"
#include "omp.h"

using namespace std;


int main (int argc, char** argv) {


  int N = 1e7;
  int NUM_THREADS = 4;
  if (argc>1) { NUM_THREADS = atoi(argv[1]); }
  if (argc>2) { N = atoi(argv[2]); }


  auto a = (double*)calloc(N,sizeof(double));
  auto b = (double*)calloc(N,sizeof(double));
  auto c = (double*)calloc(N,sizeof(double));

  for (int i=0; i<N; i++) { a[i] = 0; b[i] = 1; c[i] = 2; }

  int ntrials = 10;
  double tperformance = 0.0;


  auto benchmark = [&](auto callable) {

    ThreadPool pool(NUM_THREADS);

// -----------------------------------------------------------------------------------------------
    cout << "single thread performance\n";

    for (int j=0; j<N; j++) {
      callable(j,a,b,c);
    }

    tperformance = 0.0;
    for (int i=0; i<ntrials; i++)
    {
      Timer timer([&](int elapsed){
          //cout << "Trial " << i << ": "<< elapsed*1e-6 << " ms\n";
          tperformance+=elapsed;
      });
      for (int j=0; j<N; j++) {
        callable(j,a,b,c);
      }
    }
    cout << "Average: " << tperformance*1e-6 / ntrials << " ms\n\n";

// -----------------------------------------------------------------------------------------------
    cout << "std::thread ThreadPool performance\n";

    // cold start for timing purposes
    pool.ParallelFor(0, N, callable, a, b, c);

    tperformance = 0.0;
    for (int i=0; i<ntrials; i++)
    {
      Timer timer([&](int elapsed){
          //cout << "Trial " << i << ": "<< elapsed*1e-6 << " ms\n";
          tperformance+=elapsed;
      });
      pool.ParallelFor(0, N, callable, a, b, c);
    }
    cout << "Average: " << tperformance*1e-6 / ntrials << " ms\n\n";

// -----------------------------------------------------------------------------------------------
    cout << "OpenMP performance\n";

    // omp cold start for timing purposes
    omp_set_dynamic(0);

    #pragma omp parallel for num_threads(NUM_THREADS)
    for (int j=0; j<N; j++) {
      callable(j,a,b,c);
    }


    tperformance = 0.0;
    for (int i=0; i<ntrials; i++)
    {
      Timer timer([&](int elapsed){
          //cout << "Trial " << i << ": "<< elapsed*1e-6 << " ms\n";
          tperformance+=elapsed;
      });

      #pragma omp parallel for num_threads(NUM_THREADS)
      for (int j=0; j<N; j++) {
        callable(j,a,b,c);
      }
    }
    cout << "Average: " << tperformance*1e-6 / ntrials << " ms\n\n";

// -----------------------------------------------------------------------------------------------

  };

  cout<< endl << " - - - - - - - - - | COPY | - - - - - - - - - " << endl << endl;
  benchmark([](int k, double* _a, double* _b, double* _c) {
    return _a[k] = _b[k];
  });

  cout<< endl << " - - - - - - - - - | SCALE | - - - - - - - - - " << endl << endl;
  benchmark([](int k, double* _a, double* _b, double* _c) {
    return _a[k] = 4*_b[k];
  });

  cout<< endl << " - - - - - - - - - | ADD | - - - - - - - - - " << endl << endl;
  benchmark([](int k, double* _a, double* _b, double* _c) {
      return _a[k] = _b[k] + _c[k];
  });

  cout<< endl << " - - - - - - - - - | TRIAD | - - - - - - - - - " << endl << endl;
  benchmark([](int k, double* _a, double* _b, double* _c) {
      return _a[k] = _b[k] + 4*_c[k];
  });


  return 0;
}
