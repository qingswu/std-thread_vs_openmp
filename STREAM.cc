#include <iostream>

#include "ThreadPool.hh"
#include "Timer.hh"
#include "threadpool11/threadpool11.hpp"
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
    cout << "std::thread ThreadPool performance\n";

    {
      ThreadPool pool(NUM_THREADS);

      // i==0 is a cold start for timing purposes
      tperformance = 0.0;
      for (int i=0; i<ntrials+1; i++)
      {
        Timer timer([&](int elapsed){ if(i>0) tperformance+=elapsed; });
        pool.ParallelFor(0, N, callable, a, b, c);
      }
      cout << "Average: " << tperformance*1e-6 / ntrials << " ms\n\n";
    }


// -----------------------------------------------------------------------------------------------
    {
      threadpool11::Pool pool(NUM_THREADS);
      cout << "threadpool11 performance\n";

      tperformance = 0.0;
      for (int i=0; i<ntrials+1; i++)
      {
        // i==0 is a cold start for timing purposes
        Timer timer([&](int elapsed){ if(i>0) tperformance+=elapsed; });

        //std::array<std::future<void>, NUM_THREADS> futures;
        std::future<void>* futures = new std::future<void>[NUM_THREADS];
        auto begin = 0; auto end = N;
        int chunk = (end - begin) / NUM_THREADS;
        for (int j = 0; j < NUM_THREADS; ++j) {
          futures[j] = pool.postWork<void>([=]() {
              uint32_t threadstart = begin + j*chunk;
              uint32_t threadstop = (j == NUM_THREADS - 1) ? end : threadstart + chunk;
              for (uint32_t it = threadstart; it < threadstop; ++it) {
                callable(it,a,b,c);
              }
            });
        }
        for (int j = 0; j < NUM_THREADS; ++j) {
          futures[j].get();
        }
      }
      cout << "Average: " << tperformance*1e-6 / ntrials << " ms\n\n";
    }

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
