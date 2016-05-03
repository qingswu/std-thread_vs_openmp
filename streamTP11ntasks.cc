#include <iostream>
#include <iomanip>

#include "Timer.hh"
#include "threadpool11/threadpool11.hpp"

using namespace std;


int main (int argc, char** argv) {

  int N = 1e7;
  int NUM_THREADS = 4;
  int ntasks = 0;
  if (argc>1) { NUM_THREADS = atoi(argv[1]); }
  if (argc>2) { N = atoi(argv[2]); }
  if (argc>3) { ntasks = atoi(argv[3]); }
  ntasks = (ntasks > NUM_THREADS) ? ntasks : NUM_THREADS;

  auto a = (double*)calloc(N,sizeof(double));
  auto b = (double*)calloc(N,sizeof(double));
  auto c = (double*)calloc(N,sizeof(double));

  for (int i=0; i<N; i++) { a[i] = 0; b[i] = 1; c[i] = 2; }

  int ntrials = 10;
  double tsum = 0.0, tmin = 1e10, tmax = 0.0;

  auto print_timings = [&](auto& min, auto& sum, auto& max){
    cout << setw(7) << scientific << setprecision(3);
    cout << min*1e-6 << " ";
    cout << sum*1e-6 / ntrials << " ";
    cout << max*1e-6 << endl;
  };

  auto benchmark = [&](auto callable) {


    cout << setw(20) << std::right << "  ";
    cout << "Threads: [" << NUM_THREADS << "]  " << "Size [";
    cout << setw(2) << std::scientific << setprecision(1) << (float)N << "] ";
    cout << "NTasks: [" << ntasks << "]" << endl;
    cout << setw(20) << std::right << "  ";
    cout << setw(7) << std::internal << "Min.";
    cout << setw(10) << std::internal << "Avg.";
    cout << setw(10) << std::internal << "Max." << endl;


// -----------------------------------------------------------------------------------------------
    {
      threadpool11::Pool pool(NUM_THREADS);
      cout << setw(20) << std::right << "threadpool11(x2):  ";

      int nsize = ntasks;

      tsum = 0.0; tmax = 0.0; tmin = 1e10;
      for (int i=0; i<ntrials+1; i++)
      {
        // i==0 is a cold start for timing purposes
        Timer timer([&](int elapsed){ if (i>0) {
              tsum+=elapsed;
              if (elapsed < tmin) tmin = elapsed;
              if (elapsed > tmax) tmax = elapsed;}
          });

        std::future<void>* futures = new std::future<void>[nsize];
        auto begin = 0; auto end = N;
        int chunk = (end - begin) / nsize;
        for (int j = 0; j < nsize; ++j) {
          futures[j] = pool.postWork<void>([=]() {
              uint32_t threadstart = begin + j*chunk;
              uint32_t threadstop = (j == nsize - 1) ? end : threadstart + chunk;
              for (uint32_t it = threadstart; it < threadstop; ++it) {
                callable(it,a,b,c);
              }
            });
        }
        for (int j = 0; j < NUM_THREADS; ++j) {
          futures[j].get();
        }
      }
      print_timings(tmin,tsum,tmax);
    }


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
