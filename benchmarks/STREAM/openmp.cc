#include <iostream>
#include <iomanip>

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
    cout << setw(2) << std::scientific << setprecision(1) << (float)N << "]" <<endl;
    cout << setw(20) << std::right << "  ";
    cout << setw(7) << std::internal << "Min.";
    cout << setw(10) << std::internal << "Avg.";
    cout << setw(10) << std::internal << "Max." << endl;

// -----------------------------------------------------------------------------------------------
    cout << setw(20) << std::right << "OpenMP:  ";

    // omp cold start for timing purposes
    omp_set_dynamic(0);

    tsum = 0.0; tmax = 0.0; tmin = 1e10;
    for (int i=0; i<ntrials+1; i++)
    {
      Timer timer([&](int elapsed){ if (i>0) {
            tsum+=elapsed;
            if (elapsed < tmin) tmin = elapsed;
            if (elapsed > tmax) tmax = elapsed;
          } });

      #pragma omp parallel for num_threads(NUM_THREADS)
      for (int j=0; j<N; j++) {
        callable(j,a,b,c);
      }
    }
    print_timings(tmin,tsum,tmax);


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
