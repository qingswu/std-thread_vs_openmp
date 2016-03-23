
#include "ThreadPool.hh"
#include "Timer.hh"
#include "omp.h"

#include <iostream>


using namespace std;


void scale(int i, double* a, double* b) {
        a[i]=4*b[i];
}

int main () {

        ThreadPool pool(8);

        int N = 1e9;
        auto a = (double*)calloc(N,sizeof(double));
        auto b = (double*)calloc(N,sizeof(double));
        for (int i=0; i<N; i++) { b[i] = i; }

        int ntrials = 10;
        double tperformance = 0.0;

// -----------------------------------------------------------------------------------------------
        cout << "single thread performance\n";

        for (int j=0; j<N; j++) {
                a[j]=4*b[j];
        }

        tperformance = 0.0;
        for (int i=0; i<ntrials; i++)
        {
                Timer timer([&](int elapsed){
                                cout << "Trial " << i << ": "<< elapsed*1e-6 << " ms\n";
                                tperformance+=elapsed;
                        });
                for (int j=0; j<N; j++) {
                        a[j]=4*b[j];
                }
        }
        cout << "Average: " << tperformance*1e-6 / ntrials << " ms\n\n";

// -----------------------------------------------------------------------------------------------
        cout << "std::thread ThreadPool performance\n";

        // cold start for timing purposes
        pool.ParallelFor(0,N,scale,a,b);

        tperformance = 0.0;
        for (int i=0; i<ntrials; i++)
        {
                Timer timer([&](int elapsed){
                                cout << "Trial " << i << ": "<< elapsed*1e-6 << " ms\n";
                                tperformance+=elapsed;
                        });
                pool.ParallelFor(0,N,scale,a,b);
        }
        cout << "Average: " << tperformance*1e-6 / ntrials << " ms\n\n";

// -----------------------------------------------------------------------------------------------
        cout << "OpenMP performance\n";

        // omp cold start for timing purposes
        omp_set_dynamic(0);
        #pragma omp parallel for num_threads(8)
        for (int j=0; j<N; j++) {
                a[j]=4*b[j];
        }


        tperformance = 0.0;
        for (int i=0; i<ntrials; i++)
        {
                Timer timer([&](int elapsed){
                                cout << "Trial " << i << ": "<< elapsed*1e-6 << " ms\n";
                                tperformance+=elapsed;
                        });

                #pragma omp parallel for num_threads(8)
                for (int j=0; j<N; j++) {
                        a[j]=4*b[j];
                }
        }
        cout << "Average: " << tperformance*1e-6 / ntrials << " ms\n\n";

// -----------------------------------------------------------------------------------------------

        return 0;
}
