#include <iostream>
#include <ThreadPool.hh>
#include <vector>
#include <cassert>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <omp.h>

class Matrix
{
public:
    Matrix() = default;
    Matrix(size_t N, size_t M)
        : data(std::vector<double>(N * M, 0)), nRows(N), nCols(M) {}
    Matrix(size_t N, size_t M, double v)
        : data(std::vector<double>(N * M, v)), nRows(N), nCols(M) {}

    double& operator()(size_t i, size_t j)
    {
        assert(i < nRows);
        assert(j < nCols);
        return data[i * nRows + j];
    }

    void fill(double v)
    {
        for (double& elem : data) {
            elem = v;
        }
    }

private:
    std::vector<double> data;

public:
    size_t nRows;
    size_t nCols;
};

void dummy(Matrix*);

void ompMatrixMultiply(Matrix& matA, Matrix& matB, Matrix& out)
{
    #pragma omp parallel for
    for (size_t i = 0; i < matA.nRows; i++) {
        for (size_t j = 0; j < matA.nCols; j++) {
            double sum = 0;
            for (size_t k = 0; k < matA.nCols; k++) {
                sum += matA(i, k) * matB(k, j);
                dummy(&matA);
                dummy(&matB);
            }
            out(i, j) = sum;
        }
    }
}

void tpMatrixMultiply(ThreadPool& tp, Matrix& matA, Matrix& matB, Matrix& out)
{
    auto lambda = [&matA, &matB, &out] (size_t i) {
        for (size_t j = 0; j < matA.nCols; j++) {
            double sum = 0;
            for (size_t k = 0; k < matA.nCols; k++) {
                sum += matA(i, k) * matB(k, j);
                dummy(&matA);
                dummy(&matB);
            }
            out(i, j) = sum;
        }
    };
    tp.ParallelFor(0, matA.nRows, lambda);
}

template<class T>
double convertDuration(const T& dur)
{
    auto count = std::chrono::duration_cast<std::chrono::nanoseconds>(dur).count();
    return count * 1.0e-9;
}

template<class T>
auto processTimes(const std::vector<T>& times)
{
    auto extrema = std::minmax_element(times.begin(), times.end());
    T min = *extrema.first;
    T max = *extrema.second;
    T mean = std::accumulate(times.begin(), times.end(), T::zero()) / times.size();

    double minSec = convertDuration(min);
    double meanSec = convertDuration(mean);
    double maxSec = convertDuration(max);

    return std::make_tuple(minSec, meanSec, maxSec);
}

template<class FuncType, class ResetType, class ... ArgsType>
std::tuple<double, double, double>
timeMult(FuncType& func, ResetType& reset, int numReps, ArgsType&... args)
{
    using clock = std::chrono::high_resolution_clock;
    using duration = clock::duration;

    std::vector<duration> times;
    int innerLoopSize = 10;

    for (int rep = -1; rep < numReps; rep++) {
        reset(args...);

        auto begin = clock::now();
        for (int innerIter = 0; innerIter < innerLoopSize; innerIter++) {
            func(args...);
        }
        auto end = clock::now();

        if (rep >= 0) {
            times.push_back((end - begin) / innerLoopSize);
        }
    }

    return processTimes(times);
}

void clearMatrices(Matrix& matA, Matrix& matB, Matrix& matC)
{
    matA.fill(1);
    matB.fill(2);
    matC.fill(0);
}

void tpClearMatrices(ThreadPool& tp, Matrix& matA, Matrix& matB, Matrix& matC)
{
    clearMatrices(matA, matB, matC);
}

int main(const int argc, const char** argv)
{
    size_t matSize = 1000;
    if (argc > 1) {
        matSize = atoi(argv[1]);
    }

    int numThreads = 2;
    if (argc > 2) {
        numThreads = atoi(argv[2]);
    }
    omp_set_num_threads(numThreads);

    int numReps = 5;
    if (argc > 3) {
        numReps = atoi(argv[3]);
    }

    Matrix matA (matSize, matSize);
    Matrix matB (matSize, matSize);
    Matrix matC (matSize, matSize);

    double minTime, meanTime, maxTime;

    std::tie(minTime, meanTime, maxTime) =
        timeMult(ompMatrixMultiply, clearMatrices, numReps, matA, matB, matC);

    printf("            %10s  %10s  %10s\n", "Min", "Mean", "Max");
    printf("    OpenMP: %10.4e  %10.4e  %10.4e\n", minTime, meanTime, maxTime);

    // Use a block to limit scope of thread pool
    {
        ThreadPool tp (numThreads);

        std::tie(minTime, meanTime, maxTime) =
            timeMult(tpMatrixMultiply, tpClearMatrices, numReps, tp, matA, matB, matC);

        printf("ThreadPool: %10.4e  %10.4e  %10.4e\n", minTime, meanTime, maxTime);
    }

    return 0;
}
