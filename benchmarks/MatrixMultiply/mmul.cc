#include <ThreadPool.hh>
#include <LockfreeThreadPool.hh>
#include <ThreadPoolBase.hh>
#include <tp11helper.hh>
#include <omp.h>
#include <iostream>
#include <vector>
#include <cassert>
#include <chrono>
#include <algorithm>
#include <numeric>

// A matrix class, for convenience. This doesn't seem to slow things down appreciably
// compared to plain C-style arrays. (The speed is comparable to a similar C code.)
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

// This is defined in dummy.cc. It doesn't do anything, and exists only to fool
// the optimizer.
void dummy(Matrix*);

// This is the main function doing the matrix multiplication for the OpenMP case.
// I lifted this directly from HW 7. Blocking or parallelizing over other axes might
// be more efficient.
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

// This is the main function doing the ThreadPool-based multiplication. I pass the pool
// into the function to avoid re-spawning the threads each time.
void tpMatrixMultiply(ThreadPoolBase& tp, Matrix& matA, Matrix& matB, Matrix& out)
{
    // I don't believe capture-by-reference will cause dangling references here
    // since all work should be finished before this function returns.
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
    tp.ParallelFor(0, matA.nRows, 0, lambda);
}
// a threadpool11 version of the matrix multiply
void tp11MatrixMultiply(threadpool11::Pool& tp, Matrix& matA, Matrix& matB, Matrix& out)
{
    // I don't believe capture-by-reference will cause dangling references here
    // since all work should be finished before this function returns.
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
    ParallelFor(tp, 0, matA.nRows, 0, lambda);
}

// A convenience function to convert a std::duration into a double.
template<class T>
double convertDuration(const T& dur)
{
    auto count = std::chrono::duration_cast<std::chrono::nanoseconds>(dur).count();
    return count * 1.0e-9;
}

// This convenience function parses a vector of std::durations, finds the extrema
// and converts them to doubles.
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

// This is the core of the timing code. It runs the function `func` for the
// set of arguments `args`, and repeats it `numReps` times. The function `reset`
// is called using the same `args` list in between trials to re-initialize the
// matrices.
template<class FuncType, class ResetType, class ... ArgsType>
std::tuple<double, double, double>
timeMult(FuncType& func, ResetType& reset, int numReps, ArgsType&... args)
{
    using clock = std::chrono::high_resolution_clock;
    using duration = clock::duration;

    std::vector<duration> times;

    for (int rep = -1; rep < numReps; rep++) {
        reset(args...);

        auto begin = clock::now();
        func(args...);
        auto end = clock::now();

        // Throw out the first run (rep == -1) to avoid cold-start
        if (rep >= 0) {
            times.push_back(end - begin);
        }
    }

    return processTimes(times);
}

// Function to re-initialize the matrices
void clearMatrices(Matrix& matA, Matrix& matB, Matrix& matC)
{
    matA.fill(1);
    matB.fill(2);
    matC.fill(0);
}

// Unfortunately, since `timeMult` above calls `reset` with the same args as `func`,
// we need this overload for the ThreadPool case.
void tpClearMatrices(ThreadPoolBase& tp, Matrix& matA, Matrix& matB, Matrix& matC)
{
    clearMatrices(matA, matB, matC);
}
void tp11ClearMatrices(threadpool11::Pool& tp11, Matrix& matA, Matrix& matB, Matrix& matC)
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

    printf("              %10s  %10s  %10s\n", "Min", "Mean", "Max");
    printf("      OpenMP: %10.4e  %10.4e  %10.4e\n", minTime, meanTime, maxTime);

    // Use a block to limit scope of thread pool
    {
        ThreadPool tp (numThreads);

        std::tie(minTime, meanTime, maxTime) =
            timeMult(tpMatrixMultiply, tpClearMatrices, numReps, tp, matA, matB, matC);

        printf("  ThreadPool: %10.4e  %10.4e  %10.4e\n", minTime, meanTime, maxTime);
    }

    {
        LockfreeThreadPool lftp (numThreads);

        std::tie(minTime, meanTime, maxTime) =
            timeMult(tpMatrixMultiply, tpClearMatrices, numReps, lftp, matA, matB, matC);

        printf("LFThreadPool: %10.4e  %10.4e  %10.4e\n", minTime, meanTime, maxTime);
    }

    {
        threadpool11::Pool tp11(numThreads);
        std::tie(minTime, meanTime, maxTime) =
            timeMult(tp11MatrixMultiply, tp11ClearMatrices, numReps, tp11, matA, matB, matC);

        printf("threadpool11: %10.4e  %10.4e  %10.4e\n", minTime, meanTime, maxTime);
    }

    return 0;
}
