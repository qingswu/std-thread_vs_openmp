#include <iostream>
#include <ThreadPool.hh>
#include <vector>
#include <cassert>
#include <chrono>
#include <algorithm>
#include <numeric>

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
    size_t nRows;
    size_t nCols;
};

void dummy(Matrix*);

void ompMatrixMultiply(Matrix& matA, Matrix& matB, Matrix& out, const size_t matSize)
{
    #pragma omp parallel for
    for (size_t i = 0; i < matSize; i++) {
        for (size_t j = 0; j < matSize; j++) {
            double sum = 0;
            for (size_t k = 0; k < matSize; k++) {
                sum += matA(i, k) * matB(k, j);
                dummy(&matA);
                dummy(&matB);
            }
            out(i, j) = sum;
        }
    }
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

int main(const int argc, const char** argv)
{
    size_t matSize = 1000;
    if (argc > 1) {
        matSize = atoi(argv[1]);
    }

    int numReps = 5;
    if (argc > 2) {
        numReps = atoi(argv[2]);
    }

    int innerLoopSize = 10;

    Matrix matA (matSize, matSize);
    Matrix matB (matSize, matSize);
    Matrix matC (matSize, matSize);

    using clock = std::chrono::high_resolution_clock;
    using duration = clock::duration;

    std::vector<duration> times;

    for (int rep = -1; rep < numReps; rep++) {
        matA.fill(1);
        matB.fill(2);
        matC.fill(0);

        auto begin = clock::now();
        for (int innerIter = 0; innerIter < innerLoopSize; innerIter++) {
            ompMatrixMultiply(matA, matB, matC, matSize);
        }
        auto end = clock::now();

        if (rep >= 0) {
            times.push_back((end - begin) / innerLoopSize);
        }
    }

    double minTime, meanTime, maxTime;
    std::tie(minTime, meanTime, maxTime) = processTimes(times);

    printf("Seconds: %10.4e  %10.4e  %10.4e\n", minTime, meanTime, maxTime);

    // const double numOps = 2 * matSize * matSize * matSize;

    // printf("FLOPS: %10.4e  %10.4e  %10.4e\n", numOps / minTime, numOps / (totTime / numReps), numOps / maxTime);

    return 0;
}
