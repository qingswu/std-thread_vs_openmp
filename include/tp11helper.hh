#include <threadpool11/threadpool11.hpp>

// ParallelFor loop for threadpool11
template <typename T, typename... Params>
void ParallelFor(threadpool11::Pool& pool, uint32_t begin, uint32_t end, uint32_t n_tasks, T SerialFunction, Params&&... params) {

        n_tasks = (n_tasks > pool.getWorkerCount()) ? n_tasks : pool.getWorkerCount();
        std::future<void> futures[n_tasks];
        int chunk = (end - begin) / n_tasks;
        for (auto j = 0u; j < n_tasks; ++j) {
                futures[j] = pool.postWork<void>([=]() {
                                uint32_t threadstart = begin + j*chunk;
                                uint32_t threadstop = (j == n_tasks - 1) ? end : threadstart + chunk;
                                for (uint32_t it = threadstart; it < threadstop; ++it) {
                                        SerialFunction(it, params...);
                                }
                        });
        }
        // block until work is finished
        for (auto j = 0u; j < n_tasks; ++j) {
                futures[j].get();
        }
}
