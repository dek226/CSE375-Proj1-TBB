#include <vector>
//#include <oneapi/tbb/parallel_for.h>
//#include <oneapi/tbb/blocked_range.h>
#include "oneapi/tbb.h"
using namespace oneapi::tbb;

inline bool is_prime_naive(uint64_t n) {
    if (n <= 1){
            return false;
    }
        else {

            // Count the divisors of n
            int cnt = 0
            for (int i = 1; i <= n; i++) {
                if (n % i == 0){
                    return false;
                }
            }

        return true;
        }
    }


inline bool is_prime_TBB(uint64_t n, std::size_t grain_size) {
    if (n <= 1){
            return false;
    }
        else {

            // Count the divisors of n
            //int cnt = 0;
            //parellelizable part
            grain_size = 65536;
            std::atomic<bool> found_divisor(false);
            parallel_for(blocked_range<size_t>(1,n, grain_size), 
                [&found_divisor, n](const blocked_range<size_t>& r) {
                    for(size_t i=r.begin(); i<r.end(); i++){
                        if(found_divisor.load()){
                            return;
                        }
                        if (n % i == 0){
                            found_divisor.store(true);
                            return;
                        }
                    }
            });

            return !found_divisor.load();
        }
    }


inline bool is_prime_static(uint64_t n, int num_threads) {
    if (n <= 1){
            return false;
    }
        else {

            uint64_t start = 2;
            uint64_t end = n;
            // Count the divisors of n
            std::atomic<bool> found_divisor(false);
            uint64_t chunk = (end-start)/ num_threads;
            uint64_t remainder = (end-start) % num_threads;

            std::vector<std::thread> threads; //thread container preperation
            threads.reserve(num_threads); //prepares a container to efficiently manage a collection of std::thread objects by pre-allocating memory for a predetermined number of threads.

            for (int t=0; t<num_threads; t++){
                uint64_t curr_chunk = chunk;
                uint64_t rem_thread = 0;
                if (t<remainder){
                    curr_chunk++;
                    rem_thread = 1;
                }

                uint64_t loop_start = start + t*chunk + ((t+1)*rem_thread);
                uint64_t loop_end = loop_start + chunk_size;
                //create lambda
                threads.emplace_back([loop_start, loop_end, &found_divisor, n](){
                    for(uint64_t i= loop_start; i<loop_end; i++){
                        if(found_divisor.load()){
                            return;
                        }
                        if (n % i == 0){
                            found_divisor.store(true);
                            return;
                        }
                    }
                });

            }
        // join all threads
        for (auto &th : threads) th.join();

        return !found_divisor.load();
        }
    }


    inline bool is_prime_dynamic(uint64_t n, int num_threads, uint64_t batch_size) {
    if (n <= 1)
            return false;
        else {

            batch_size = 64;
            uint64_t start = 2;
            uint64_t end = n;
            
            std::atomic<uint64_t> next(start);
            // Count the divisors of n
            std::atomic<bool> found_divisor(false);

            std::vector<std::thread> threads; //thread container preperation
            threads.reserve(num_threads); //prepares a container to efficiently manage a collection of std::thread objects by pre-allocating memory for a predetermined number of threads.

            for (int t=0; t<num_threads; t++){
                threads.emplace_back([&next, &found_divisor, batch_size, end, n](){
                    while(true){
                        //std::atomic::fetch_add is a C++ function that performs an atomic read-modify-write operation on 
                        // an std::atomic variable. This means the entire operation—reading the current value, adding an argument to it, 
                        // and writing the new value back—is guaranteed to be indivisible and cannot be interrupted by other threads.
                        // Atomically claim a batch: get starting number for this worker
                        uint64_t loop_start = next.fetch_add(batch_size);
                        // If my_start is already >= end_exclusive, no work remains
                        if (loop_start >= end){
                            return;
                        }
                        // Compute my exclusive end (bounded by end_exclusive)
                        uint64_t loop_end = loop_start + batch_size;
                        if (loop_end > end) {
                            loop_end = end;
                        }
                        for(uint64_t i= loop_start; i<loop_end; i++){
                        if(found_divisor.load()){
                            return;
                        }
                        if (n % i == 0){
                            found_divisor.store(true);
                            return;
                        }
                        }

                    }
                });
            }
        // join all threads
        for (auto &th : threads) th.join();

        return !found_divisor.load();
        }
    }