#include <iostream>
#include <thread>
#include <atomic>
#include <vector>
#include <chrono>
#include <cmath> // For std::sqrt and std::pow
#include <numeric> // For std::accumulate
#include <fstream> // For CSV output
#include <iomanip> // For std::setprecision
#include "oneapi/tbb.h"
using namespace oneapi::tbb;
//g++ -O3 -std=c++17 prime.cpp -o prime_benchmark -pthread $(pkg-config --cflags --libs tbb)
//./prime_benchmark

bool is_prime_full_check(uint64_t n) {
    if (n <= 1) return false;
    for (uint64_t i = 2; i < n; i++) { 
        if (n % i == 0) {
            return false;
        }
    }
    return true;
}




bool is_prime_TBB(uint64_t n, std::size_t grain_size) {
    if (n <= 1){
            return false;
    }
        else {

            // Count the divisors of n
            //int cnt = 0;
            //parellelizable part
            //grain_size = 65536;
            std::atomic<bool> found_divisor(false);
            parallel_for(blocked_range<size_t>(2,n, grain_size), 
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


bool is_prime_static(uint64_t n, int num_threads) {
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
                uint64_t loop_end = loop_start + curr_chunk;
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


    bool is_prime_dynamic(uint64_t n, int num_threads, uint64_t batch_size) {
    if (n <= 1)
            return false;
        else {

            //batch_size = 64;
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



// --- TOP-LEVEL BENCHMARK WRAPPERS (SEQUENTIAL OUTER LOOP) ---


std::vector<uint64_t> find_primes_naive(uint64_t N_limit) {
    std::vector<uint64_t> primes;
    for (uint64_t i = 2; i <= N_limit; ++i) {
        if (is_prime_full_check(i)) {
            primes.push_back(i);
        }
    }
    return primes;
}


uint64_t find_primes_tbb_checker(uint64_t N_limit, std::size_t grain_size) {
    uint64_t prime_count = 0;
    for (uint64_t i = 2; i <= N_limit; ++i) {
        // TBB parallelization happens *inside* the checker function
        if (is_prime_TBB(i, grain_size)) {
            prime_count++;
        }
    }
    return prime_count;
}


uint64_t find_primes_static_checker(uint64_t N_limit, int num_threads) {
    uint64_t prime_count = 0;
    for (uint64_t i = 2; i <= N_limit; ++i) {
        // Static parallelization happens *inside* the checker function
        if (is_prime_static(i, num_threads)) {
            prime_count++;
        }
    }
    return prime_count;
}


uint64_t find_primes_dynamic_checker(uint64_t N_limit, int num_threads, uint64_t batch_size) {
    uint64_t prime_count = 0;
    for (uint64_t i = 2; i <= N_limit; ++i) {
        // Dynamic parallelization happens *inside* the checker function
        if (is_prime_dynamic(i, num_threads, batch_size)) {
            prime_count++;
        }
    }
    return prime_count;
}


// --- TIMING AND STATISTICS UTILITIES ---

template<typename F>
double time_it(F func) {
    auto start = std::chrono::high_resolution_clock::now();
    func();
    auto end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double>(end - start).count();
}

struct Stats {
    double mean = 0.0;
    double stdev = 0.0;
};

Stats calculate_stats(const std::vector<double>& times) {
    Stats s;
    if (times.empty()) return s;

    // Calculate Mean
    double sum = std::accumulate(times.begin(), times.end(), 0.0);
    s.mean = sum / times.size();

    // Calculate Standard Deviation (Population StDev)
    double sq_sum = 0.0;
    for (double t : times) {
        sq_sum += (t - s.mean) * (t - s.mean);
    }
    s.stdev = std::sqrt(sq_sum / times.size()); 
    
    return s;
}

// --- MAIN BENCHMARK HARNESS ---

int main() {
    // Define the N limits to test (find all primes up to this number)
    
    std::vector<uint64_t> limits = {
        1000,          // N=10^3
        1000000,         // N=10^6
        1500000          // N=10^6 * 1.5
    };

    // Define thread counts to test
    std::vector<int> thread_counts = {2, 4, 8};

    const int NUM_TRIALS = 5; // Number of trials per test
    const uint64_t TBB_GRAIN_SIZE = 65536; // Grain size for TBB
    const uint64_t DYNAMIC_BATCH_SIZE = 64; // Batch size for dynamic scheduling
    
    // Setup CSV file
    std::ofstream csv_file("prime_finder_performance_O_N_checkers.csv");
    csv_file << "N_Limit,Function,Threads,Avg_Time_s,StDev_s,Primes_Found\n";
    std::cout << "Starting Prime Finder Benchmark (O(N^2) complexity: sequential outer loop calling parallel O(N) primality checks)...\n";
    std::cout << "Results will be saved to prime_finder_performance_O_N_checkers.csv\n";
    
    // Set console precision
    std::cout << std::fixed << std::setprecision(6);

    for (uint64_t N_limit : limits) {
        std::cout << "\n========================================\n";
        std::cout << "Testing Range [2 to N] for N = " << N_limit << "\n";
        std::cout << "========================================\n";

        // --- 1. SEQUENTIAL RUN (Naive O(N) Checker) ---
        std::vector<double> naive_times;
        std::vector<uint64_t> result_primes; // Naive still collects the vector
        std::string func_name = "Naive_O_N";

        std::cout << "\n-- " << func_name << " (Sequential Outer/Inner, T=1) --\n";
        for (int t = 1; t <= NUM_TRIALS; ++t) {
            // Time the full 2 to N loop
            double time = time_it([&]() { result_primes = find_primes_naive(N_limit); });
            naive_times.push_back(time);
            std::cout << "Trial " << t << ": Time = " << time << "s\n";
        }
        Stats naive_stats = calculate_stats(naive_times);
        // Use result_primes.size() for the Naive version
        std::cout << " -> Primes Found: " << result_primes.size() << ".\n"; 
        std::cout << " -> AVG: " << naive_stats.mean << "s, StDev: " << naive_stats.stdev << "s.\n";
        
        // Write to CSV
        csv_file << N_limit << "," << func_name << "," << 1 << "," << naive_stats.mean << "," << naive_stats.stdev << "," << result_primes.size() << "\n";


        // --- 2. PARALLEL CHECKERS (TBB, Static, Dynamic) ---
        
        for (int num_threads : thread_counts) {
            std::cout << "\n-- Threads: " << num_threads << " --\n";
            
            // Note: Function signature updated to return uint64_t (the count)
            std::vector<std::pair<std::string, std::function<uint64_t()>>> parallel_checker_funcs = {
                {"TBB_Checker", [&]() { return find_primes_tbb_checker(N_limit, TBB_GRAIN_SIZE); }},
                {"Static_Checker", [&]() { return find_primes_static_checker(N_limit, num_threads); }},
                {"Dynamic_Checker", [&]() { return find_primes_dynamic_checker(N_limit, num_threads, DYNAMIC_BATCH_SIZE); }}
            };

            for (const auto& func_pair : parallel_checker_funcs) {
                const std::string& current_func_name = func_pair.first;
                const auto& func = func_pair.second;
                
                std::vector<double> current_times;
                uint64_t current_prime_count = 0; // Use simple count for non-naive functions
                
                std::cout << "  " << current_func_name << ":\n";
                for (int t = 1; t <= NUM_TRIALS; ++t) {
                    // Time the full 2 to N loop, storing only the count
                    double time = time_it([&]() { current_prime_count = func(); });
                    current_times.push_back(time);
                    
                    std::cout << "    Trial " << t << ": Time = " << time << "s, Threads = " << num_threads << ", N_Limit = " << N_limit << "\n";
                }
                
                Stats current_stats = calculate_stats(current_times);
                
                // Print average stats to console
                std::cout << "    -> Primes Found: " << current_prime_count << ".\n";
                std::cout << "    -> AVG: " << current_stats.mean << "s, StDev: " << current_stats.stdev << "s.\n";
                
                // Write to CSV
                csv_file << N_limit << "," << current_func_name << "," << num_threads << "," << current_stats.mean << "," << current_stats.stdev << "," << current_prime_count << "\n";
            }
        }
    }

    csv_file.close();
    std::cout << "\nBenchmark complete. Data written to prime_finder_performance_O_N_checkers.csv\n";
    return 0;
}