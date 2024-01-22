/*!
 * # responsibility_calcs.cpp
 *
 * Perform the key steps in generating responsibilities
 * (the E-step in the EM algorithm). Also calculates
 * statistics of all the sequences assigned to each cluster
 * for the "hard assignment statistics" routine.
 *
 * + getProbsCExt
 * Updates the responsibilities array using multithreading.
 *
 * + getProbsCExt_worker
 * A single thread of the updates launched by
 * getProbsCExt
 *
 * + getProbsCExt_terminal_masked
 * Updates the responsibilities array using multithreading,
 * but with n- and c-terminal masking supplied by caller.
 *
 * + getProbsCExt_terminal_masked_worker
 * A single thread of the updates launched by
 * getProbsCExt_terminal_masked
 *
 * + getProbsCExt_gapped
 * Updates the responsibilities array using multithreading,
 * but ignoring any element of x where x[i]==20 (corresponding
 * to gaps in amino acid sequences).
 *
 * + getProbsCExt_gapped_worker
 * A single thread of the updates launched by
 * getProbsCExt_masked
 */
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <thread>
#include "responsibility_calcs.h"


/*!
 * # getProbsCExt
 *
 * Calculates the updated responsibilities (the E-step
 * in the EM algorithm) for a batch of input data.
 * This function does not do any bounds checking, so it
 * is important for caller to do so. This function
 * is multithreaded and divides the work up into groups
 * of clusters each thread will handle.
 *
 * ## Args
 *
 * + `x` Pointer to the first element of the array x containing
 * input data. Should be an (N x C) array for N datapoints, sequence
 * length C. Each element indicates the item chosen at that position
 * in the raw data.
 * `mu` The current set of parameters of the model, in a (K x C x D)
 * array for K clusters, C sequence length, D options per sequence
 * element.
 * `resp` The (K x N) array of cluster responsibilities, for K clusters
 * and N datapoints.
 * `n_threads` Number of threads to launch.
 *
 * All operations are in place, nothing is returned.
 */
void getProbsCExt(py::array_t<uint8_t, py::array::c_style | py::array::forcecast> x, 
        py::array_t<double, py::array::c_style | py::array::forcecast> mu,
        py::array_t<double, py::array::c_style | py::array::forcecast> resp,
        int n_threads){
    py::ssize_t startRow, endRow;
    if (n_threads > mu.shape(0))
        n_threads = mu.shape(0);

    int chunkSize = (mu.shape(0) + n_threads - 1) / n_threads;
    std::vector<std::thread> threads(n_threads);


    for (int i=0; i < n_threads; i++){
        startRow = i * chunkSize;
        endRow = (i + 1) * chunkSize;
        if (endRow > mu.shape(0))
            endRow = mu.shape(0);
        threads[i] = std::thread(&getProbsCExt_worker,
                x, resp, mu, startRow, endRow);
    }

    for (auto& th : threads)
        th.join();
}


/*!
 * # getProbsCExt_worker
 *
 * Performs the E-step responsibility calculations for a subset
 * of the K available clusters.
 *
 * ## Args
 *
 * + `x` Pointer to the first element of the array x containing
 * input data. Should be an (N x C) array for N datapoints, sequence
 * length C. Each element indicates the item chosen at that position
 * in the raw data.
 * `resp` The (K x N) array of cluster responsibilities, for K clusters
 * and N datapoints.
 * `mu` The current set of parameters of the model, in a (K x C x D)
 * array for K clusters, C sequence length, D options per sequence
 * element.
 * `startRow` The first row of the resp and mu arrays to use for this
 * thread since this thread will only update for some clusters.
 * `endRow` The last row of the resp and mu arrays to use for this
 * thread.
 *
 * All operations are in place, nothing is returned.
 */
void *getProbsCExt_worker(py::array_t<uint8_t> x,
        py::array_t<double> resp, py::array_t<double> mu,
        py::ssize_t startRow, py::ssize_t endRow){
    auto resp_current = resp.mutable_unchecked<2>();
    auto x_current = x.unchecked<2>();
    auto mu_current = mu.unchecked<3>();
    int aa;
    double resp_value;

    for (py::ssize_t k=startRow; k < endRow; k++){
        for (py::ssize_t i=0; i < x.shape(0); i++){
            resp_value = 0;
            for (py::ssize_t j=0; j < x.shape(1); j++){
                aa = x_current(i,j);
                resp_value += mu_current(k,j,aa);
            }
            resp_current(k,i) = resp_value;
        }
    }
    return NULL;
}


/*!
 * # getProbsCExt_terminal_masked
 *
 * Calculates the updated responsibilities (the E-step
 * in the EM algorithm) for a batch of input data,
 * but with masking so that some n- and c-terminal columns
 * are excluded.
 * This function does not do any bounds checking, so it
 * is important for caller to do so. This function
 * is multithreaded and divides the work up into groups
 * of clusters each thread will handle.
 *
 * ## Args
 *
 * + `x` Pointer to the first element of the array x containing
 * input data. Should be an (N x C) array for N datapoints, sequence
 * length C. Each element indicates the item chosen at that position
 * in the raw data.
 * `mu` The current set of parameters of the model, in a (K x C x D)
 * array for K clusters, C sequence length, D options per sequence
 * element.
 * `resp` The (K x N) array of cluster responsibilities, for K clusters
 * and N datapoints.
 * `n_threads` Number of threads to launch.
 * `startCol` the first column (all previous are masked).
 * `endCol` the last column (all previous are masked).
 *
 * All operations are in place, nothing is returned.
 */
void getProbsCExt_terminal_masked(py::array_t<uint8_t, py::array::c_style | py::array::forcecast> x, 
        py::array_t<double, py::array::c_style | py::array::forcecast> mu,
        py::array_t<double, py::array::c_style | py::array::forcecast> resp,
        int n_threads, size_t startCol, size_t endCol){
    py::ssize_t startRow, endRow;
    if (n_threads > mu.shape(0))
        n_threads = mu.shape(0);

    int chunkSize = (mu.shape(0) + n_threads - 1) / n_threads;
    std::vector<std::thread> threads(n_threads);


    for (int i=0; i < n_threads; i++){
        startRow = i * chunkSize;
        endRow = (i + 1) * chunkSize;
        if (endRow > mu.shape(0))
            endRow = mu.shape(0);
        threads[i] = std::thread(&getProbsCExt_terminal_masked_worker,
                x, resp, mu, startRow, endRow,
                startCol, endCol);
    }

    for (auto& th : threads)
        th.join();
}


/*!
 * # getProbsCExt_terminal_masked_worker
 *
 * Performs the E-step responsibility calculations for a subset
 * of the K available clusters using a "mask" to exclude some
 * columns at the start and end of the sequence.
 *
 * ## Args
 *
 * + `x` Pointer to the first element of the array x containing
 * input data. Should be an (N x C) array for N datapoints, sequence
 * length C. Each element indicates the item chosen at that position
 * in the raw data.
 * `resp` The (K x N) array of cluster responsibilities, for K clusters
 * and N datapoints.
 * `mu` The current set of parameters of the model, in a (K x C x D)
 * array for K clusters, C sequence length, D options per sequence
 * element.
 * `startRow` The first row of the resp and mu arrays to use for this
 * thread since this thread will only update for some clusters.
 * `endRow` The last row of the resp and mu arrays to use for this
 * thread.
 * `startCol` the first column (all previous are masked).
 * `endCol` the last column (all previous are masked).
 *
 * All operations are in place, nothing is returned.
 */
void *getProbsCExt_terminal_masked_worker(py::array_t<uint8_t> x,
        py::array_t<double> resp, py::array_t<double> mu,
        py::ssize_t startRow, py::ssize_t endRow, py::ssize_t startCol,
        py::ssize_t endCol){
    auto resp_current = resp.mutable_unchecked<2>();
    auto x_current = x.unchecked<2>();
    auto mu_current = mu.unchecked<3>();
    int aa;
    double resp_value;

    for (py::ssize_t k=startRow; k < endRow; k++){
        for (py::ssize_t i=0; i < x.shape(0); i++){
            resp_value = 0;
            for (py::ssize_t j=startCol; j < endCol; j++){
                aa = x_current(i,j);
                resp_value += mu_current(k,j,aa);
            }
            resp_current(k,i) = resp_value;
        }
    }
    return NULL;
}



/*!
 * # getProbsCExt_gapped
 *
 * Calculates the updated responsibilities (the E-step
 * in the EM algorithm) for a batch of input data,
 * but ignoring gaps (cases where x[i]==20).
 * This function does not do any bounds checking, so it
 * is important for caller to do so. This function
 * is multithreaded and divides the work up into groups
 * of clusters each thread will handle.
 *
 * ## Args
 *
 * + `x` Pointer to the first element of the array x containing
 * input data. Should be an (N x C) array for N datapoints, sequence
 * length C. Each element indicates the item chosen at that position
 * in the raw data.
 * `mu` The current set of parameters of the model, in a (K x C x D)
 * array for K clusters, C sequence length, D options per sequence
 * element.
 * `resp` The (K x N) array of cluster responsibilities, for K clusters
 * and N datapoints.
 * `n_threads` Number of threads to launch.
 * `startCol` the first column (all previous are masked).
 * `endCol` the last column (all previous are masked).
 *
 * All operations are in place, nothing is returned.
 */
void getProbsCExt_gapped(py::array_t<uint8_t, py::array::c_style | py::array::forcecast> x, 
        py::array_t<double, py::array::c_style | py::array::forcecast> mu,
        py::array_t<double, py::array::c_style | py::array::forcecast> resp,
        int n_threads){
    py::ssize_t startRow, endRow;
    if (n_threads > x.shape(0))
        n_threads = x.shape(0);

    int chunkSize = (mu.shape(0) + n_threads - 1) / n_threads;
    std::vector<std::thread> threads(n_threads);


    for (int i=0; i < n_threads; i++){
        startRow = i * chunkSize;
        endRow = (i + 1) * chunkSize;
        if (endRow > mu.shape(0))
            endRow = mu.shape(0);
        threads[i] = std::thread(&getProbsCExt_gapped_worker,
                x, resp, mu, startRow, endRow);
    }

    for (auto& th : threads)
        th.join();
}


/*!
 * # getProbsCExt_gapped_worker
 *
 * Performs the E-step responsibility calculations for a subset
 * of the K available clusters excluding all gaps.
 *
 * ## Args
 *
 * + `x` Pointer to the first element of the array x containing
 * input data. Should be an (N x C) array for N datapoints, sequence
 * length C. Each element indicates the item chosen at that position
 * in the raw data.
 * `resp` The (K x N) array of cluster responsibilities, for K clusters
 * and N datapoints.
 * `mu` The current set of parameters of the model, in a (K x C x D)
 * array for K clusters, C sequence length, D options per sequence
 * element.
 * `startRow` The first row of the resp and mu arrays to use for this
 * thread since this thread will only update for some clusters.
 * `endRow` The last row of the resp and mu arrays to use for this
 * thread.
 *
 * All operations are in place, nothing is returned.
 */
void *getProbsCExt_gapped_worker(py::array_t<uint8_t> x, py::array_t<double> resp,
        py::array_t<double> mu, py::ssize_t startRow, py::ssize_t endRow){
    auto resp_current = resp.mutable_unchecked<2>();
    auto x_current = x.unchecked<2>();
    auto mu_current = mu.unchecked<3>();
    int aa;
    double resp_value;

    for (py::ssize_t k=startRow; k < endRow; k++){
        for (py::ssize_t i=0; i < x.shape(0); i++){
            resp_value = 0;
            for (py::ssize_t j=0; j < x.shape(1); j++){
                aa = x_current(i,j);
                if (aa == 20){
                    continue;
                }
                resp_value += mu_current(k,j,aa);
            }
            resp_current(k,i) = resp_value;
        }
    }
    return NULL;
}
