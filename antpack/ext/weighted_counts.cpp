/*!
 * # weighted_counts.cpp
 *
 * Perform the key steps involved in generating weighted counts
 * for the M-step in the EM algorithm.
 *
 * + getWeightedCountCExt_main
 * Updates the weighted counts (used to generate the updated mu
 * values) for an input data array, using multithreading.
 *
 * + getWeightedCountCExt_worker
 * A single thread of the wcount updates launched by
 * getWeightedCountCExt_main
 *
 * + getWeightedCountCExt_single_thread
 * A single thread version of getWeightedCountCExt_main.
 */
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <thread>
#include "weighted_counts.h"

#define MEMORY_ERROR 1
#define NO_ERROR 0
#define THREAD_ERROR 2




/*!
 * # getWeightedCountCExt_main
 *
 * Updates the wcount array containing the responsibility-
 * weighted counts. In the M step of EM optimization, this
 * array will become the new mu values. This function does
 * no bounds checking, so it is very important that
 * caller check bounds.
 *
 * ## Args
 *
 * + `x` Pointer to the first element of the array x containing
 * input data. Should be an (N x C) array for N datapoints, sequence
 * length C. Each element indicates the item chosen at that position
 * in the raw data.
 * `wcount` Pointer to first element of weighted count array; will
 * be updated in place. Array is of shape (K x C x D) for K clusters,
 * C sequence length, D options per sequence element.
 * `resp` The (K x N) array of cluster responsibilities, for K clusters
 * and N datapoints.
 * `wcount_dim0` shape[0] of wcount
 * `wcount_dim1` shape[1] of wcount
 * `wcount_dim2` shape[2] of wcount
 * `x_dim0` shape 0 of x
 * `x_dim1` shape 1 of x
 * `n_threads` Number of threads to launch.
 *
 * All operations are in place, nothing is returned.
 */
int getWeightedCountCExt_main(uint8_t *x, double *wcount, double *resp,
                   int wcount_dim0, int wcount_dim1, int wcount_dim2,
                   int x_dim0, int x_dim1, int n_threads){

    int startRow, endRow;
    int chunkSize = (wcount_dim0 + n_threads - 1) / n_threads;
    std::vector<std::thread> threads(n_threads);

    if (n_threads > x_dim0)
        n_threads = x_dim0;

    for (int i=0; i < n_threads; i++){
        startRow = i * chunkSize;
        endRow = (i + 1) * chunkSize;
        if (endRow > wcount_dim0)
            endRow = wcount_dim0;
        threads[i] = std::thread(&getWeightedCountCExt_worker,
                x, resp, wcount, wcount_dim1, wcount_dim2,
                x_dim0, x_dim1, startRow, endRow);
    }

    for (auto& th : threads)
        th.join();
    return NO_ERROR;
}



/*!
 * # getWeightedCountCExt_single_thread
 *
 * Updates the wcount array containing the responsibility-
 * weighted counts. It only performs this operation on
 * clusters in (startRow, endRow); each thread is assigned
 * to complete some subset of the total clusters.
 *
 * ## Args
 *
 * + `x` Pointer to the first element of the array x containing
 * input data. Should be an (N x C) array for N datapoints, sequence
 * length C. Each element indicates the item chosen at that position
 * in the raw data.
 * `resp` The (K x N) array of cluster responsibilities, for K clusters
 * and N datapoints.
 * `wcount` Pointer to first element of weighted count array; will
 * be updated in place.
 * `wcount_dim1` shape[1] of wcount
 * `wcount_dim2` shape[2] of wcount
 * `x_dim0` shape 0 of x
 * `x_dim1` shape 1 of x
 * `startRow` The first row of wcount to update for this thread.
 * `endRow` The last row of wcount to update for this thread.
 *
 * All operations are in place, nothing is returned.
 */
void *getWeightedCountCExt_worker(uint8_t *x, double *resp,
        double *wcount, int wcount_dim1, int wcount_dim2,
        int x_dim0, int x_dim1, int startRow, int endRow){

    int i, j, k, wcount_row;
    uint8_t *x_current;
    double *resp_current, *wcount_current, *wcount_marker;
    double resp_value;
    int wcount_row_size = wcount_dim1 * wcount_dim2;


    for (k=startRow; k < endRow; k++){
        x_current = x;
        resp_current = resp + k * x_dim0;
        wcount_row = k * wcount_row_size;
        for (i=0; i < x_dim0; i++){
            resp_value = *resp_current;
            wcount_marker = wcount + wcount_row;
            for (j=0; j < x_dim1; j++){
                wcount_current = wcount_marker + *x_current;
                *wcount_current += resp_value;
                x_current++;
                wcount_marker += wcount_dim2;
            }
            resp_current++;
        }
    }
    return NULL;
}





/*!
 * # getWeightedCountCExt_single_thread
 *
 * Updates the wcount array containing the responsibility-
 * weighted counts. In the M step of EM optimization, this
 * array will become the new mu values. This function does
 * no bounds checking, so it is very important that
 * caller check bounds. This is a single-threaded version
 * (calling the multi-threaded version if only one thread
 * is desired would incur a little extra unnecessary
 * overhead).
 *
 * ## Args
 *
 * + `x` Pointer to the first element of the array x containing
 * input data. Should be an (N x C) array for N datapoints, sequence
 * length C. Each element indicates the item chosen at that position
 * in the raw data.
 * `wcount` Pointer to first element of weighted count array; will
 * be updated in place.
 * `resp` The (K x N) array of cluster responsibilities, for K clusters
 * and N datapoints.
 * `wcount_dim0` shape[0] of wcount
 * `wcount_dim1` shape[1] of wcount
 * `wcount_dim2` shape[2] of wcount
 * `x_dim0` shape 0 of x
 * `x_dim1` shape 1 of x
 *
 * All operations are in place, nothing is returned.
 */
int getWeightedCountCExt_single_thread(uint8_t *x, double *wcount, double *resp,
                   int wc_dim0, int wc_dim1, int wc_dim2,
                   int x_dim0, int x_dim1){
    int i, j, k, wcount_row;
    uint8_t *x_current;
    double *resp_current, *wcount_current, *wcount_marker;
    double resp_value;
    int wcount_row_size = wc_dim1 * wc_dim2;

    for (k=0; k < wc_dim0; k++){
        x_current = x;
        resp_current = resp + k * x_dim0;
        wcount_row = k * wcount_row_size;
        for (i=0; i < x_dim0; i++){
            resp_value = *resp_current;
            wcount_marker = wcount + wcount_row;
            for (j=0; j < x_dim1; j++){
                wcount_current = wcount_marker + *x_current;
                *wcount_current += resp_value;
                x_current++;
                wcount_marker += wc_dim2;
            }
            resp_current++;
        }
    }
    return NO_ERROR;
}
