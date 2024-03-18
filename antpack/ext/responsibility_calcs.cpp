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
 * but ignoring any element of x where x[i] == 20 (corresponding
 * to gaps in amino acid sequences).
 *
 * + getProbsCExt_gapped_worker
 * A single thread of the updates launched by
 * getProbsCExt_masked
 *
 * + getProbsCExt_masked
 * Updates the responsibilities array using multithreading,
 * but ignoring any element of x where x[i] > 20 (corresponding
 * to masked sites in amino acid sequences).
 *
 * + getProbsCExt_masked_worker
 * A single thread of the updates launched by
 * getProbsCExt_masked
 */
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <thread>
#include <iostream>
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
void getProbsCExt(py::array_t<uint8_t, py::array::c_style> x, 
        py::array_t<double, py::array::c_style> mu,
        py::array_t<double, py::array::c_style> resp,
        py::ssize_t n_threads){

    if (n_threads > mu.shape(0))
        n_threads = mu.shape(0);
    if (n_threads > x.shape(0))
        n_threads = x.shape(0);

    int nClusters = mu.shape(0), seqLen = mu.shape(1), muDim2 = mu.shape(2);
    int ndatapoints = x.shape(0);
    int startRow, endRow;
    uint8_t *xref = (uint8_t*) x.data();
    double *muref = (double*) mu.data(), *respref = (double*) resp.data();

    int chunkSize = (mu.shape(0) + n_threads - 1) / n_threads;
    std::vector<std::thread> threads(n_threads);

    for (int i=0; i < n_threads; i++){
        startRow = i * chunkSize;
        endRow = (i + 1) * chunkSize;
        if (endRow > mu.shape(0))
            endRow = mu.shape(0);

        threads[i] = std::thread(&getProbsCExt_worker,
                xref, respref, muref,
                startRow, endRow, nClusters, seqLen,
                muDim2, ndatapoints);
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
 * `nClusters` dim0 of mu.
 * `seqLen` dim1 of mu.
 * `muDim2` dim2 of mu.
 * `ndatapoints` number of datapoints, dim0 of x.
 *
 * All operations are in place, nothing is returned.
 */
void *getProbsCExt_worker(uint8_t *x, double *resp,
            double *mu, int startRow, int endRow,
            int nClusters, int seqLen, int muDim2,
            int ndatapoints){        

    double *resp_current, *mu_current, *mu_marker;
    uint8_t *x_current;
    int mu_row, mu_row_size = seqLen * muDim2;
    double resp_value;

    for (int k=startRow; k < endRow; k++){
        resp_current = resp + k * ndatapoints;
        x_current = x;
        mu_row = k * mu_row_size;

        for (int i=0; i < ndatapoints; i++){
            resp_value = 0;
            mu_marker = mu + mu_row;

            for (int j=0; j < seqLen; j++){
                mu_current = mu_marker + *x_current;
                resp_value += *mu_current;
                x_current++;
                mu_marker += muDim2;
            }
            *resp_current = resp_value;
            resp_current++;
        }
    }
    return NULL;
}


/*!
 * # mask_terminal_deletions
 *
 * Converts all gaps at the N- and C-terminal ends of each sequence
 * into value 21. This can then be passed to a masked scoring
 * function. The operation is performed in place. Once this
 * conversion has been performed, the result should not under
 * any circumstances be passed to a non-masked scoring function
 * (the non-masked scoring function does not use value 21).
 *
 * ## Args
 *
 * + `x` Pointer to the first element of the array x containing
 * input data. Should be an (N x C) array for N datapoints, sequence
 * length C. Each element indicates the item chosen at that position
 * in the raw data.
 *
 * All operations are in place, nothing is returned.
 */
void mask_terminal_deletions(py::array_t<uint8_t, py::array::c_style> x){
    int nDatapoints = x.shape(0);
    int nCols = x.shape(1);
    uint8_t *xref = (uint8_t*) x.data();

    for (int i=0; i < nDatapoints; i++){
        for (int k=0; k < nCols; k++){
            if (xref[k] != 20)
                break;
            xref[k] = 21;
        }
        for (int k=nCols - 1; k > 0; k--){
            if (xref[k] != 20)
                break;
            xref[k] = 21;
        }
        xref += nCols;
    }
}





/*!
 * # getProbsCExt_masked
 *
 * Calculates the updated responsibilities (the E-step
 * in the EM algorithm) for a batch of input data,
 * but ignoring masked positions (cases where x[i] > 20).
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
void getProbsCExt_masked(py::array_t<uint8_t, py::array::c_style> x, 
        py::array_t<double, py::array::c_style> mu,
        py::array_t<double, py::array::c_style> resp,
        int n_threads){
    if (n_threads > x.shape(0))
        n_threads = x.shape(0);
    if (n_threads > x.shape(0))
        n_threads = x.shape(0);

    int nClusters = mu.shape(0), seqLen = mu.shape(1), muDim2 = mu.shape(2);
    int ndatapoints = x.shape(0);
    int startRow, endRow;
    uint8_t *xref = (uint8_t*) x.data();
    double *muref = (double*) mu.data(), *respref = (double*) resp.data();

    int chunkSize = (mu.shape(0) + n_threads - 1) / n_threads;
    std::vector<std::thread> threads(n_threads);


    for (int i=0; i < n_threads; i++){
        startRow = i * chunkSize;
        endRow = (i + 1) * chunkSize;
        if (endRow > mu.shape(0))
            endRow = mu.shape(0);
        threads[i] = std::thread(&getProbsCExt_masked_worker,
                xref, respref, muref, startRow, endRow,
                nClusters, seqLen, muDim2, ndatapoints);
    }

    for (auto& th : threads)
        th.join();
}


/*!
 * # getProbsCExt_masked_worker
 *
 * Performs the E-step responsibility calculations for a subset
 * of the K available clusters excluding masked positions.
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
 * `nClusters` dim0 of mu.
 * `seqLen` dim1 of mu.
 * `muDim2` dim2 of mu.
 * `ndatapoints` dim0 of x.
 *
 * All operations are in place, nothing is returned.
 */
void *getProbsCExt_masked_worker(uint8_t *x, double *resp,
        double *mu, int startRow, int endRow,
        int nClusters, int seqLen, int muDim2,
        int ndatapoints){
    int i, j, k, mu_row;
    uint8_t *x_current;
    double *resp_current, *mu_current, *mu_marker;
    double resp_value;
    int mu_row_size = seqLen * muDim2;

    for (k=startRow; k < endRow; k++){
        resp_current = resp + k * ndatapoints;
        x_current = x;
        mu_row = k * mu_row_size;

        for (i=0; i < ndatapoints; i++){
            resp_value = 0;
            mu_marker = mu + mu_row;

            for (j=0; j < seqLen; j++){
                if (*x_current > 20){
                    x_current++;
                    mu_marker += muDim2;
                    continue;
                }
                mu_current = mu_marker + *x_current;
                resp_value += *mu_current;
                x_current++;
                mu_marker += muDim2;
            }
            *resp_current = resp_value;
            resp_current++;
        }
    }
    return NULL;
}
