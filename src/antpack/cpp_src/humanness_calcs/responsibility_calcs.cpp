/* Copyright (C) 2025 Jonathan Parkinson
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "responsibility_calcs.h"



namespace HumannessCalculations{

/// Calculates the updated responsibilities (the E-step
/// in the EM algorithm) for a batch of input data.
/// This function does not do any bounds checking, so it
/// is important for caller to do so. This function
/// is multithreaded and divides the work up into groups
/// of clusters each thread will handle.
void getProbsCExt(nb::ndarray<uint8_t, nb::shape<-1,-1>, nb::device::cpu, nb::c_contig> x, 
        nb::ndarray<double, nb::shape<-1,-1,-1>, nb::device::cpu, nb::c_contig> mu,
        nb::ndarray<double, nb::shape<-1,-1>, nb::device::cpu, nb::c_contig> resp,
        size_t n_threads) {

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

    for (int i=0; i < n_threads; i++) {
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


/// Performs the E-step responsibility calculations for a subset
/// of the K available clusters.
void *getProbsCExt_worker(uint8_t *x, double *resp,
            double *mu, int startRow, int endRow,
            int nClusters, int seqLen, int muDim2,
            int ndatapoints) {        

    double *resp_current, *mu_current, *mu_marker;
    uint8_t *x_current;
    int mu_row, mu_row_size = seqLen * muDim2;
    double resp_value;

    for (int k=startRow; k < endRow; k++) {
        resp_current = resp + k * ndatapoints;
        x_current = x;
        mu_row = k * mu_row_size;

        for (int i=0; i < ndatapoints; i++) {
            resp_value = 0;
            mu_marker = mu + mu_row;

            for (int j=0; j < seqLen; j++) {
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


/// Converts all gaps at the N- and C-terminal ends of each sequence
/// into value 21. This can then be passed to a masked scoring
/// function. The operation is performed in place. Once this
/// conversion has been performed, the result should not under
/// any circumstances be passed to a non-masked scoring function
/// (the non-masked scoring function does not use value 21).
void mask_terminal_deletions(nb::ndarray<uint8_t, nb::shape<-1,-1>, nb::device::cpu,
        nb::c_contig> x) {
    int nDatapoints = x.shape(0);
    int nCols = x.shape(1);
    uint8_t *xref = (uint8_t*) x.data();

    for (int i=0; i < nDatapoints; i++) {
        for (int k=0; k < nCols; k++) {
            if (xref[k] != 20)
                break;
            xref[k] = 21;
        }
        for (int k=nCols - 1; k > 0; k--) {
            if (xref[k] != 20)
                break;
            xref[k] = 21;
        }
        xref += nCols;
    }
}





/// Calculates the updated responsibilities (the E-step
/// in the EM algorithm) for a batch of input data,
/// but ignoring masked positions (cases where x[i] > 20).
/// This function does not do any bounds checking, so it
/// is important for caller to do so. This function
/// is multithreaded and divides the work up into groups
/// of clusters each thread will handle.
void getProbsCExt_masked(nb::ndarray<uint8_t, nb::shape<-1,-1>, nb::device::cpu,
                nb::c_contig> x, 
        nb::ndarray<double, nb::shape<-1,-1,-1>, nb::device::cpu, nb::c_contig> mu,
        nb::ndarray<double, nb::shape<-1,-1>, nb::device::cpu, nb::c_contig> resp,
        int n_threads) {
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


    for (int i=0; i < n_threads; i++) {
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


/// Performs the E-step responsibility calculations for a subset
/// of the K available clusters excluding masked positions.
void *getProbsCExt_masked_worker(uint8_t *x, double *resp,
        double *mu, int startRow, int endRow,
        int nClusters, int seqLen, int muDim2,
        int ndatapoints) {
    int i, j, k, mu_row;
    uint8_t *x_current;
    double *resp_current, *mu_current, *mu_marker;
    double resp_value;
    int mu_row_size = seqLen * muDim2;

    for (k=startRow; k < endRow; k++) {
        resp_current = resp + k * ndatapoints;
        x_current = x;
        mu_row = k * mu_row_size;

        for (i=0; i < ndatapoints; i++) {
            resp_value = 0;
            mu_marker = mu + mu_row;

            for (j=0; j < seqLen; j++) {
                if (*x_current > 20) {
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

}  // namespace HumannessCalculations
