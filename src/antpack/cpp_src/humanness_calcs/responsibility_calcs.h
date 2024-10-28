#ifndef RESPONSIBILITY_CALCS_H
#define RESPONSIBILITY_CALCS_H

// C++ headers
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <thread>
#include <iostream>

// Library headers
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

// Project headers


namespace nb = nanobind;



namespace HumannessCalculations{


/// @brief Calculates the cluster responsibilities and sequence probabilities.
void getProbsCExt(nb::ndarray<uint8_t, nb::shape<-1,-1>, nb::device::cpu, nb::c_contig> x, 
        nb::ndarray<double, nb::shape<-1,-1,-1>, nb::device::cpu, nb::c_contig> mu,
        nb::ndarray<double, nb::shape<-1,-1>, nb::device::cpu, nb::c_contig> resp,
        size_t n_threads);

/// @brief Executes a single thread of calculations for getProbsCExt.
void *getProbsCExt_worker(uint8_t *x, double *resp,
            double *mu, int startRow, int endRow,
            int nClusters, int seqLen, int muDim2,
            int ndatapoints);    

/// @brief Masks gaps on the n- and c-terminals of an input array by converting them
/// @brief to a "mask" value.
void mask_terminal_deletions(nb::ndarray<uint8_t, nb::shape<-1,-1>, nb::device::cpu,
        nb::c_contig> x);

/// @brief Calculates responsibilities and sequence probabilities for an encoded
///        input sequence.
void getProbsCExt_masked(nb::ndarray<uint8_t, nb::shape<-1,-1>, nb::device::cpu,
                nb::c_contig> x, 
        nb::ndarray<double, nb::shape<-1,-1,-1>, nb::device::cpu, nb::c_contig> mu,
        nb::ndarray<double, nb::shape<-1,-1>, nb::device::cpu, nb::c_contig> resp,
        int n_threads);

/// @brief Executes one thread of calculations for getProbsCExt_masked.
void *getProbsCExt_masked_worker(uint8_t *x, double *resp,
        double *mu, int startRow, int endRow,
        int nClusters, int seqLen, int muDim2,
        int ndatapoints);


}  // namespace HumannessCalculations

#endif
