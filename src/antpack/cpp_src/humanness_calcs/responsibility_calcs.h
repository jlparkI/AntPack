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
