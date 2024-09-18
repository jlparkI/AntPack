#ifndef RESPONSIBILITY_CALCS_H
#define RESPONSIBILITY_CALCS_H

#include <stdint.h>
#include <stdlib.h>
#include "nanobind/nanobind.h"
#include <nanobind/ndarray.h>


namespace nb = nanobind;


void getProbsCExt(nb::ndarray<uint8_t, nb::shape<-1,-1>, nb::device::cpu, nb::c_contig> x, 
        nb::ndarray<double, nb::shape<-1,-1,-1>, nb::device::cpu, nb::c_contig> mu,
        nb::ndarray<double, nb::shape<-1,-1>, nb::device::cpu, nb::c_contig> resp,
        size_t n_threads);
void *getProbsCExt_worker(uint8_t *x, double *resp,
            double *mu, int startRow, int endRow,
            int nClusters, int seqLen, int muDim2,
            int ndatapoints);    


void mask_terminal_deletions(nb::ndarray<uint8_t, nb::shape<-1,-1>, nb::device::cpu,
        nb::c_contig> x);


void getProbsCExt_masked(nb::ndarray<uint8_t, nb::shape<-1,-1>, nb::device::cpu,
                nb::c_contig> x, 
        nb::ndarray<double, nb::shape<-1,-1,-1>, nb::device::cpu, nb::c_contig> mu,
        nb::ndarray<double, nb::shape<-1,-1>, nb::device::cpu, nb::c_contig> resp,
        int n_threads);
void *getProbsCExt_masked_worker(uint8_t *x, double *resp,
        double *mu, int startRow, int endRow,
        int nClusters, int seqLen, int muDim2,
        int ndatapoints);

#endif
