#ifndef RESPONSIBILITY_CALCS_H
#define RESPONSIBILITY_CALCS_H

#include <stdint.h>
#include <stdlib.h>
#include <pybind11/numpy.h>


namespace py = pybind11;


void getProbsCExt(py::array_t<uint8_t, py::array::c_style> x, 
        py::array_t<double, py::array::c_style> mu,
        py::array_t<double, py::array::c_style> resp,
        py::ssize_t n_threads);
void *getProbsCExt_worker(uint8_t *x, double *resp,
            double *mu, int startRow, int endRow,
            int nClusters, int seqLen, int muDim2,
            int ndatapoints);    


void mask_terminal_deletions(py::array_t<uint8_t, py::array::c_style> x);


void getProbsCExt_masked(py::array_t<uint8_t, py::array::c_style> x, 
        py::array_t<double, py::array::c_style> mu,
        py::array_t<double, py::array::c_style> resp,
        int n_threads);
void *getProbsCExt_masked_worker(uint8_t *x, double *resp,
        double *mu, int startRow, int endRow,
        int nClusters, int seqLen, int muDim2,
        int ndatapoints);

#endif
