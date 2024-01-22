#ifndef RESPONSIBILITY_CALCS_H
#define RESPONSIBILITY_CALCS_H

#include <stdint.h>
#include <stdlib.h>
#include <pybind11/numpy.h>


namespace py = pybind11;


void getProbsCExt(py::array_t<uint8_t, py::array::c_style | py::array::forcecast> x, 
        py::array_t<double, py::array::c_style | py::array::forcecast> mu,
        py::array_t<double, py::array::c_style | py::array::forcecast> resp,
        int n_threads);
void *getProbsCExt_worker(py::array_t<uint8_t> x,
        py::array_t<double> resp, py::array_t<double> mu,
        py::ssize_t startRow, py::ssize_t endRow);



void getProbsCExt_terminal_masked(py::array_t<uint8_t, py::array::c_style | py::array::forcecast> x, 
        py::array_t<double, py::array::c_style | py::array::forcecast> mu,
        py::array_t<double, py::array::c_style | py::array::forcecast> resp,
        int n_threads, size_t startCol, size_t endCol);
void *getProbsCExt_terminal_masked_worker(py::array_t<uint8_t> x,
        py::array_t<double> resp, py::array_t<double> mu,
        py::ssize_t startRow, py::ssize_t endRow, py::ssize_t startCol,
        py::ssize_t endCol);


void getProbsCExt_gapped(py::array_t<uint8_t, py::array::c_style | py::array::forcecast> x, 
        py::array_t<double, py::array::c_style | py::array::forcecast> mu,
        py::array_t<double, py::array::c_style | py::array::forcecast> resp,
        int n_threads);
void *getProbsCExt_gapped_worker(py::array_t<uint8_t> x, py::array_t<double> resp,
        py::array_t<double> mu, py::ssize_t startRow, py::ssize_t endRow);


#endif
