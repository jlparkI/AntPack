// Note that this code is intended to set up the GUI as 
// a C++ extension and is not yet fully implemented since
// the build process will need some significant updtes
// to make this work. This is therefore in development
// at this time.

// C++ headers
#include <string>

// Library headers
#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/unordered_set.h>

// Project headers
#include "gui_cpp_src/application_runner.h"

namespace nb = nanobind;



NB_MODULE(antpack_gui_ext, m) {
    m.def("run_antpack_gui", &GUIApplicationRunner::run_antpack_gui);
}
