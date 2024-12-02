// Copyright (c) 2009-2024 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// TODO: Include the header files of classes that will be exported to Python.

#include "CapillarySoftShell.h"
#include <pybind11/pybind11.h>

namespace hoomd
    {
namespace md
    {

// TODO: Set the name of the python module to match ${COMPONENT_NAME} (set in
// CMakeLists.txt), prefixed with an underscore.
PYBIND11_MODULE(_capillary_core_shell_colloids, m)
    {
        detail::export_PotentialPair<CapillaryInteraction> (m,"CapillaryInteraction");
        detail::export_PotentialPair<SoftShell> (m,"SoftShell");
        // TODO: Call export_Class(m) for each C++ class to be exported to Python.

#ifdef ENABLE_HIP
        // TODO: Call export_ClassGPU(m) for each GPU enabled C++ class to be exported
        // to Python.
#endif
    }

    } // end namespace md
    } // end namespace hoomd
