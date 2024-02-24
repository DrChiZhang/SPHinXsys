
/* -------------------------------------------------------------------------*
*								pySPHinXsys									*
* --------------------------------------------------------------------------*
* pySPHinXsys provides a python binding for SPHinXsys library.				*
*																			*
* Copyright(c) 2017-2024 The authors and the authors' affiliations.		    *
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
* @file 	sph_system_module.cpp
* @brief 	Python binding for SPH_System in SPHinxsys library.
* @author	Chi ZHANG
* @version  1.0 
*			Start the python binding work.
*			-- Chi ZHANG
*/
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "dambreak_module.hpp"

namespace py = pybind11;

PYBIND11_MODULE(pysphinxsys, m)
{
    py::class_<Environment>(m, "Dambreak")
        .def(py::init<const int &>())
        .def("CmakeTest", &Environment::cmakeTest)
        .def("RunCase", &Environment::runCase);
}