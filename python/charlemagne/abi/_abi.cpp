// Yeah, I'm not going to lie: this is more of an incantation than prose.
// There's probably not many people apart from the folks who wrote Python
// and the C parts of NumPy who know how this works.
//
// I am proud that I can make and understand this, even with lots of refs.
//
// Python / Numpy modules
#pragma once
#include "sg03.hpp"
#include "tqm05.hpp"
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <cstring>
#include <exception>

// Time for some error catching :sob:
#define CHARLES_TRY try {
#define CHARLES_CATCH                                                        \
    } catch (const std::exception &e) {                                      \
        PyErr_SetString(PyExc_RuntimeError, e.what());                       \
        return nullptr;                                                      \
    } catch (...) {                                                          \
        PyErr_SetString(PyExc_RuntimeError, "unknown C++ exception");        \
        return nullptr;                                                      \
    }

/* Define docstrings */
static char module_docstring[] =
    "Python bindings for CharlemAGNe";
static char abi_solve_sg03_docstring[] =
    "Solve the SG03 model using the C++ module; do not persist the disk";
static char abi_solve_tqm05_docstring[] =
    "Solve the TQM05 model using the C++ module; do not persist the disk";

/* Forward declarations */
static PyObject *abi_solve_sg03(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *abi_solve_tqm05(PyObject *self, PyObject *args, PyObject *kwargs);

/* Define the methods */
static PyMethodDef charlemagne_methods[] = {
    {
    "abi_solve_sg03",
    (PyCFunction)(void(*)(void))abi_solve_sg03, METH_VARARGS | METH_KEYWORDS,
    abi_solve_sg03_docstring,
    },
    {
    "abi_solve_tqm05",
    (PyCFunction)(void(*)(void))abi_solve_tqm05, METH_VARARGS | METH_KEYWORDS,
    abi_solve_tqm05_docstring,
    },
    {NULL, NULL, 0, NULL}
};

// Define the Python module initializer
static struct PyModuleDef charlemagne_module = {
    PyModuleDef_HEAD_INIT, "charlemagne._abi",
    "CharlemAGNe C++ core (CPython limited API)",
    -1, charlemagne_methods,
    nullptr, nullptr, nullptr, nullptr,
};

PyMODINIT_FUNC PyInit__abi(void) {
    // Fail here when NumPy C-API is not available
    import_array();
    return PyModule_Create(&charlemagne_module);
}

// Profiles to dicts
static PyObject *profile_to_dict(const Profile &p) {
    // Create a new empty dictionary
    PyObject *dict = PyDict_New();
    // Check that allocating it succeeded
    if (!dict) return nullptr;
    // Get the number of elements in the profiles
    const npy_intp npts = (npy_intp)p.r_cyl.size();
    // assign it as the dimension of the array
    npy_intp dims[1] = {npts};

    // This makes a NumPy array from a vector of doubles
    auto add_doubles = [&](const char *name, const std::vector<double> &v) -> bool {
        // Make the array or die trying
        PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
        if (!arr) return false;
        // Copy data into array
        // memcpy should work for doubles
        std::memcpy(PyArray_DATA((PyArrayObject *)arr), v.data(), npts * sizeof(double));
        // Assign array to dictionary
        PyDict_SetItemString(dict, name, arr);
        // Remove the reference to the array
        Py_DECREF(arr);
        return true;
    };
    // This makes a NumPy array from a vector of booleans
    auto add_bools = [&](const char *name, const std::vector<bool> &v) -> bool {
        // Make the array or die trying
        PyObject *arr = PyArray_SimpleNew(1, dims, NPY_BOOL);
        if (!arr) return false;
        // Copy data into array
        npy_bool *d = (npy_bool *)PyArray_DATA((PyArrayObject *)arr);
        for (npy_intp i = 0; i < npts; i++) { d[i] = v[i] ? 1 : 0;}
        // Assign array to dictionary
        PyDict_SetItemString(dict, name, arr);
        // Remove the reference to the array
        Py_DECREF(arr);
        return true;
    };
    // This makes a NumPy array from a vector of ints
    auto add_ints = [&](const char *name, const std::vector<int> &v) -> bool {
        // Make the array or die trying
        PyObject *arr = PyArray_SimpleNew(1, dims, NPY_INT);
        if (!arr) return false;
        // Copy data into array
        npy_int *d = (npy_int *)PyArray_DATA((PyArrayObject *)arr);
        for (npy_intp i = 0; i < npts; i++) { d[i] = v[i];}
        // Assign array to dictionary
        PyDict_SetItemString(dict, name, arr);
        // Remove the reference to the array
        Py_DECREF(arr);
        return true;
    };
    // Apply to every item in profile
    if (
        !add_bools("converged", p.converged) ||
        !add_doubles("max_residual", p.max_residual) ||
        !add_ints("steps", p.steps) ||
        !add_bools("star_formation", p.star_formation) ||
        !add_doubles("r_cyl", p.r_cyl) ||
        !add_doubles("omega", p.omega) ||
        !add_doubles("Teff", p.Teff) ||
        !add_doubles("T", p.T) ||
        !add_doubles("cs", p.cs) ||
        !add_doubles("rho", p.rho) ||
        !add_doubles("kappa", p.kappa) ||
        !add_doubles("tauV", p.tauV) ||
        !add_doubles("Q", p.Q) ||
        !add_doubles("h", p.h) ||
        !add_doubles("eta", p.eta) ||
        !add_doubles("Mdot", p.Mdot) ||
        !add_doubles("Mdotprime", p.Mdotprime) ||
        !add_doubles("PPgas", p.PPgas) ||
        !add_doubles("PPrad", p.PPrad) ||
        !add_doubles("opacfac", p.opacfac) ||
        !add_doubles("dr_cyl", p.dr_cyl) ||
        !add_doubles("dA", p.dA) ||
        !add_doubles("dV", p.dV) ||
        !add_doubles("mass_disk", p.mass_disk) ||
        !add_doubles("Mdot_sf", p.Mdot_sf)
    ) {
        Py_DECREF(dict);
        return nullptr;
    }
    return dict;
}

// 1D SG03 disk
static PyObject *abi_solve_sg03(PyObject *, PyObject *args, PyObject *kwargs) {
    // args
    std::string opac_name;
    char *opac_name_ = nullptr;
    std::string opac_path;
    char *opac_path_ = nullptr;
    int opac_id;
    // kwargs
    double rmin = 6.0, rmax = 1.0e7;
    int npts = 1000;
    double smbh_mass = 1.0e8;
    double eddington_ratio = 0.5;
    double alpha = 0.01;
    double X_mass_fraction = 0.7;
    double radiative_efficiency = 0.1;
    bool viscosity_switch = false;
    // These are 
    static const char *kwlist[] = {
        "opac_name", "opac_path", "opac_id",
        "rmin", "rmax", "npts", 
        "smbh_mass", "eddington_ratio", "alpha", "X_mass_fraction",
        "radiative_efficiency",
        "viscosity_switch",
        nullptr};
    // We're looking for the arguments we want to use in this block,
    //   and assigning their values to the variables we defined.
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ssi|ddidddddp", (char **)kwlist,
        &opac_name_, &opac_path_, &opac_id,
        &rmin, &rmax, &npts,
        &smbh_mass, &eddington_ratio, &alpha, &X_mass_fraction,
        &radiative_efficiency, &viscosity_switch
      ))
      return nullptr;
    if (!opac_name_) return nullptr;
    if (!opac_path_) return nullptr;
    opac_name = (std::string)opac_name_;
    opac_path = (std::string)opac_path_;
    
    // Give it the old charles try
    CHARLES_TRY
    // Check for SG03
    SG03 disk(
        opac_name, opac_path, opac_id,
        smbh_mass, eddington_ratio, alpha, X_mass_fraction,
        radiative_efficiency, viscosity_switch
    );
    auto prof = disk.SolveDisk(rmin, rmax, npts, SolverOpts{});
    return profile_to_dict(prof);
    CHARLES_CATCH
}

// 1D TQM05 disk
static PyObject *abi_solve_tqm05(PyObject *, PyObject *args, PyObject *kwargs) {
    // args
    std::string opac_name;
    char *opac_name_ = nullptr;
    std::string opac_path;
    char *opac_path_ = nullptr;
    int opac_id;
    // kwargs
    double rmin = 6.0, rmax = 1.0e7;
    int npts = 1000;
    double smbh_mass = 1.e8;
    double velocity_dispersion = 0.0;
    double radiative_efficiency = 1.e-3;
    double transport_efficiency = 0.2;
    double supernova_pressure = 1.0;
    double X_mass_fraction = 0.7;
    double radius_inner = 6.0;
    double radius_outer = 0.0;
    double mdot_outer = 2.0169e+25;
    static const char *kwlist[] = {
        "opac_name", "opac_path", "opac_id",
        "rmin", "rmax", "npts", 
        "smbh_mass", "velocity_dispersion", "radiative_efficiency",
        "transport_efficiency", "supernova_pressure",
        "X_mass_fraction", "radius_inner", "radius_outer", "mdot_outer",
        nullptr};
    // We're looking for the arguments we want to use in this block,
    //   and assigning their values to the variables we defined.
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ssi|ddiddddddddd", (char **)kwlist,
        &opac_name_, &opac_path_, &opac_id,
        &rmin, &rmax, &npts,
        &smbh_mass,
        &velocity_dispersion,
        &radiative_efficiency,
        &transport_efficiency,
        &supernova_pressure,
        &X_mass_fraction,
        &radius_inner,
        &radius_outer,
        &mdot_outer
      ))
      return nullptr;
    if (!opac_name_) return nullptr;
    if (!opac_path_) return nullptr;
    opac_name = (std::string)opac_name_;
    opac_path = (std::string)opac_path_;
    
    // Give it the old charles try
    CHARLES_TRY
    // Check for TQM05
    TQM05 disk(
        opac_name, opac_path, opac_id,
        smbh_mass,
        velocity_dispersion,
        radiative_efficiency,
        transport_efficiency,
        supernova_pressure,
        X_mass_fraction,
        radius_inner,
        radius_outer,
        mdot_outer
    );
    auto prof = disk.SolveDisk(rmin, rmax, npts, SolverOpts{});
    return profile_to_dict(prof);
    CHARLES_CATCH
}
