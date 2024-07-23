/*
 * Python interface for LinearPartition
 *
 * Copyright 2024 Hyeshik Chang
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * “Software”), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
 * NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "Python.h"
#include "numpy/arrayobject.h"

#include <fstream>
int
trap_fprintf(FILE *fp, const char *fmt, ...)
{
    /* Block outputs to stderr */
    return 0;
}
int
trap_printf(const char *fmt, ...)
{
    /* Block outputs to stderr */
    return 0;
}

#define private protected
#undef fprintf
#define fprintf trap_fprintf
#undef printf
#define printf trap_printf

// Intercept the MEA structure output
#define __mea_hook__ (1) { threshknot_file_index = structure; } if (1)

// Monkey patch LinearPartition symbols to allow double-linking of E and V

/* EternaFold model */
#define main _linearpartition_e_main
#define BeamCKYParser LPE_BeamCKYParser
#define multi_base LPE_multi_base
#define multi_unpaired LPE_multi_unpaired
#define multi_paired LPE_multi_paired
#define external_unpaired LPE_external_unpaired
#define external_paired LPE_external_paired
#define base_pair LPE_base_pair
#define internal_1x1_nucleotides LPE_internal_1x1_nucleotides
#define helix_stacking LPE_helix_stacking
#define terminal_mismatch LPE_terminal_mismatch
#define bulge_0x1_nucleotides LPE_bulge_0x1_nucleotides
#define helix_closing LPE_helix_closing
#define dangle_left LPE_dangle_left
#define dangle_right LPE_dangle_right
#define internal_explicit LPE_internal_explicit
#define hairpin_length LPE_hairpin_length
#define bulge_length LPE_bulge_length
#define internal_length LPE_internal_length
#define internal_symmetric_length LPE_internal_symmetric_length
#define internal_asymmetry LPE_internal_asymmetry
#define hairpin_length_at_least LPE_hairpin_length_at_least
#define bulge_length_at_least LPE_bulge_length_at_least
#define internal_length_at_least LPE_internal_length_at_least
#define internal_symmetric_length_at_least LPE_internal_symmetric_length_at_least
#define internal_asymmetry_at_least LPE_internal_asymmetry_at_least
#include "contrib/feature_weight_e.h"
#include "LinearPartition.cpp"
#undef multi_base
#undef multi_unpaired
#undef multi_paired
#undef external_unpaired
#undef external_paired
#undef base_pair
#undef internal_1x1_nucleotides
#undef helix_stacking
#undef terminal_mismatch
#undef bulge_0x1_nucleotides
#undef helix_closing
#undef dangle_left
#undef dangle_right
#undef internal_explicit
#undef hairpin_length
#undef bulge_length
#undef internal_length
#undef internal_symmetric_length
#undef internal_asymmetry
#undef hairpin_length_at_least
#undef bulge_length_at_least
#undef internal_length_at_least
#undef internal_symmetric_length_at_least
#undef internal_asymmetry_at_least

#define lpv /* Vienna model */
#undef FASTCKY_BEAMCKYPAR_H
#undef FASTCKY_W
#undef main
#undef BeamCKYParser
#define main _linearpartition_v_main
#define BeamCKYParser LPV_BeamCKYParser
#define hash_pair LPV_hash_pair
#define comp LPV_comp
#define State LPV_State
#define Fast_Exp LPV_Fast_Exp
#define Fast_LogExpPlusOne LPV_Fast_LogExpPlusOne
#define Fast_LogPlusEquals LPV_Fast_LogPlusEquals
#define quickselect LPV_quickselect
#define quickselect_partition LPV_quickselect_partition
#define rtrim LPV_rtrim
#define pf_type LPV_pf_type
#define value_type LPV_value_type
#include "LinearPartition.cpp"

#undef State
#undef main
#undef private
#undef fprintf
#undef printf

struct basepair_prob {
    int32_t i;
    int32_t j;
    double prob;
};

static PyArray_Descr *partition_return_descr;

#define GET_BASEPAIR_PROB \
    PyObject * \
    get_basepair_prob(void) \
    { \
        PyArrayObject *res; \
        npy_intp dim; \
        struct basepair_prob *bpp; \
    \
        dim = Pij.size(); \
    \
        res = (PyArrayObject *)PyArray_SimpleNewFromDescr(1, &dim, partition_return_descr); \
        if (res == NULL) \
            return NULL; \
        Py_INCREF(partition_return_descr); \
    \
        assert(partition_return_descr->elsize == sizeof(struct basepair_prob)); \
        bpp = (struct basepair_prob *)PyArray_DATA(res); \
    \
        for (auto it = Pij.begin(); it != Pij.end(); ++it) { \
            bpp->i = it->first.first - 1; \
            bpp->j = it->first.second - 1; \
            bpp->prob = it->second; \
            bpp++; \
        } \
    \
        return (PyObject *)res; \
    }

class EternaBeamCKYParser : public LPE_BeamCKYParser {
public:
    using LPE_BeamCKYParser::LPE_BeamCKYParser;

    GET_BASEPAIR_PROB

    double
    get_free_energy(void)
    {
        State& viterbi=bestC[seq_length - 1];
        return viterbi.alpha;
    }
};

class ViennaBeamCKYParser : public LPV_BeamCKYParser {
public:
    using LPV_BeamCKYParser::LPV_BeamCKYParser;

    GET_BASEPAIR_PROB

    double
    get_free_energy(void)
    {
        LPV_State& viterbi=bestC[seq_length - 1];
        return -kT * viterbi.alpha / 100.0;
    }
};

PyDoc_STRVAR(linearpartition_partition_doc,
"partition(seq)\n\
\n\
Return the base-pairing probability matrix and ensembl free energy \
predicted by LinearPartition.");

static PyObject *
linearpartition_partition(PyObject *self, PyObject *args, PyObject *kwds)
{
    const char *seq, *mode="vienna";
    int beamsize=100, dangles=2;
    static const char *kwlist[] = {"seq", "mode", "beamsize", "dangles", NULL};
    enum { ETERNA, VIENNA } mode_enum;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|sii:partition",
                                     (char**)kwlist, &seq, &mode,
                                     &beamsize, &dangles))
        return NULL;

    if (strcmp(mode, "eterna") == 0)
        mode_enum = ETERNA;
    else if (strcmp(mode, "vienna") == 0)
        mode_enum = VIENNA;
    else {
        PyErr_SetString(PyExc_ValueError,
                        "mode must be either 'eterna' or 'vienna'");
        return NULL;
    }

    string rna_seq(seq);
    PyObject *probmtx;
    double free_energy;
    string mea_structure;

    /* Call LinearPartition */
    switch (mode_enum) {
    case ETERNA: {
        EternaBeamCKYParser parser(beamsize, true, false, "", "", false, 0.0,
            "", true, 3.0, "", false, false, 0.3, "", "", false, dangles);
        Py_BEGIN_ALLOW_THREADS
        parser.parse(rna_seq);
        Py_END_ALLOW_THREADS

        probmtx = parser.get_basepair_prob();
        if (probmtx == NULL)
            return NULL;
        free_energy = parser.get_free_energy();
        mea_structure = parser.threshknot_file_index;
        break;
    }

    case VIENNA: {
        ViennaBeamCKYParser parser(beamsize, true, false, "", "", false, 0.0,
            "", true, 3.0, "", false, false, 0.3, "", "", false, dangles);
        Py_BEGIN_ALLOW_THREADS
        parser.parse(rna_seq);
        Py_END_ALLOW_THREADS

        probmtx = parser.get_basepair_prob();
        if (probmtx == NULL)
            return NULL;
        free_energy = parser.get_free_energy();
        mea_structure = parser.threshknot_file_index;
        break;
    }

    default:
        PyErr_SetString(PyExc_RuntimeError, "unknown mode");
        return NULL;
    }

    PyObject *ret=NULL, *py_free_energy=NULL, *py_mea_structure=NULL;
    bool failed=true;

    do {
        ret = PyDict_New();
        if (ret == NULL)
            break;

        py_free_energy = PyFloat_FromDouble(free_energy);
        if (py_free_energy == NULL)
            break;

        py_mea_structure = PyUnicode_FromString(mea_structure.c_str());
        if (py_mea_structure == NULL)
            break;

        PyDict_SetItemString(ret, "structure", py_mea_structure);
        PyDict_SetItemString(ret, "free_energy", py_free_energy);
        PyDict_SetItemString(ret, "bpp", probmtx);

        failed = false;
    } while (0);

    Py_XDECREF(py_free_energy);
    Py_XDECREF(py_mea_structure);
    Py_DECREF(probmtx);

    if (failed) {
        Py_XDECREF(ret);
        return NULL;
    }
    else
        return ret;
}

static PyMethodDef linearpartition_methods[] = {
    {"partition",   (PyCFunction)linearpartition_partition,
     METH_VARARGS | METH_KEYWORDS, linearpartition_partition_doc},
    {NULL,          NULL} /* sentinel */
};

PyDoc_STRVAR(module_doc,
"CPython interface to LinearPartition");

static struct PyModuleDef linearpartitionmodule = {
    PyModuleDef_HEAD_INIT,
    "linearpartition",
    module_doc,
    0,
    linearpartition_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit_linearpartition(void)
{

    import_array();

    PyObject *op=Py_BuildValue("[(s, s), (s, s), (s, s)]",
                               "i", "i4", "j", "i4", "prob", "f8");
    if (op == NULL)
        return NULL;

    if (PyArray_DescrConverter(op, &partition_return_descr) == NPY_FAIL) {
        Py_DECREF(op);
        return NULL;
    }
    Py_DECREF(op);

    return PyModuleDef_Init(&linearpartitionmodule);
}
