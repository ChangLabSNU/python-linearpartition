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

#define main _linearpartition_main
#define lpv /* Vienna model */
#define private protected
#undef fprintf
#define fprintf trap_fprintf
#include "LinearPartition.cpp"
#undef main
#undef private
#undef fprintf

class MyBeamCKYParser : public BeamCKYParser {
public:
    MyBeamCKYParser(int beamsize, bool no_sharpturn, bool verbose,
                    const string &bpp_file, const string &bpp_file_index,
                    bool pf_only, float bpp_cutoff, const string &forest_file,
                    bool mea, float MEA_gamma, const string &MEA_file_index,
                    bool MEA_bpseq, bool ThreshKnot, float ThreshKnot_threshold,
                    const string &ThreshKnot_file_index, const string &shape_file_path,
                    bool fasta, int dangles)
        : BeamCKYParser(beamsize, no_sharpturn, verbose, bpp_file, bpp_file_index,
                        pf_only, bpp_cutoff, forest_file, mea, MEA_gamma, MEA_file_index,
                        MEA_bpseq, ThreshKnot, ThreshKnot_threshold, ThreshKnot_file_index,
                        shape_file_path, fasta, dangles)
    {
    }

    PyObject *
    get_basepair_prob(void)
    {
        PyArrayObject *res;
        const npy_intp dims[2]={seq_length, seq_length};
        double *buf;

        res = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
        if (res == NULL)
            return NULL;

        PyArray_FILLWBYTE(res, 0);

        buf = (double*)PyArray_DATA(res);

        // Taken from bpp.cpp of LinearPartion
        int turn = no_sharp_turn ? 3 : 0;
        for (unsigned int i = 1; i <= seq_length; i++) {
            for (unsigned int j = i + turn + 1; j <= seq_length; j++) {
                pair<int, int> key = make_pair(i,j);
                auto got = Pij.find(key);

                if (got != Pij.end()) {
                    buf[(i - 1) * seq_length + j - 1] = got->second;
                    buf[(j - 1) * seq_length + i - 1] = got->second;
                }
            }
        }

        return (PyObject *)res;
    }

    double
    get_free_energy(void)
    {
        State& viterbi=bestC[seq_length - 1];
        return -kT * viterbi.alpha / 100.0;
    }
};

PyDoc_STRVAR(linearpartition_partition_doc,
"partition(seq)\n\
\n\
Return the base-pairing probability matrix and ensembl free energy \
predicted by LinearPartition.");

static PyObject *
linearpartition_partition(PyObject *self, PyObject *args)
{
    const char *seq;
    Py_ssize_t len;

    if (!PyArg_ParseTuple(args, "s#:partition", &seq, &len))
        return NULL;

    /* LinearPartition arguments */
    int beamsize = 100;
    bool sharpturn = false;
    bool pf_only = false;
    float bpp_cutoff = 0.0;

    float MEA_gamma = 3.0;
    float ThreshKnot_threshold = 0.3;
    int dangles = 2;
    // End of LinearPartion parameters

    string rna_seq(seq);

    /* Call LinearPartition */
    MyBeamCKYParser parser(beamsize, !sharpturn, false, "", "",
                           pf_only, bpp_cutoff, "", false, MEA_gamma, "",
                           false, false, ThreshKnot_threshold, "",
                           "", false, dangles);
    Py_BEGIN_ALLOW_THREADS
    parser.parse(rna_seq);
    Py_END_ALLOW_THREADS

    PyObject *ret, *probmtx;

    probmtx = parser.get_basepair_prob();
    if (probmtx == NULL)
        return NULL;

    ret = Py_BuildValue("Od", probmtx, parser.get_free_energy());
    Py_DECREF(probmtx);

    return ret;
}

static PyMethodDef linearpartition_methods[] = {
    {"partition",   linearpartition_partition,  METH_VARARGS,
     linearpartition_partition_doc},
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

    return PyModuleDef_Init(&linearpartitionmodule);
}
