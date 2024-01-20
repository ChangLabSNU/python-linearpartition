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

struct basepair_prob {
    int32_t i;
    int32_t j;
    double prob;
};

static PyArray_Descr *partition_return_descr;

class MyBeamCKYParser : public BeamCKYParser {
public:
    using BeamCKYParser::BeamCKYParser;

    PyObject *
    get_basepair_prob(void)
    {
        PyArrayObject *res;
        npy_intp dim;
        struct basepair_prob *bpp;

        dim = Pij.size();

        res = (PyArrayObject *)PyArray_SimpleNewFromDescr(1, &dim, partition_return_descr);
        if (res == NULL)
            return NULL;
        Py_INCREF(partition_return_descr);

        assert(partition_return_descr->elsize == sizeof(struct basepair_prob));
        bpp = (struct basepair_prob *)PyArray_DATA(res);

        for (auto it = Pij.begin(); it != Pij.end(); ++it) {
            bpp->i = it->first.first - 1;
            bpp->j = it->first.second - 1;
            bpp->prob = it->second;
            bpp++;
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
linearpartition_partition(PyObject *self, PyObject *args, PyObject *kwds)
{
    const char *seq;
    Py_ssize_t len;
    int beamsize=100, dangles=2;
    static const char *kwlist[] = {"seq", "beamsize", "dangles", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s#|ii:partition", (char**)kwlist,
                                     &seq, &len, &beamsize, &dangles))
        return NULL;

    string rna_seq(seq);

    /* Call LinearPartition */
    MyBeamCKYParser parser(beamsize, true, false, "", "", false, 0.0, "",
                           false, 3.0, "", false, false, 0.3, "", "", false,
                           dangles);
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
