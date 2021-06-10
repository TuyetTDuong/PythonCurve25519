/* tell python that PyArg_ParseTuple(t#) means Py_ssize_t, not int */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#if (PY_VERSION_HEX < 0x02050000)
	typedef int Py_ssize_t;
#endif

/* This is required for compatibility with Python 2. */
#if PY_MAJOR_VERSION >= 3
	#include <bytesobject.h>
	#define y "y"
#else
	#define PyBytes_FromStringAndSize PyString_FromStringAndSize
	#define y "t"
#endif

extern void curve25519_dh_CalculatePublicKey(uint8_t *mypublic,const uint8_t *secret);
extern void curve25519_dh_CalculateSumTwoPublicKeys(uint8_t *mysumpublic, uint8_t *publicP, uint8_t *publicQ);
extern void curve25519_dh_CreateSharedKey(uint8_t *mysumpublic, const uint8_t *publicP, uint8_t *publicQ);
extern void curve25519_dh_CalculateSumFieldElements(uint8_t *sum, uint8_t *element1, uint8_t *element2);
extern void curve25519_dh_CalculateProductFieldElements(uint8_t *sum, uint8_t *element1, uint8_t *element2);
static PyObject *
generatePrivateKey(PyObject *self, PyObject *args)
{
    char *random;
    Py_ssize_t randomlen;

    if(!PyArg_ParseTuple(args, y"#:clamp", &random, &randomlen)) {
        return NULL;
    }

    if(randomlen != 32) {
        PyErr_SetString(PyExc_ValueError, "random must be 32-byte string");
        return NULL;
    }
    random[0] &= 248;
    random[31] &= 127;
    random[31] |= 64;

    return PyBytes_FromStringAndSize((char *)random, 32);
}

static PyObject *
generatePublicKey(PyObject *self, PyObject *args)
{
    const char *private;
    char mypublic[32];
    char basepoint[32] = {9};
    Py_ssize_t privatelen;
    if (!PyArg_ParseTuple(args, y"#:makepublic", &private, &privatelen))
        return NULL;
    if (privatelen != 32) {
        PyErr_SetString(PyExc_ValueError, "input must be 32-byte string");
        return NULL;
    }
    curve25519_dh_CalculatePublicKey(mypublic, private);
    return PyBytes_FromStringAndSize((char *)mypublic, 32);
}


static PyObject *
addTwoPubKeys(PyObject *self, PyObject *args)
{
    const char *leftbytes, *parentbytes;

    char sum[32];
    Py_ssize_t leftbyteslen, parentbyteslen;
    if (!PyArg_ParseTuple(args, y"#"y"#:generate",
                          &leftbytes, &leftbyteslen, &parentbytes, &parentbyteslen))
        return NULL;
    if (leftbyteslen != 32) {
        PyErr_SetString(PyExc_ValueError, "input must be 32-byte string");
        return NULL;
    }
    if (parentbyteslen != 32) {
        PyErr_SetString(PyExc_ValueError, "input must be 32-byte string");
        return NULL;
    }
    curve25519_dh_CalculateSumTwoPublicKeys(sum, leftbytes, parentbytes);

    return PyBytes_FromStringAndSize((char *)sum, 32);
}
static PyObject *addTwoFieldElements(PyObject *self, PyObject *args)
{
    const char *element1, *element2;

    char sum[32];
    Py_ssize_t element1len, element2len;
    if (!PyArg_ParseTuple(args, y"#"y"#:generate",
                          &element1, &element1len, &element2, &element2len))
        return NULL;
    if (element1len != 32) {
        PyErr_SetString(PyExc_ValueError, "input must be 32-byte string");
        return NULL;
    }
    if (element2len != 32) {
        PyErr_SetString(PyExc_ValueError, "input must be 32-byte string");
        return NULL;
    }
    curve25519_dh_CalculateSumFieldElements(sum, element1, element2);

    return PyBytes_FromStringAndSize((char *)sum, 32);
}
static PyObject *mulTwoFieldElements(PyObject *self, PyObject *args)
{
    const char *element1, *element2;

    char sum[32];
    Py_ssize_t element1len, element2len;
    if (!PyArg_ParseTuple(args, y"#"y"#:generate",
                          &element1, &element1len, &element2, &element2len))
        return NULL;
    if (element1len != 32) {
        PyErr_SetString(PyExc_ValueError, "input must be 32-byte string");
        return NULL;
    }
    if (element2len != 32) {
        PyErr_SetString(PyExc_ValueError, "input must be 32-byte string");
        return NULL;
    }
    curve25519_dh_CalculateProductFieldElements(sum, element1, element2);

    return PyBytes_FromStringAndSize((char *)sum, 32);
}
static PyMethodDef
curve25519_functions[] = {
    {"private", generatePrivateKey, METH_VARARGS, "data->private"},
    {"public", generatePublicKey, METH_VARARGS, "private->public"},
    {"addTwoPubKeys", addTwoPubKeys, METH_VARARGS, "pointA+pointB->sum"},
    {"addTwoFieldElements", addTwoFieldElements, METH_VARARGS, "pointA+pointB->sum"},
    {"mulTwoFieldElements", mulTwoFieldElements, METH_VARARGS, "pointA+pointB->sum"},
    {NULL, NULL, 0, NULL},
};


#if PY_MAJOR_VERSION >= 3
    static struct PyModuleDef
    curve25519_module = {
        PyModuleDef_HEAD_INIT,
        "curve25519_topl",
        NULL,
        NULL,
        curve25519_functions,
    };

    PyObject *
    PyInit_curve25519_topl(void)
    {
        return PyModule_Create(&curve25519_module);
    }
#else

    PyMODINIT_FUNC
    initcurve25519_topl(void)
    {
        (void)Py_InitModule("curve25519_topl", curve25519_functions);
    }

#endif
