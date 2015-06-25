/*
Copyright (C) 2008-2012 Association of Universities for Research in Astronomy (AURA)

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    3. The name of AURA and its representatives may not be used to
      endorse or promote products derived from this software without
      specific prior written permission.

THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#define NO_IMPORT_ARRAY

#include "pyutil.h"

/***************************************************************************
 * List-of-strings proxy object
 ***************************************************************************/

static PyTypeObject PyStrListProxyType;

typedef struct {
  PyObject_HEAD
  /*@null@*/ /*@shared@*/ PyObject* pyobject;
  Py_ssize_t size;
  Py_ssize_t maxsize;
  char (*array)[72];
  str_verify_fn verify;
} PyStrListProxy;

static void
PyStrListProxy_dealloc(
    PyStrListProxy* self) {

  Py_XDECREF(self->pyobject);
  Py_TYPE(self)->tp_free((PyObject*)self);
}

/*@null@*/ static PyObject *
PyStrListProxy_new(
    PyTypeObject* type,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  PyStrListProxy* self = NULL;

  self = (PyStrListProxy*)type->tp_alloc(type, 0);
  if (self != NULL) {
    self->pyobject = NULL;
  }
  return (PyObject*)self;
}

static int
PyStrListProxy_traverse(
    PyStrListProxy* self,
    visitproc visit,
    void *arg) {

  int vret;

  if (self->pyobject) {
    vret = visit(self->pyobject, arg);
    if (vret != 0) {
      return vret;
    }
  }

  return 0;
}

static int
PyStrListProxy_clear(
    PyStrListProxy *self) {

  PyObject *tmp;

  tmp = self->pyobject;
  self->pyobject = NULL;
  Py_XDECREF(tmp);

  return 0;
}

/*@null@*/ PyObject *
PyStrListProxy_New(
    /*@shared@*/ PyObject* owner,
    Py_ssize_t size,
    Py_ssize_t maxsize,
    char (*array)[72],
    str_verify_fn verify) {

  PyStrListProxy* self = NULL;

  if (maxsize == 0) {
    maxsize = 68;
  }

  self = (PyStrListProxy*)PyStrListProxyType.tp_alloc(&PyStrListProxyType, 0);
  if (self == NULL) {
    return NULL;
  }

  Py_XINCREF(owner);
  self->pyobject = owner;
  self->size = size;
  self->maxsize = maxsize;
  self->array = array;
  self->verify = verify;
  return (PyObject*)self;
}

static Py_ssize_t
PyStrListProxy_len(
    PyStrListProxy* self) {

  return self->size;
}

/*@null@*/ static PyObject*
PyStrListProxy_getitem(
    PyStrListProxy* self,
    Py_ssize_t index) {

  if (index >= self->size) {
    PyErr_SetString(PyExc_IndexError, "index out of range");
    return NULL;
  }

  #if PY3K
  return PyBytes_FromString(self->array[index]);
  #else
  return PyString_FromString(self->array[index]);
  #endif
}

static int
PyStrListProxy_setitem(
    PyStrListProxy* self,
    Py_ssize_t index,
    PyObject* arg) {

  char* value;
  Py_ssize_t value_length;

  if (index > self->size) {
    PyErr_SetString(PyExc_IndexError, "index out of range");
    return -1;
  }

  #if PY3K
  if (PyBytes_AsStringAndSize(arg, &value, &value_length)) {
  #else
  if (PyString_AsStringAndSize(arg, &value, &value_length)) {
  #endif
      return -1;
  }

  if (value_length >= self->maxsize) {
    PyErr_Format(PyExc_ValueError,
                 "string must be less than %zd characters", self->maxsize);
    return -1;
  }

  if (self->verify && !self->verify(value)) {
    return -1;
  }

  strncpy(self->array[index], value, self->maxsize);

  return 0;
}

/*@null@*/ static PyObject*
PyStrListProxy_repr(
    PyStrListProxy* self) {

  char*       buffer  = NULL;
  char*       wp      = NULL;
  char*       rp      = NULL;
  Py_ssize_t  i       = 0;
  Py_ssize_t  j       = 0;
  PyObject*   result  = NULL;
  /* These are in descending order, so we can exit the loop quickly.  They
     are in pairs: (char_to_escape, char_escaped) */
  const char* escapes   = "\\\\''\rr\ff\vv\nn\tt\bb\aa";
  const char* e         = NULL;
  char        next_char = '\0';

  /* Overallocating to allow for escaped characters */
  buffer = malloc((size_t)self->size*self->maxsize*2 + 2);
  if (buffer == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Could not allocate memory.");
    return NULL;
  }

  wp = buffer;
  *wp++ = '[';

  for (i = 0; i < self->size; ++i) {
    *wp++ = '\'';
    rp = self->array[i];
    for (j = 0; j < self->maxsize && *rp != '\0'; ++j) {
      /* Check if this character should be escaped */
      e = escapes;
      next_char = *rp++;
      do {
        if (next_char > *e) {
          break;
        } else if (next_char == *e) {
          *wp++ = '\\';
          next_char = *(++e);
          break;
        } else {
          e += 2;
        }
      } while (*e != '\0');

      *wp++ = next_char;
    }
    *wp++ = '\'';

    /* Add a comma for all but the last one */
    if (i != self->size - 1) {
      *wp++ = ',';
      *wp++ = ' ';
    }
  }

  *wp++ = ']';
  *wp++ = '\0';

  #if PY3K
  result = PyUnicode_FromString(buffer);
  #else
  result = PyString_FromString(buffer);
  #endif
  free(buffer);
  return result;
}

static PySequenceMethods PyStrListProxy_sequence_methods = {
  (lenfunc)PyStrListProxy_len,
  NULL,
  NULL,
  (ssizeargfunc)PyStrListProxy_getitem,
  NULL,
  (ssizeobjargproc)PyStrListProxy_setitem,
  NULL,
  NULL,
  NULL,
  NULL
};

static PyTypeObject PyStrListProxyType = {
  #if PY3K
  PyVarObject_HEAD_INIT(NULL, 0)
  #else
  PyObject_HEAD_INIT(NULL)
  0,                          /*ob_size*/
  #endif
  "pywcs.StrListProxy",       /*tp_name*/
  sizeof(PyStrListProxy),  /*tp_basicsize*/
  0,                          /*tp_itemsize*/
  (destructor)PyStrListProxy_dealloc, /*tp_dealloc*/
  0,                          /*tp_print*/
  0,                          /*tp_getattr*/
  0,                          /*tp_setattr*/
  0,                          /*tp_compare*/
  (reprfunc)PyStrListProxy_repr, /*tp_repr*/
  0,                          /*tp_as_number*/
  &PyStrListProxy_sequence_methods, /*tp_as_sequence*/
  0,                          /*tp_as_mapping*/
  0,                          /*tp_hash */
  0,                          /*tp_call*/
  (reprfunc)PyStrListProxy_repr, /*tp_str*/
  0,                          /*tp_getattro*/
  0,                          /*tp_setattro*/
  0,                          /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC, /*tp_flags*/
  0,                          /* tp_doc */
  (traverseproc)PyStrListProxy_traverse, /* tp_traverse */
  (inquiry)PyStrListProxy_clear, /* tp_clear */
  0,                          /* tp_richcompare */
  0,                          /* tp_weaklistoffset */
  0,                          /* tp_iter */
  0,                          /* tp_iternext */
  0,                          /* tp_methods */
  0,                          /* tp_members */
  0,                          /* tp_getset */
  0,                          /* tp_base */
  0,                          /* tp_dict */
  0,                          /* tp_descr_get */
  0,                          /* tp_descr_set */
  0,                          /* tp_dictoffset */
  0,                          /* tp_init */
  0,                          /* tp_alloc */
  PyStrListProxy_new,      /* tp_new */
};

int
_setup_str_list_proxy_type(
    /*@unused@*/ PyObject* m) {

  if (PyType_Ready(&PyStrListProxyType) < 0) {
    return 1;
  }

  return 0;
}