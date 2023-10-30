#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <structure_paralleltestbed.h>
#include <python_interface.h>

int init_python(const char * string1, const char * string2, void *data){

     PyObject *pModule, *pFunc, *pName;

     experiment_total *exp1;

     exp1 = (experiment_total *) data;


     Py_Initialize();
     PyRun_SimpleString("import sys");
     PyRun_SimpleString("sys.path.append( \"./benchmarks/python\" )");

     pName = PyUnicode_DecodeFSDefault(string1);

     if (pName == NULL) printf("error null pname \n");

     pModule = PyImport_Import(pName);

     Py_DECREF(pName);

    if (pModule != NULL) {

       pFunc = PyObject_GetAttrString(pModule, string2);
    }
    else {

        PyErr_Print();
        fprintf(stderr, "Failed to load \"%s\"\n", string1);
    }
    
    exp1->execution.pFunc = pFunc;
    exp1->execution.pModule = pModule;
//    exp1->execution.gstate = (PyGILState_STATE *) malloc(1 * sizeof(PyGILState_STATE));

    return 1;
}


void* call_to_obj_function_python(double *x, void *data) {

    double fx = 1e10;
    int i;
    PyObject *pArgs, *pValue;
    PyObject *pFunc;
    int size_x;
    output_function *res;
    experiment_total *exp1;

//    printf("call to obj function python\n");
    exp1 = (experiment_total *) data;
    pFunc = exp1->execution.pFunc;

//    PyGILState_STATE *tstate;

//    *tstate= PyGILState_Ensure();

    size_x=(*exp1).test.bench.dim;
    res = NULL;
    res = (output_function *) calloc(1,sizeof(output_function));

    /* pFunc is a new reference */
    if (pFunc && PyCallable_Check(pFunc)) {
        pArgs = PyTuple_New(1);
	    // X
	    PyObject *plist1 = PyList_New(size_x);
	    for (i=0;i<size_x;i++) {
            	pValue = PyFloat_FromDouble(x[i]);
	    	PyList_SetItem(plist1, i, pValue);
	    }
        PyTuple_SetItem(pArgs, 0, plist1);
		
        pValue = PyObject_CallObject(pFunc, pArgs);
        Py_DECREF(pArgs);
	Py_DECREF(plist1);

        if (pValue != NULL) {
	    fx = PyFloat_AsDouble(pValue);
            Py_DECREF(pValue);
        }

     }
     else {
        if (PyErr_Occurred())
            PyErr_Print();
        return res;
    }

//    PyGILState_Release(*tstate);

    res->value = fx;
//    printf("%lf\n", fx);
    return res;
}



int end_python( void *data ){
    PyObject *pModule, *pFunc;

    experiment_total *exp1;

    exp1 = (experiment_total *) data;

    pFunc = exp1->execution.pFunc;
    pModule = exp1->execution.pModule;

    Py_XDECREF(pFunc);
    Py_DECREF(pModule);

    if (Py_FinalizeEx() < 0) {
       return 120;
    }
     return 0;
}

