#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _isodistrreg_pavaCorrect(SEXP);
extern SEXP _isodistrreg_pavaDec(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_isodistrreg_pavaCorrect", (DL_FUNC) &_isodistrreg_pavaCorrect, 1},
    {"_isodistrreg_pavaDec",     (DL_FUNC) &_isodistrreg_pavaDec,     3},
    {NULL, NULL, 0}
};

void R_init_isodistrreg(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
