#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void ascov(void *, void *, void *, void *, void *, void *, void *);
extern void ascov_all(void *, void *, void *, void *, void *, void *);
extern void ascov_deflij(void *, void *, void *, void *, void *, void *);
extern void ascov_deflji(void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"ascov",        (DL_FUNC) &ascov,        7},
    {"ascov_all",    (DL_FUNC) &ascov_all,    6},
    {"ascov_deflij", (DL_FUNC) &ascov_deflij, 6},
    {"ascov_deflji", (DL_FUNC) &ascov_deflji, 6},
    {NULL, NULL, 0}
};

void R_init_BSSasymp(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
