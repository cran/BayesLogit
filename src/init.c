#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .C calls */
extern void ar1_llh(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void combine(void *, void *, void *, void *, void *);
extern void draw_indicators_generic(void *, void *, void *, void *, void *, void *, void *);
extern void draw_indicators_logistic(void *, void *, void *, void *, void *, void *, void *);
extern void draw_omega(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void draw_stc_beta(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void EM(void *, void *, void *, void *, void *, void *, void *, void *);
extern void ffbs(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gibbs(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void mult_combine(void *, void *, void *, void *, void *, void *);
extern void mult_gibbs(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rpg_alt(void *, void *, void *, void *);
extern void rpg_devroye(void *, void *, void *, void *);
extern void rpg_gamma(void *, void *, void *, void *, void *);
extern void rpg_hybrid(void *, void *, void *, void *);
extern void rpg_sp(void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    // {"ar1_llh",                  (DL_FUNC) &ar1_llh,                  11}, // not included
    {"combine",                  (DL_FUNC) &combine,                   5},
    {"draw_indicators_generic",  (DL_FUNC) &draw_indicators_generic,   7},
    {"draw_indicators_logistic", (DL_FUNC) &draw_indicators_logistic,  7},
    // {"draw_omega",               (DL_FUNC) &draw_omega,               22}, // not included
    // {"draw_stc_beta",            (DL_FUNC) &draw_stc_beta,            19}, // not included
    {"EM",                       (DL_FUNC) &EM,                        8},
    // {"ffbs",                     (DL_FUNC) &ffbs,                     14}, // not included
    {"gibbs",                    (DL_FUNC) &gibbs,                    11},
    {"mult_combine",             (DL_FUNC) &mult_combine,              6},
    {"mult_gibbs",               (DL_FUNC) &mult_gibbs,               12},
    {"rpg_alt",                  (DL_FUNC) &rpg_alt,                   4},
    {"rpg_devroye",              (DL_FUNC) &rpg_devroye,               4},
    {"rpg_gamma",                (DL_FUNC) &rpg_gamma,                 5},
    {"rpg_hybrid",               (DL_FUNC) &rpg_hybrid,                4},
    {"rpg_sp",                   (DL_FUNC) &rpg_sp,                    5},
    {NULL, NULL, 0}
};

void R_init_BayesLogit(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
