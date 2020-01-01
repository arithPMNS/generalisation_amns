#ifndef STRUCTS_DATA
#define STRUCTS_DATA


//~ IMPORTANT : We take 'phi = 1 << WORD_SIZE'
#define WORD_SIZE 32
#define POLY_DEG 12
#define NB_COEFF 13
#define NB_ADD_MAX 0

#define RHO_LOG2 25  // rho = 1 << RHO_LOG2.

typedef long long llong;

//~ representations of the polynomials P0 and P1, used for conversion into the AMNS
int poly_P0[NB_COEFF] = {1378454, 1172524, 1271469, 1417701, 1762033, 1822691, 1765092, 998244, 1155434, 1269198, 954449, 927611, 815732};
int poly_P1[NB_COEFF] = {1452962, 1412803, 1575334, 1403022, 1887167, 1717534, 1778213, 1189339, 927788, 1395248, 1172295, 614777, 840824};

//~ representations of polynomials Pi, for i=2,...,n-1
int polys_P[(NB_COEFF - 2)][NB_COEFF];

mpz_t modul_p;
mpz_t gama_pow[POLY_DEG];

#endif

