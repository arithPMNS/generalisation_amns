#ifndef STRUCTS_DATA
#define STRUCTS_DATA


//~ IMPORTANT : We take 'phi = 1 << WORD_SIZE'
#define WORD_SIZE 64
#define POLY_DEG 4
#define NB_COEFF 5
#define NB_ADD_MAX 4

#define RHO_LOG2 55  // rho = 1 << RHO_LOG2.

typedef __int128 int128;

//~ representations of the polynomials P0 and P1, used for conversion into the AMNS
int64_t poly_P0[NB_COEFF] = {2417154222381215, 1717597245556995, 2230246623216472, 1342091056533950, 1871634882583606};
int64_t poly_P1[NB_COEFF] = {3124006495515283, 1604926006855910, 2387750801860216, 1034119891388689, 1305308913034411};

//~ representations of polynomials Pi, for i=2,...,n-1
int64_t polys_P[(NB_COEFF - 2)][NB_COEFF];

mpz_t modul_p;
mpz_t gama_pow[POLY_DEG];

#endif

