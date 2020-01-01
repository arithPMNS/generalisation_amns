#include "amns_init.h"


void init_data(){
	int i;
	for(i=0; i<POLY_DEG; i++)
		mpz_init (gama_pow[i]);

	mpz_init (modul_p);


	mpz_set_str (modul_p, "78724239600625480298692691298392864706707889970105486261087202076552731179457", 10);

	mpz_set_str (gama_pow[0], "67121246142329733953157532273927750633628991646367117191031794526296140945147", 10);
	for(i=1; i<POLY_DEG; i++){
		mpz_mul (gama_pow[i], gama_pow[i-1], gama_pow[0]);
		mpz_mod (gama_pow[i], gama_pow[i], modul_p);
	}

	//~ IMPORTANT : initialisations above must be done before those below.
	compute_polys_P();
}


//~ computes representations of the polynomials P, used for conversion into the AMNS
void compute_polys_P(){
	int i, l;
	int64_t tmp_poly[NB_COEFF];

	//~ computation of a representation of 'phi*rho'
	from_mont_domain(tmp_poly, poly_P1);

	l = NB_COEFF - 2;
	if (l > 0){
		mult_mod_poly(polys_P[0], poly_P1, tmp_poly);
		for(i=1; i<l; i++)
			mult_mod_poly(polys_P[i], polys_P[i-1], tmp_poly);
	}
}


void free_data(){
	int i;
	for(i=0; i<POLY_DEG; i++)
		mpz_clear (gama_pow[i]);

	mpz_clear (modul_p);
}

