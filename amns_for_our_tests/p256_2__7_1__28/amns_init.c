#include "amns_init.h"


void init_data(){
	int i;
	for(i=0; i<POLY_DEG; i++)
		mpz_init (gama_pow[i]);

	mpz_init (modul_p);


	mpz_set_str (modul_p, "90625167144997146313007430523377040135861388266225474166759963580345530807309", 10);

	mpz_set_str (gama_pow[0], "72548496643481765512992903671589743454070297352744731396602810466087967718302", 10);
	for(i=1; i<POLY_DEG; i++){
		mpz_mul (gama_pow[i], gama_pow[i-1], gama_pow[0]);
		mpz_mod (gama_pow[i], gama_pow[i], modul_p);
	}

	//~ Note : lambda = 1 (see : modular multiplication in 'add_mult_poly.c').

	amns_rho = 1L << RHO_LOG2;

	//~ IMPORTANT : initialisations above must be done before those below.
	compute_rho_pows();
}


//~ computes representatives of powers of 'phi' in the AMNS.
void compute_rho_pows(){
	int i, l;
	int64_t tmp_rho[NB_COEFF];

	for(i=0; i<NB_COEFF; i++)
		tmp_rho[i] = rho_rep[i];

	//~ computation of a representative of 'rho'
	from_mont_domain(rho_rep, rho_rep);

	//~ computation of representatives of (rho)^i (for i=2,3,...)
	l = NB_COEFF - 2;
	if (l > 0){
		mult_mod_poly(RHO_POWS[0], rho_rep, tmp_rho);
		for(i=1; i<l; i++)
			mult_mod_poly(RHO_POWS[i], RHO_POWS[i-1], tmp_rho);
	}
}


void free_data(){
	int i;
	for(i=0; i<POLY_DEG; i++)
		mpz_clear (gama_pow[i]);

	mpz_clear (modul_p);
}

