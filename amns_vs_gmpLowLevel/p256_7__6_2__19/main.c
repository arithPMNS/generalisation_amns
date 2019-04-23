#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <gmp.h>

#include "structs_data.h"
#include "add_mult_poly.c"
#include "useful_functs.c"
#include "amns_init.c"

#define BILLION 1000000000L


//~ Compilation command : gcc -Wall -O3 main.c -o main -lgmp
//~ Execution command : ./main

//~ Important : polynomials representations form is P(X) = a0 + ... + an.X^n = (a0, ..., an).


int main(void){

	int i, nbiter;
	mpz_t A, B, C;
	mpz_inits (A, B, C, NULL);


	int nb_limbs;
	int64_t pa[NB_COEFF];
	int64_t pb[NB_COEFF];
	mp_limb_t *p_limbs, *a_limbs, *b_limbs, *c_limbs, *q_limbs;

	struct timespec start1, end1;
	struct timespec start2, end2;
	struct timespec start3, end3;
	uint64_t diff1, diff2, diff3;

	unsigned long seed = time(NULL);
	gmp_randstate_t r;
	gmp_randinit_default(r);
	gmp_randseed_ui(r, seed);


	init_data();
	
	nb_limbs = mpz_size (modul_p);
	
	q_limbs = (mp_limb_t*) calloc ((nb_limbs+1), sizeof(mp_limb_t));
	c_limbs = (mp_limb_t*) calloc ((nb_limbs*2), sizeof(mp_limb_t));


	mpz_urandomm(A, r, modul_p);
	mpz_urandomm(B, r, modul_p);
	
	p_limbs = mpz_limbs_modify (modul_p, nb_limbs);
	a_limbs = mpz_limbs_modify (A, nb_limbs);
	b_limbs = mpz_limbs_modify (B, nb_limbs);


	nbiter = 1 << 25;

	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start1);
	for (i=0; i<nbiter; i++) {
		mpn_mul_n (c_limbs, a_limbs, b_limbs, nb_limbs); // compute: z = y*x
		mpn_tdiv_qr (q_limbs, a_limbs, 0, c_limbs, (nb_limbs*2), p_limbs, nb_limbs); // compute: y = z%p
	}
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end1);
	diff1 = BILLION * (end1.tv_sec - start1.tv_sec) + (end1.tv_nsec - start1.tv_nsec);
	
	
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start3);
	from_int_to_amns(pa, A);
	from_int_to_amns(pb, B);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end3);
	diff3 = BILLION * (end3.tv_sec - start3.tv_sec) + (end3.tv_nsec - start3.tv_nsec);


	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start2);
	for (i=0; i<nbiter; i++) {
		mult_mod_poly(pa, pa, pb);
	}
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end2);
	diff2 = BILLION * (end2.tv_sec - start2.tv_sec) + (end2.tv_nsec - start2.tv_nsec);


	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start3);
	from_amns_to_int(C, pa);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end3);
	diff3 += BILLION * (end3.tv_sec - start3.tv_sec) + (end3.tv_nsec - start3.tv_nsec);


	printf("%lu %d %lu %lu %lf\n\n", mpz_sizeinbase (modul_p, 2), NB_COEFF, (diff2+diff3), diff1, (double)(diff2+diff3)/diff1);

	printf("nbiter = %d\n\n", nbiter);
	
	printf("\ntime using gmp_low	= %lu nanoseconds\n", diff1);
	printf("\ntime using amns prod	= %lu nanoseconds\n", (diff2+diff3));
	printf("\nratio			= %lf\n", (double)(diff2+diff3)/diff1);

	
	free(c_limbs);
	mpz_clears (A, B, C, NULL);
	gmp_randclear(r);


	free_data();
	return 0;
}

