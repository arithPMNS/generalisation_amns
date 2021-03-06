#include "add_mult_poly.h"


void add_poly(int64_t *rop, int64_t *pa, int64_t *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] + pb[j];
}

void sub_poly(int64_t *rop, int64_t *pa, int64_t *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] - pb[j];
}

void neg_poly(int64_t *rop, int64_t *op){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = -op[j];
}

//~ assumes 'scalar' and/or coeffs of 'op' small enough to avoid an overflow.
void scalar_mult_poly(int64_t *rop, int64_t *op, int64_t scalar){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = scalar * op[j];
}

//~ Computes pa(X)*pb(X) mod(X^n - c)
void mult_mod_poly(int64_t *rop, int64_t *pa, int64_t *pb){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 41927707469162263UL) + ((((uint64_t)op[1] * 8645687795459869980UL) + ((uint64_t)op[2] * 18418081090531272131UL) + ((uint64_t)op[3] * 14000976652896641067UL) + ((uint64_t)op[4] * 12510145157123786978UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 12510145157123786978UL) + ((uint64_t)op[1] * 41927707469162263UL) + ((((uint64_t)op[2] * 8645687795459869980UL) + ((uint64_t)op[3] * 18418081090531272131UL) + ((uint64_t)op[4] * 14000976652896641067UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 14000976652896641067UL) + ((uint64_t)op[1] * 12510145157123786978UL) + ((uint64_t)op[2] * 41927707469162263UL) + ((((uint64_t)op[3] * 8645687795459869980UL) + ((uint64_t)op[4] * 18418081090531272131UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 18418081090531272131UL) + ((uint64_t)op[1] * 14000976652896641067UL) + ((uint64_t)op[2] * 12510145157123786978UL) + ((uint64_t)op[3] * 41927707469162263UL) + ((uint64_t)op[4] * 7490319312670058324UL);
	tmp_q[4] = ((uint64_t)op[0] * 8645687795459869980UL) + ((uint64_t)op[1] * 18418081090531272131UL) + ((uint64_t)op[2] * 14000976652896641067UL) + ((uint64_t)op[3] * 12510145157123786978UL) + ((uint64_t)op[4] * 41927707469162263UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1928769942123551L) + ((((int128)tmp_q[1] * 578815723136903L) - ((int128)tmp_q[2] * 336716278798018L) + ((int128)tmp_q[3] * 543544617305826L) - ((int128)tmp_q[4] * 226511629592445L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 226511629592445L) + ((int128)tmp_q[1] * 1928769942123551L) + ((((int128)tmp_q[2] * 578815723136903L) - ((int128)tmp_q[3] * 336716278798018L) + ((int128)tmp_q[4] * 543544617305826L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 543544617305826L) - ((int128)tmp_q[1] * 226511629592445L) + ((int128)tmp_q[2] * 1928769942123551L) + ((((int128)tmp_q[3] * 578815723136903L) - ((int128)tmp_q[4] * 336716278798018L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 336716278798018L) + ((int128)tmp_q[1] * 543544617305826L) - ((int128)tmp_q[2] * 226511629592445L) + ((int128)tmp_q[3] * 1928769942123551L) + ((int128)tmp_q[4] * 1736447169410709L);
	tmp_zero[4] = ((int128)tmp_q[0] * 578815723136903L) - ((int128)tmp_q[1] * 336716278798018L) + ((int128)tmp_q[2] * 543544617305826L) - ((int128)tmp_q[3] * 226511629592445L) + ((int128)tmp_q[4] * 1928769942123551L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

