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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[3] + (int128)pa[3] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[3] * pa[1]) << 1) + (int128)pa[2] * pa[2]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[3] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 2123063455193643631UL) + ((((uint64_t)op[1] * 3251580039755050679UL) + ((uint64_t)op[2] * 7601314161907600500UL) + ((uint64_t)op[3] * 8021143697540587289UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 8021143697540587289UL) + ((uint64_t)op[1] * 2123063455193643631UL) + ((((uint64_t)op[2] * 3251580039755050679UL) + ((uint64_t)op[3] * 7601314161907600500UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 7601314161907600500UL) + ((uint64_t)op[1] * 8021143697540587289UL) + ((uint64_t)op[2] * 2123063455193643631UL) + ((uint64_t)op[3] * 4314316204575803137UL);
	tmp_q[3] = ((uint64_t)op[0] * 3251580039755050679UL) + ((uint64_t)op[1] * 7601314161907600500UL) + ((uint64_t)op[2] * 8021143697540587289UL) + ((uint64_t)op[3] * 2123063455193643631UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 83616374387919L) + ((-((int128)tmp_q[1] * 101258035838919L) - ((int128)tmp_q[2] * 27881365505404L) + ((int128)tmp_q[3] * 178915845697687L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 178915845697687L) + ((int128)tmp_q[1] * 83616374387919L) + ((-((int128)tmp_q[2] * 101258035838919L) - ((int128)tmp_q[3] * 27881365505404L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 27881365505404L) + ((int128)tmp_q[1] * 178915845697687L) + ((int128)tmp_q[2] * 83616374387919L) - ((int128)tmp_q[3] * 708806250872433L);
	tmp_zero[3] = -((int128)tmp_q[0] * 101258035838919L) - ((int128)tmp_q[1] * 27881365505404L) + ((int128)tmp_q[2] * 178915845697687L) + ((int128)tmp_q[3] * 83616374387919L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
}

