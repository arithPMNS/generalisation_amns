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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 12);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 12);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17457719723006371951UL) + ((((uint64_t)op[1] * 4411737149905776881UL) + ((uint64_t)op[2] * 12141044045065805467UL) + ((uint64_t)op[3] * 3085019962685461437UL) + ((uint64_t)op[4] * 12554935446243464728UL) + ((uint64_t)op[5] * 14479965021734112814UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 14479965021734112814UL) + ((uint64_t)op[1] * 17457719723006371951UL) + ((((uint64_t)op[2] * 4411737149905776881UL) + ((uint64_t)op[3] * 12141044045065805467UL) + ((uint64_t)op[4] * 3085019962685461437UL) + ((uint64_t)op[5] * 12554935446243464728UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 12554935446243464728UL) + ((uint64_t)op[1] * 14479965021734112814UL) + ((uint64_t)op[2] * 17457719723006371951UL) + ((((uint64_t)op[3] * 4411737149905776881UL) + ((uint64_t)op[4] * 12141044045065805467UL) + ((uint64_t)op[5] * 3085019962685461437UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 3085019962685461437UL) + ((uint64_t)op[1] * 12554935446243464728UL) + ((uint64_t)op[2] * 14479965021734112814UL) + ((uint64_t)op[3] * 17457719723006371951UL) + ((((uint64_t)op[4] * 4411737149905776881UL) + ((uint64_t)op[5] * 12141044045065805467UL)) * 18446744073709551610);
	tmp_q[4] = ((uint64_t)op[0] * 12141044045065805467UL) + ((uint64_t)op[1] * 3085019962685461437UL) + ((uint64_t)op[2] * 12554935446243464728UL) + ((uint64_t)op[3] * 14479965021734112814UL) + ((uint64_t)op[4] * 17457719723006371951UL) + ((uint64_t)op[5] * 10423065247984441946UL);
	tmp_q[5] = ((uint64_t)op[0] * 4411737149905776881UL) + ((uint64_t)op[1] * 12141044045065805467UL) + ((uint64_t)op[2] * 3085019962685461437UL) + ((uint64_t)op[3] * 12554935446243464728UL) + ((uint64_t)op[4] * 14479965021734112814UL) + ((uint64_t)op[5] * 17457719723006371951UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 520552517L) - ((((int128)tmp_q[1] * 1854106277L) + ((int128)tmp_q[2] * 436854579L) - ((int128)tmp_q[3] * 1529184345L) + ((int128)tmp_q[4] * 2424554422L) + ((int128)tmp_q[5] * 1922647630L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 1922647630L) - ((int128)tmp_q[1] * 520552517L) - ((((int128)tmp_q[2] * 1854106277L) + ((int128)tmp_q[3] * 436854579L) - ((int128)tmp_q[4] * 1529184345L) + ((int128)tmp_q[5] * 2424554422L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 2424554422L) + ((int128)tmp_q[1] * 1922647630L) - ((int128)tmp_q[2] * 520552517L) - ((((int128)tmp_q[3] * 1854106277L) + ((int128)tmp_q[4] * 436854579L) - ((int128)tmp_q[5] * 1529184345L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 1529184345L) + ((int128)tmp_q[1] * 2424554422L) + ((int128)tmp_q[2] * 1922647630L) - ((int128)tmp_q[3] * 520552517L) - ((((int128)tmp_q[4] * 1854106277L) + ((int128)tmp_q[5] * 436854579L)) * 6);
	tmp_zero[4] = ((int128)tmp_q[0] * 436854579L) - ((int128)tmp_q[1] * 1529184345L) + ((int128)tmp_q[2] * 2424554422L) + ((int128)tmp_q[3] * 1922647630L) - ((int128)tmp_q[4] * 520552517L) - ((int128)tmp_q[5] * 11124637662L);
	tmp_zero[5] = ((int128)tmp_q[0] * 1854106277L) + ((int128)tmp_q[1] * 436854579L) - ((int128)tmp_q[2] * 1529184345L) + ((int128)tmp_q[3] * 2424554422L) + ((int128)tmp_q[4] * 1922647630L) - ((int128)tmp_q[5] * 520552517L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

