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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 12964538047280455440UL) + ((((uint64_t)op[1] * 756347366088876916UL) + ((uint64_t)op[2] * 11934202902924023557UL) + ((uint64_t)op[3] * 4154485431099718184UL) + ((uint64_t)op[4] * 1709890508680016728UL) + ((uint64_t)op[5] * 3635370807865044720UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 3635370807865044720UL) + ((uint64_t)op[1] * 12964538047280455440UL) + ((((uint64_t)op[2] * 756347366088876916UL) + ((uint64_t)op[3] * 11934202902924023557UL) + ((uint64_t)op[4] * 4154485431099718184UL) + ((uint64_t)op[5] * 1709890508680016728UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 1709890508680016728UL) + ((uint64_t)op[1] * 3635370807865044720UL) + ((uint64_t)op[2] * 12964538047280455440UL) + ((((uint64_t)op[3] * 756347366088876916UL) + ((uint64_t)op[4] * 11934202902924023557UL) + ((uint64_t)op[5] * 4154485431099718184UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 4154485431099718184UL) + ((uint64_t)op[1] * 1709890508680016728UL) + ((uint64_t)op[2] * 3635370807865044720UL) + ((uint64_t)op[3] * 12964538047280455440UL) + ((((uint64_t)op[4] * 756347366088876916UL) + ((uint64_t)op[5] * 11934202902924023557UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 11934202902924023557UL) + ((uint64_t)op[1] * 4154485431099718184UL) + ((uint64_t)op[2] * 1709890508680016728UL) + ((uint64_t)op[3] * 3635370807865044720UL) + ((uint64_t)op[4] * 12964538047280455440UL) + ((uint64_t)op[5] * 16177701975442920868UL);
	tmp_q[5] = ((uint64_t)op[0] * 756347366088876916UL) + ((uint64_t)op[1] * 11934202902924023557UL) + ((uint64_t)op[2] * 4154485431099718184UL) + ((uint64_t)op[3] * 1709890508680016728UL) + ((uint64_t)op[4] * 3635370807865044720UL) + ((uint64_t)op[5] * 12964538047280455440UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 863878712L) - ((((int128)tmp_q[1] * 1463514672L) - ((int128)tmp_q[2] * 489091392L) + ((int128)tmp_q[3] * 1274017700L) - ((int128)tmp_q[4] * 1594745809L) + ((int128)tmp_q[5] * 1184269000L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 1184269000L) + ((int128)tmp_q[1] * 863878712L) - ((((int128)tmp_q[2] * 1463514672L) - ((int128)tmp_q[3] * 489091392L) + ((int128)tmp_q[4] * 1274017700L) - ((int128)tmp_q[5] * 1594745809L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 1594745809L) + ((int128)tmp_q[1] * 1184269000L) + ((int128)tmp_q[2] * 863878712L) - ((((int128)tmp_q[3] * 1463514672L) - ((int128)tmp_q[4] * 489091392L) + ((int128)tmp_q[5] * 1274017700L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 1274017700L) - ((int128)tmp_q[1] * 1594745809L) + ((int128)tmp_q[2] * 1184269000L) + ((int128)tmp_q[3] * 863878712L) - ((((int128)tmp_q[4] * 1463514672L) - ((int128)tmp_q[5] * 489091392L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 489091392L) + ((int128)tmp_q[1] * 1274017700L) - ((int128)tmp_q[2] * 1594745809L) + ((int128)tmp_q[3] * 1184269000L) + ((int128)tmp_q[4] * 863878712L) - ((int128)tmp_q[5] * 4390544016L);
	tmp_zero[5] = ((int128)tmp_q[0] * 1463514672L) - ((int128)tmp_q[1] * 489091392L) + ((int128)tmp_q[2] * 1274017700L) - ((int128)tmp_q[3] * 1594745809L) + ((int128)tmp_q[4] * 1184269000L) + ((int128)tmp_q[5] * 863878712L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

