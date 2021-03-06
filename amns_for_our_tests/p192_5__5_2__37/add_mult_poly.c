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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9532067335759982199UL) + ((((uint64_t)op[1] * 4629830729692043233UL) + ((uint64_t)op[2] * 2397011310387355323UL) + ((uint64_t)op[3] * 10658673912707400852UL) + ((uint64_t)op[4] * 17978903265931653112UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 17978903265931653112UL) + ((uint64_t)op[1] * 9532067335759982199UL) + ((((uint64_t)op[2] * 4629830729692043233UL) + ((uint64_t)op[3] * 2397011310387355323UL) + ((uint64_t)op[4] * 10658673912707400852UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 10658673912707400852UL) + ((uint64_t)op[1] * 17978903265931653112UL) + ((uint64_t)op[2] * 9532067335759982199UL) + ((((uint64_t)op[3] * 4629830729692043233UL) + ((uint64_t)op[4] * 2397011310387355323UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 2397011310387355323UL) + ((uint64_t)op[1] * 10658673912707400852UL) + ((uint64_t)op[2] * 17978903265931653112UL) + ((uint64_t)op[3] * 9532067335759982199UL) + ((uint64_t)op[4] * 9259661459384086466UL);
	tmp_q[4] = ((uint64_t)op[0] * 4629830729692043233UL) + ((uint64_t)op[1] * 2397011310387355323UL) + ((uint64_t)op[2] * 10658673912707400852UL) + ((uint64_t)op[3] * 17978903265931653112UL) + ((uint64_t)op[4] * 9532067335759982199UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 157925223427L) + ((((int128)tmp_q[1] * 210399272471L) + ((int128)tmp_q[2] * 117394992429L) + ((int128)tmp_q[3] * 89728184392L) + ((int128)tmp_q[4] * 61694466958L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 61694466958L) - ((int128)tmp_q[1] * 157925223427L) + ((((int128)tmp_q[2] * 210399272471L) + ((int128)tmp_q[3] * 117394992429L) + ((int128)tmp_q[4] * 89728184392L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 89728184392L) + ((int128)tmp_q[1] * 61694466958L) - ((int128)tmp_q[2] * 157925223427L) + ((((int128)tmp_q[3] * 210399272471L) + ((int128)tmp_q[4] * 117394992429L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 117394992429L) + ((int128)tmp_q[1] * 89728184392L) + ((int128)tmp_q[2] * 61694466958L) - ((int128)tmp_q[3] * 157925223427L) + ((int128)tmp_q[4] * 420798544942L);
	tmp_zero[4] = ((int128)tmp_q[0] * 210399272471L) + ((int128)tmp_q[1] * 117394992429L) + ((int128)tmp_q[2] * 89728184392L) + ((int128)tmp_q[3] * 61694466958L) - ((int128)tmp_q[4] * 157925223427L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

