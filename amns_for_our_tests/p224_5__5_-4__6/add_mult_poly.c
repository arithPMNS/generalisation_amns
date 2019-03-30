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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 16398651064455239047UL) + ((((uint64_t)op[1] * 17557939556453240517UL) + ((uint64_t)op[2] * 11865656144459585195UL) + ((uint64_t)op[3] * 7307527870069863074UL) + ((uint64_t)op[4] * 560116234617608061UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 560116234617608061UL) + ((uint64_t)op[1] * 16398651064455239047UL) + ((((uint64_t)op[2] * 17557939556453240517UL) + ((uint64_t)op[3] * 11865656144459585195UL) + ((uint64_t)op[4] * 7307527870069863074UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 7307527870069863074UL) + ((uint64_t)op[1] * 560116234617608061UL) + ((uint64_t)op[2] * 16398651064455239047UL) + ((((uint64_t)op[3] * 17557939556453240517UL) + ((uint64_t)op[4] * 11865656144459585195UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 11865656144459585195UL) + ((uint64_t)op[1] * 7307527870069863074UL) + ((uint64_t)op[2] * 560116234617608061UL) + ((uint64_t)op[3] * 16398651064455239047UL) + ((uint64_t)op[4] * 3555218069025244396UL);
	tmp_q[4] = ((uint64_t)op[0] * 17557939556453240517UL) + ((uint64_t)op[1] * 11865656144459585195UL) + ((uint64_t)op[2] * 7307527870069863074UL) + ((uint64_t)op[3] * 560116234617608061UL) + ((uint64_t)op[4] * 16398651064455239047UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 3627804929385L) - ((-((int128)tmp_q[1] * 12882736031498L) + ((int128)tmp_q[2] * 13538367139072L) - ((int128)tmp_q[3] * 2290271083953L) + ((int128)tmp_q[4] * 16289326583409L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 16289326583409L) + ((int128)tmp_q[1] * 3627804929385L) - ((-((int128)tmp_q[2] * 12882736031498L) + ((int128)tmp_q[3] * 13538367139072L) - ((int128)tmp_q[4] * 2290271083953L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 2290271083953L) + ((int128)tmp_q[1] * 16289326583409L) + ((int128)tmp_q[2] * 3627804929385L) - ((-((int128)tmp_q[3] * 12882736031498L) + ((int128)tmp_q[4] * 13538367139072L)) * 4);
	tmp_zero[3] = ((int128)tmp_q[0] * 13538367139072L) - ((int128)tmp_q[1] * 2290271083953L) + ((int128)tmp_q[2] * 16289326583409L) + ((int128)tmp_q[3] * 3627804929385L) + ((int128)tmp_q[4] * 51530944125992L);
	tmp_zero[4] = -((int128)tmp_q[0] * 12882736031498L) + ((int128)tmp_q[1] * 13538367139072L) - ((int128)tmp_q[2] * 2290271083953L) + ((int128)tmp_q[3] * 16289326583409L) + ((int128)tmp_q[4] * 3627804929385L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

