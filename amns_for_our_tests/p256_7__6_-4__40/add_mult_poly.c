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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 15115987933955654561UL) + ((((uint64_t)op[1] * 17742519405405275086UL) + ((uint64_t)op[2] * 14306888772530416663UL) + ((uint64_t)op[3] * 15903573978590718595UL) + ((uint64_t)op[4] * 9976239351840727943UL) + ((uint64_t)op[5] * 10759206284277806979UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 10759206284277806979UL) + ((uint64_t)op[1] * 15115987933955654561UL) + ((((uint64_t)op[2] * 17742519405405275086UL) + ((uint64_t)op[3] * 14306888772530416663UL) + ((uint64_t)op[4] * 15903573978590718595UL) + ((uint64_t)op[5] * 9976239351840727943UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 9976239351840727943UL) + ((uint64_t)op[1] * 10759206284277806979UL) + ((uint64_t)op[2] * 15115987933955654561UL) + ((((uint64_t)op[3] * 17742519405405275086UL) + ((uint64_t)op[4] * 14306888772530416663UL) + ((uint64_t)op[5] * 15903573978590718595UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 15903573978590718595UL) + ((uint64_t)op[1] * 9976239351840727943UL) + ((uint64_t)op[2] * 10759206284277806979UL) + ((uint64_t)op[3] * 15115987933955654561UL) + ((((uint64_t)op[4] * 17742519405405275086UL) + ((uint64_t)op[5] * 14306888772530416663UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 14306888772530416663UL) + ((uint64_t)op[1] * 15903573978590718595UL) + ((uint64_t)op[2] * 9976239351840727943UL) + ((uint64_t)op[3] * 10759206284277806979UL) + ((uint64_t)op[4] * 15115987933955654561UL) + ((uint64_t)op[5] * 2816898673217106120UL);
	tmp_q[5] = ((uint64_t)op[0] * 17742519405405275086UL) + ((uint64_t)op[1] * 14306888772530416663UL) + ((uint64_t)op[2] * 15903573978590718595UL) + ((uint64_t)op[3] * 9976239351840727943UL) + ((uint64_t)op[4] * 10759206284277806979UL) + ((uint64_t)op[5] * 15115987933955654561UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 920775363459L) - ((((int128)tmp_q[1] * 2512523328151L) + ((int128)tmp_q[2] * 3980909414548L) + ((int128)tmp_q[3] * 2652892840516L) + ((int128)tmp_q[4] * 2193274373550L) + ((int128)tmp_q[5] * 2090478455051L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 2090478455051L) + ((int128)tmp_q[1] * 920775363459L) - ((((int128)tmp_q[2] * 2512523328151L) + ((int128)tmp_q[3] * 3980909414548L) + ((int128)tmp_q[4] * 2652892840516L) + ((int128)tmp_q[5] * 2193274373550L)) * 4);
	tmp_zero[2] = ((int128)tmp_q[0] * 2193274373550L) + ((int128)tmp_q[1] * 2090478455051L) + ((int128)tmp_q[2] * 920775363459L) - ((((int128)tmp_q[3] * 2512523328151L) + ((int128)tmp_q[4] * 3980909414548L) + ((int128)tmp_q[5] * 2652892840516L)) * 4);
	tmp_zero[3] = ((int128)tmp_q[0] * 2652892840516L) + ((int128)tmp_q[1] * 2193274373550L) + ((int128)tmp_q[2] * 2090478455051L) + ((int128)tmp_q[3] * 920775363459L) - ((((int128)tmp_q[4] * 2512523328151L) + ((int128)tmp_q[5] * 3980909414548L)) * 4);
	tmp_zero[4] = ((int128)tmp_q[0] * 3980909414548L) + ((int128)tmp_q[1] * 2652892840516L) + ((int128)tmp_q[2] * 2193274373550L) + ((int128)tmp_q[3] * 2090478455051L) + ((int128)tmp_q[4] * 920775363459L) - ((int128)tmp_q[5] * 10050093312604L);
	tmp_zero[5] = ((int128)tmp_q[0] * 2512523328151L) + ((int128)tmp_q[1] * 3980909414548L) + ((int128)tmp_q[2] * 2652892840516L) + ((int128)tmp_q[3] * 2193274373550L) + ((int128)tmp_q[4] * 2090478455051L) + ((int128)tmp_q[5] * 920775363459L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

