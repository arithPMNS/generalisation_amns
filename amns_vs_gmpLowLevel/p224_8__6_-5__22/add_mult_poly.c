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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 523237931511268145UL) + ((((uint64_t)op[1] * 4138601077840036877UL) + ((uint64_t)op[2] * 1578004128485388527UL) + ((uint64_t)op[3] * 16696583577785229245UL) + ((uint64_t)op[4] * 2162567755331498464UL) + ((uint64_t)op[5] * 3469539936740169139UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 3469539936740169139UL) + ((uint64_t)op[1] * 523237931511268145UL) + ((((uint64_t)op[2] * 4138601077840036877UL) + ((uint64_t)op[3] * 1578004128485388527UL) + ((uint64_t)op[4] * 16696583577785229245UL) + ((uint64_t)op[5] * 2162567755331498464UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 2162567755331498464UL) + ((uint64_t)op[1] * 3469539936740169139UL) + ((uint64_t)op[2] * 523237931511268145UL) + ((((uint64_t)op[3] * 4138601077840036877UL) + ((uint64_t)op[4] * 1578004128485388527UL) + ((uint64_t)op[5] * 16696583577785229245UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 16696583577785229245UL) + ((uint64_t)op[1] * 2162567755331498464UL) + ((uint64_t)op[2] * 3469539936740169139UL) + ((uint64_t)op[3] * 523237931511268145UL) + ((((uint64_t)op[4] * 4138601077840036877UL) + ((uint64_t)op[5] * 1578004128485388527UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 1578004128485388527UL) + ((uint64_t)op[1] * 16696583577785229245UL) + ((uint64_t)op[2] * 2162567755331498464UL) + ((uint64_t)op[3] * 3469539936740169139UL) + ((uint64_t)op[4] * 523237931511268145UL) + ((uint64_t)op[5] * 16200482758218918847UL);
	tmp_q[5] = ((uint64_t)op[0] * 4138601077840036877UL) + ((uint64_t)op[1] * 1578004128485388527UL) + ((uint64_t)op[2] * 16696583577785229245UL) + ((uint64_t)op[3] * 2162567755331498464UL) + ((uint64_t)op[4] * 3469539936740169139UL) + ((uint64_t)op[5] * 523237931511268145UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 15990742489L) - ((((int128)tmp_q[1] * 27309018979L) + ((int128)tmp_q[2] * 13389885616L) + ((int128)tmp_q[3] * 80672704221L) - ((int128)tmp_q[4] * 9032816811L) - ((int128)tmp_q[5] * 410328935L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 410328935L) + ((int128)tmp_q[1] * 15990742489L) - ((((int128)tmp_q[2] * 27309018979L) + ((int128)tmp_q[3] * 13389885616L) + ((int128)tmp_q[4] * 80672704221L) - ((int128)tmp_q[5] * 9032816811L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 9032816811L) - ((int128)tmp_q[1] * 410328935L) + ((int128)tmp_q[2] * 15990742489L) - ((((int128)tmp_q[3] * 27309018979L) + ((int128)tmp_q[4] * 13389885616L) + ((int128)tmp_q[5] * 80672704221L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 80672704221L) - ((int128)tmp_q[1] * 9032816811L) - ((int128)tmp_q[2] * 410328935L) + ((int128)tmp_q[3] * 15990742489L) - ((((int128)tmp_q[4] * 27309018979L) + ((int128)tmp_q[5] * 13389885616L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 13389885616L) + ((int128)tmp_q[1] * 80672704221L) - ((int128)tmp_q[2] * 9032816811L) - ((int128)tmp_q[3] * 410328935L) + ((int128)tmp_q[4] * 15990742489L) - ((int128)tmp_q[5] * 136545094895L);
	tmp_zero[5] = ((int128)tmp_q[0] * 27309018979L) + ((int128)tmp_q[1] * 13389885616L) + ((int128)tmp_q[2] * 80672704221L) - ((int128)tmp_q[3] * 9032816811L) - ((int128)tmp_q[4] * 410328935L) + ((int128)tmp_q[5] * 15990742489L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

