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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11278657610187139869UL) + ((((uint64_t)op[1] * 2865830637571099385UL) + ((uint64_t)op[2] * 2667697366537990637UL) + ((uint64_t)op[3] * 14857454676137577653UL) + ((uint64_t)op[4] * 1576080279011637355UL) + ((uint64_t)op[5] * 12058308162857972107UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 12058308162857972107UL) + ((uint64_t)op[1] * 11278657610187139869UL) + ((((uint64_t)op[2] * 2865830637571099385UL) + ((uint64_t)op[3] * 2667697366537990637UL) + ((uint64_t)op[4] * 14857454676137577653UL) + ((uint64_t)op[5] * 1576080279011637355UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 1576080279011637355UL) + ((uint64_t)op[1] * 12058308162857972107UL) + ((uint64_t)op[2] * 11278657610187139869UL) + ((((uint64_t)op[3] * 2865830637571099385UL) + ((uint64_t)op[4] * 2667697366537990637UL) + ((uint64_t)op[5] * 14857454676137577653UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 14857454676137577653UL) + ((uint64_t)op[1] * 1576080279011637355UL) + ((uint64_t)op[2] * 12058308162857972107UL) + ((uint64_t)op[3] * 11278657610187139869UL) + ((((uint64_t)op[4] * 2865830637571099385UL) + ((uint64_t)op[5] * 2667697366537990637UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 2667697366537990637UL) + ((uint64_t)op[1] * 14857454676137577653UL) + ((uint64_t)op[2] * 1576080279011637355UL) + ((uint64_t)op[3] * 12058308162857972107UL) + ((uint64_t)op[4] * 11278657610187139869UL) + ((uint64_t)op[5] * 5731661275142198770UL);
	tmp_q[5] = ((uint64_t)op[0] * 2865830637571099385UL) + ((uint64_t)op[1] * 2667697366537990637UL) + ((uint64_t)op[2] * 14857454676137577653UL) + ((uint64_t)op[3] * 1576080279011637355UL) + ((uint64_t)op[4] * 12058308162857972107UL) + ((uint64_t)op[5] * 11278657610187139869UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1295096451L) + ((-((int128)tmp_q[1] * 695661660L) - ((int128)tmp_q[2] * 1111365954L) + ((int128)tmp_q[3] * 2519168654L) + ((int128)tmp_q[4] * 1010436614L) + ((int128)tmp_q[5] * 2060387529L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 2060387529L) - ((int128)tmp_q[1] * 1295096451L) + ((-((int128)tmp_q[2] * 695661660L) - ((int128)tmp_q[3] * 1111365954L) + ((int128)tmp_q[4] * 2519168654L) + ((int128)tmp_q[5] * 1010436614L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 1010436614L) + ((int128)tmp_q[1] * 2060387529L) - ((int128)tmp_q[2] * 1295096451L) + ((-((int128)tmp_q[3] * 695661660L) - ((int128)tmp_q[4] * 1111365954L) + ((int128)tmp_q[5] * 2519168654L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 2519168654L) + ((int128)tmp_q[1] * 1010436614L) + ((int128)tmp_q[2] * 2060387529L) - ((int128)tmp_q[3] * 1295096451L) + ((-((int128)tmp_q[4] * 695661660L) - ((int128)tmp_q[5] * 1111365954L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 1111365954L) + ((int128)tmp_q[1] * 2519168654L) + ((int128)tmp_q[2] * 1010436614L) + ((int128)tmp_q[3] * 2060387529L) - ((int128)tmp_q[4] * 1295096451L) - ((int128)tmp_q[5] * 1391323320L);
	tmp_zero[5] = -((int128)tmp_q[0] * 695661660L) - ((int128)tmp_q[1] * 1111365954L) + ((int128)tmp_q[2] * 2519168654L) + ((int128)tmp_q[3] * 1010436614L) + ((int128)tmp_q[4] * 2060387529L) - ((int128)tmp_q[5] * 1295096451L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

