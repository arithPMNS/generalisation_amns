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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 12687561557760135985UL) + ((((uint64_t)op[1] * 8211973297694143997UL) + ((uint64_t)op[2] * 2282297830457177227UL) + ((uint64_t)op[3] * 14058784841623740212UL) + ((uint64_t)op[4] * 6410084851094084583UL) + ((uint64_t)op[5] * 9065461491947900103UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 9065461491947900103UL) + ((uint64_t)op[1] * 12687561557760135985UL) + ((((uint64_t)op[2] * 8211973297694143997UL) + ((uint64_t)op[3] * 2282297830457177227UL) + ((uint64_t)op[4] * 14058784841623740212UL) + ((uint64_t)op[5] * 6410084851094084583UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 6410084851094084583UL) + ((uint64_t)op[1] * 9065461491947900103UL) + ((uint64_t)op[2] * 12687561557760135985UL) + ((((uint64_t)op[3] * 8211973297694143997UL) + ((uint64_t)op[4] * 2282297830457177227UL) + ((uint64_t)op[5] * 14058784841623740212UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 14058784841623740212UL) + ((uint64_t)op[1] * 6410084851094084583UL) + ((uint64_t)op[2] * 9065461491947900103UL) + ((uint64_t)op[3] * 12687561557760135985UL) + ((((uint64_t)op[4] * 8211973297694143997UL) + ((uint64_t)op[5] * 2282297830457177227UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 2282297830457177227UL) + ((uint64_t)op[1] * 14058784841623740212UL) + ((uint64_t)op[2] * 6410084851094084583UL) + ((uint64_t)op[3] * 9065461491947900103UL) + ((uint64_t)op[4] * 12687561557760135985UL) + ((uint64_t)op[5] * 2022797478321263622UL);
	tmp_q[5] = ((uint64_t)op[0] * 8211973297694143997UL) + ((uint64_t)op[1] * 2282297830457177227UL) + ((uint64_t)op[2] * 14058784841623740212UL) + ((uint64_t)op[3] * 6410084851094084583UL) + ((uint64_t)op[4] * 9065461491947900103UL) + ((uint64_t)op[5] * 12687561557760135985UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1458336960269L) - ((((int128)tmp_q[1] * 2666293583807L) + ((int128)tmp_q[2] * 468407816070L) - ((int128)tmp_q[3] * 1686966865711L) + ((int128)tmp_q[4] * 1820830798948L) - ((int128)tmp_q[5] * 3558601687093L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 3558601687093L) - ((int128)tmp_q[1] * 1458336960269L) - ((((int128)tmp_q[2] * 2666293583807L) + ((int128)tmp_q[3] * 468407816070L) - ((int128)tmp_q[4] * 1686966865711L) + ((int128)tmp_q[5] * 1820830798948L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 1820830798948L) - ((int128)tmp_q[1] * 3558601687093L) - ((int128)tmp_q[2] * 1458336960269L) - ((((int128)tmp_q[3] * 2666293583807L) + ((int128)tmp_q[4] * 468407816070L) - ((int128)tmp_q[5] * 1686966865711L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 1686966865711L) + ((int128)tmp_q[1] * 1820830798948L) - ((int128)tmp_q[2] * 3558601687093L) - ((int128)tmp_q[3] * 1458336960269L) - ((((int128)tmp_q[4] * 2666293583807L) + ((int128)tmp_q[5] * 468407816070L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 468407816070L) - ((int128)tmp_q[1] * 1686966865711L) + ((int128)tmp_q[2] * 1820830798948L) - ((int128)tmp_q[3] * 3558601687093L) - ((int128)tmp_q[4] * 1458336960269L) - ((int128)tmp_q[5] * 5332587167614L);
	tmp_zero[5] = ((int128)tmp_q[0] * 2666293583807L) + ((int128)tmp_q[1] * 468407816070L) - ((int128)tmp_q[2] * 1686966865711L) + ((int128)tmp_q[3] * 1820830798948L) - ((int128)tmp_q[4] * 3558601687093L) - ((int128)tmp_q[5] * 1458336960269L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

