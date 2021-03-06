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
	tmp_q[0] = ((uint64_t)op[0] * 16669024764891734359UL) + ((((uint64_t)op[1] * 7914124323988538737UL) + ((uint64_t)op[2] * 14526944410034908394UL) + ((uint64_t)op[3] * 3721202232199656432UL) + ((uint64_t)op[4] * 15980713680128620233UL) + ((uint64_t)op[5] * 1513002437964802420UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 1513002437964802420UL) + ((uint64_t)op[1] * 16669024764891734359UL) + ((((uint64_t)op[2] * 7914124323988538737UL) + ((uint64_t)op[3] * 14526944410034908394UL) + ((uint64_t)op[4] * 3721202232199656432UL) + ((uint64_t)op[5] * 15980713680128620233UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 15980713680128620233UL) + ((uint64_t)op[1] * 1513002437964802420UL) + ((uint64_t)op[2] * 16669024764891734359UL) + ((((uint64_t)op[3] * 7914124323988538737UL) + ((uint64_t)op[4] * 14526944410034908394UL) + ((uint64_t)op[5] * 3721202232199656432UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 3721202232199656432UL) + ((uint64_t)op[1] * 15980713680128620233UL) + ((uint64_t)op[2] * 1513002437964802420UL) + ((uint64_t)op[3] * 16669024764891734359UL) + ((((uint64_t)op[4] * 7914124323988538737UL) + ((uint64_t)op[5] * 14526944410034908394UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 14526944410034908394UL) + ((uint64_t)op[1] * 3721202232199656432UL) + ((uint64_t)op[2] * 15980713680128620233UL) + ((uint64_t)op[3] * 1513002437964802420UL) + ((uint64_t)op[4] * 16669024764891734359UL) + ((uint64_t)op[5] * 5236990851464948284UL);
	tmp_q[5] = ((uint64_t)op[0] * 7914124323988538737UL) + ((uint64_t)op[1] * 14526944410034908394UL) + ((uint64_t)op[2] * 3721202232199656432UL) + ((uint64_t)op[3] * 15980713680128620233UL) + ((uint64_t)op[4] * 1513002437964802420UL) + ((uint64_t)op[5] * 16669024764891734359UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 20063782875L) - ((-((int128)tmp_q[1] * 78678711731L) - ((int128)tmp_q[2] * 50366749173L) - ((int128)tmp_q[3] * 74041783780L) + ((int128)tmp_q[4] * 79275749965L) + ((int128)tmp_q[5] * 2880411388L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 2880411388L) - ((int128)tmp_q[1] * 20063782875L) - ((-((int128)tmp_q[2] * 78678711731L) - ((int128)tmp_q[3] * 50366749173L) - ((int128)tmp_q[4] * 74041783780L) + ((int128)tmp_q[5] * 79275749965L)) * 4);
	tmp_zero[2] = ((int128)tmp_q[0] * 79275749965L) + ((int128)tmp_q[1] * 2880411388L) - ((int128)tmp_q[2] * 20063782875L) - ((-((int128)tmp_q[3] * 78678711731L) - ((int128)tmp_q[4] * 50366749173L) - ((int128)tmp_q[5] * 74041783780L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 74041783780L) + ((int128)tmp_q[1] * 79275749965L) + ((int128)tmp_q[2] * 2880411388L) - ((int128)tmp_q[3] * 20063782875L) - ((-((int128)tmp_q[4] * 78678711731L) - ((int128)tmp_q[5] * 50366749173L)) * 4);
	tmp_zero[4] = -((int128)tmp_q[0] * 50366749173L) - ((int128)tmp_q[1] * 74041783780L) + ((int128)tmp_q[2] * 79275749965L) + ((int128)tmp_q[3] * 2880411388L) - ((int128)tmp_q[4] * 20063782875L) + ((int128)tmp_q[5] * 314714846924L);
	tmp_zero[5] = -((int128)tmp_q[0] * 78678711731L) - ((int128)tmp_q[1] * 50366749173L) - ((int128)tmp_q[2] * 74041783780L) + ((int128)tmp_q[3] * 79275749965L) + ((int128)tmp_q[4] * 2880411388L) - ((int128)tmp_q[5] * 20063782875L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

