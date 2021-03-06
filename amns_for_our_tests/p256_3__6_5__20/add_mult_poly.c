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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 13919076647584448179UL) + ((((uint64_t)op[1] * 12702490863403369954UL) + ((uint64_t)op[2] * 5913686514441635457UL) + ((uint64_t)op[3] * 10231470964538927109UL) + ((uint64_t)op[4] * 16711766198510297840UL) + ((uint64_t)op[5] * 2840063755512190688UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 2840063755512190688UL) + ((uint64_t)op[1] * 13919076647584448179UL) + ((((uint64_t)op[2] * 12702490863403369954UL) + ((uint64_t)op[3] * 5913686514441635457UL) + ((uint64_t)op[4] * 10231470964538927109UL) + ((uint64_t)op[5] * 16711766198510297840UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 16711766198510297840UL) + ((uint64_t)op[1] * 2840063755512190688UL) + ((uint64_t)op[2] * 13919076647584448179UL) + ((((uint64_t)op[3] * 12702490863403369954UL) + ((uint64_t)op[4] * 5913686514441635457UL) + ((uint64_t)op[5] * 10231470964538927109UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 10231470964538927109UL) + ((uint64_t)op[1] * 16711766198510297840UL) + ((uint64_t)op[2] * 2840063755512190688UL) + ((uint64_t)op[3] * 13919076647584448179UL) + ((((uint64_t)op[4] * 12702490863403369954UL) + ((uint64_t)op[5] * 5913686514441635457UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 5913686514441635457UL) + ((uint64_t)op[1] * 10231470964538927109UL) + ((uint64_t)op[2] * 16711766198510297840UL) + ((uint64_t)op[3] * 2840063755512190688UL) + ((uint64_t)op[4] * 13919076647584448179UL) + ((uint64_t)op[5] * 8172222095888194922UL);
	tmp_q[5] = ((uint64_t)op[0] * 12702490863403369954UL) + ((uint64_t)op[1] * 5913686514441635457UL) + ((uint64_t)op[2] * 10231470964538927109UL) + ((uint64_t)op[3] * 16711766198510297840UL) + ((uint64_t)op[4] * 2840063755512190688UL) + ((uint64_t)op[5] * 13919076647584448179UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1379340792658L) + ((-((int128)tmp_q[1] * 2629088917404L) + ((int128)tmp_q[2] * 2907619501939L) - ((int128)tmp_q[3] * 361994755736L) - ((int128)tmp_q[4] * 5333338706209L) - ((int128)tmp_q[5] * 3719095313539L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 3719095313539L) - ((int128)tmp_q[1] * 1379340792658L) + ((-((int128)tmp_q[2] * 2629088917404L) + ((int128)tmp_q[3] * 2907619501939L) - ((int128)tmp_q[4] * 361994755736L) - ((int128)tmp_q[5] * 5333338706209L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 5333338706209L) - ((int128)tmp_q[1] * 3719095313539L) - ((int128)tmp_q[2] * 1379340792658L) + ((-((int128)tmp_q[3] * 2629088917404L) + ((int128)tmp_q[4] * 2907619501939L) - ((int128)tmp_q[5] * 361994755736L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 361994755736L) - ((int128)tmp_q[1] * 5333338706209L) - ((int128)tmp_q[2] * 3719095313539L) - ((int128)tmp_q[3] * 1379340792658L) + ((-((int128)tmp_q[4] * 2629088917404L) + ((int128)tmp_q[5] * 2907619501939L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 2907619501939L) - ((int128)tmp_q[1] * 361994755736L) - ((int128)tmp_q[2] * 5333338706209L) - ((int128)tmp_q[3] * 3719095313539L) - ((int128)tmp_q[4] * 1379340792658L) - ((int128)tmp_q[5] * 13145444587020L);
	tmp_zero[5] = -((int128)tmp_q[0] * 2629088917404L) + ((int128)tmp_q[1] * 2907619501939L) - ((int128)tmp_q[2] * 361994755736L) - ((int128)tmp_q[3] * 5333338706209L) - ((int128)tmp_q[4] * 3719095313539L) - ((int128)tmp_q[5] * 1379340792658L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

