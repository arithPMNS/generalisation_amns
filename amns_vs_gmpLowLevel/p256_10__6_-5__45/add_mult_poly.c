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
	tmp_q[0] = ((uint64_t)op[0] * 11060305342176205678UL) + ((((uint64_t)op[1] * 8325519388666372500UL) + ((uint64_t)op[2] * 7152994296373255593UL) + ((uint64_t)op[3] * 2087204625910139544UL) + ((uint64_t)op[4] * 17536542135617051365UL) + ((uint64_t)op[5] * 6263951399347846837UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 6263951399347846837UL) + ((uint64_t)op[1] * 11060305342176205678UL) + ((((uint64_t)op[2] * 8325519388666372500UL) + ((uint64_t)op[3] * 7152994296373255593UL) + ((uint64_t)op[4] * 2087204625910139544UL) + ((uint64_t)op[5] * 17536542135617051365UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 17536542135617051365UL) + ((uint64_t)op[1] * 6263951399347846837UL) + ((uint64_t)op[2] * 11060305342176205678UL) + ((((uint64_t)op[3] * 8325519388666372500UL) + ((uint64_t)op[4] * 7152994296373255593UL) + ((uint64_t)op[5] * 2087204625910139544UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 2087204625910139544UL) + ((uint64_t)op[1] * 17536542135617051365UL) + ((uint64_t)op[2] * 6263951399347846837UL) + ((uint64_t)op[3] * 11060305342176205678UL) + ((((uint64_t)op[4] * 8325519388666372500UL) + ((uint64_t)op[5] * 7152994296373255593UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 7152994296373255593UL) + ((uint64_t)op[1] * 2087204625910139544UL) + ((uint64_t)op[2] * 17536542135617051365UL) + ((uint64_t)op[3] * 6263951399347846837UL) + ((uint64_t)op[4] * 11060305342176205678UL) + ((uint64_t)op[5] * 13712635277796792348UL);
	tmp_q[5] = ((uint64_t)op[0] * 8325519388666372500UL) + ((uint64_t)op[1] * 7152994296373255593UL) + ((uint64_t)op[2] * 2087204625910139544UL) + ((uint64_t)op[3] * 17536542135617051365UL) + ((uint64_t)op[4] * 6263951399347846837UL) + ((uint64_t)op[5] * 11060305342176205678UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 3952223501529L) - ((-((int128)tmp_q[1] * 1398757316390L) - ((int128)tmp_q[2] * 1372567733785L) - ((int128)tmp_q[3] * 2129776524765L) + ((int128)tmp_q[4] * 125007874142L) + ((int128)tmp_q[5] * 4144166176190L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 4144166176190L) + ((int128)tmp_q[1] * 3952223501529L) - ((-((int128)tmp_q[2] * 1398757316390L) - ((int128)tmp_q[3] * 1372567733785L) - ((int128)tmp_q[4] * 2129776524765L) + ((int128)tmp_q[5] * 125007874142L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 125007874142L) + ((int128)tmp_q[1] * 4144166176190L) + ((int128)tmp_q[2] * 3952223501529L) - ((-((int128)tmp_q[3] * 1398757316390L) - ((int128)tmp_q[4] * 1372567733785L) - ((int128)tmp_q[5] * 2129776524765L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 2129776524765L) + ((int128)tmp_q[1] * 125007874142L) + ((int128)tmp_q[2] * 4144166176190L) + ((int128)tmp_q[3] * 3952223501529L) - ((-((int128)tmp_q[4] * 1398757316390L) - ((int128)tmp_q[5] * 1372567733785L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 1372567733785L) - ((int128)tmp_q[1] * 2129776524765L) + ((int128)tmp_q[2] * 125007874142L) + ((int128)tmp_q[3] * 4144166176190L) + ((int128)tmp_q[4] * 3952223501529L) + ((int128)tmp_q[5] * 6993786581950L);
	tmp_zero[5] = -((int128)tmp_q[0] * 1398757316390L) - ((int128)tmp_q[1] * 1372567733785L) - ((int128)tmp_q[2] * 2129776524765L) + ((int128)tmp_q[3] * 125007874142L) + ((int128)tmp_q[4] * 4144166176190L) + ((int128)tmp_q[5] * 3952223501529L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

