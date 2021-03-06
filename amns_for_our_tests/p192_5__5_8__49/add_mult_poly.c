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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 6979754601800116337UL) + ((((uint64_t)op[1] * 12268654575493002125UL) + ((uint64_t)op[2] * 11599715444872488152UL) + ((uint64_t)op[3] * 5624979254503744693UL) + ((uint64_t)op[4] * 15704748041887679244UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 15704748041887679244UL) + ((uint64_t)op[1] * 6979754601800116337UL) + ((((uint64_t)op[2] * 12268654575493002125UL) + ((uint64_t)op[3] * 11599715444872488152UL) + ((uint64_t)op[4] * 5624979254503744693UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 5624979254503744693UL) + ((uint64_t)op[1] * 15704748041887679244UL) + ((uint64_t)op[2] * 6979754601800116337UL) + ((((uint64_t)op[3] * 12268654575493002125UL) + ((uint64_t)op[4] * 11599715444872488152UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 11599715444872488152UL) + ((uint64_t)op[1] * 5624979254503744693UL) + ((uint64_t)op[2] * 15704748041887679244UL) + ((uint64_t)op[3] * 6979754601800116337UL) + ((uint64_t)op[4] * 5915516235396258920UL);
	tmp_q[4] = ((uint64_t)op[0] * 12268654575493002125UL) + ((uint64_t)op[1] * 11599715444872488152UL) + ((uint64_t)op[2] * 5624979254503744693UL) + ((uint64_t)op[3] * 15704748041887679244UL) + ((uint64_t)op[4] * 6979754601800116337UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 95068509263L) + ((((int128)tmp_q[1] * 84648364660L) + ((int128)tmp_q[2] * 848659912L) - ((int128)tmp_q[3] * 153582550779L) - ((int128)tmp_q[4] * 215858141340L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 215858141340L) + ((int128)tmp_q[1] * 95068509263L) + ((((int128)tmp_q[2] * 84648364660L) + ((int128)tmp_q[3] * 848659912L) - ((int128)tmp_q[4] * 153582550779L)) * 8);
	tmp_zero[2] = -((int128)tmp_q[0] * 153582550779L) - ((int128)tmp_q[1] * 215858141340L) + ((int128)tmp_q[2] * 95068509263L) + ((((int128)tmp_q[3] * 84648364660L) + ((int128)tmp_q[4] * 848659912L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 848659912L) - ((int128)tmp_q[1] * 153582550779L) - ((int128)tmp_q[2] * 215858141340L) + ((int128)tmp_q[3] * 95068509263L) + ((int128)tmp_q[4] * 677186917280L);
	tmp_zero[4] = ((int128)tmp_q[0] * 84648364660L) + ((int128)tmp_q[1] * 848659912L) - ((int128)tmp_q[2] * 153582550779L) - ((int128)tmp_q[3] * 215858141340L) + ((int128)tmp_q[4] * 95068509263L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

