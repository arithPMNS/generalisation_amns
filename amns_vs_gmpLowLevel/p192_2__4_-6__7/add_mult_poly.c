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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[3] + (int128)pa[3] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[3] * pa[1]) << 1) + (int128)pa[2] * pa[2]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[3] * pa[2]) * 12);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 7399373302649539965UL) + ((((uint64_t)op[1] * 9696381619417806478UL) + ((uint64_t)op[2] * 556469831145731517UL) + ((uint64_t)op[3] * 10173615455381375160UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 10173615455381375160UL) + ((uint64_t)op[1] * 7399373302649539965UL) + ((((uint64_t)op[2] * 9696381619417806478UL) + ((uint64_t)op[3] * 556469831145731517UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 556469831145731517UL) + ((uint64_t)op[1] * 10173615455381375160UL) + ((uint64_t)op[2] * 7399373302649539965UL) + ((uint64_t)op[3] * 15608686578331367596UL);
	tmp_q[3] = ((uint64_t)op[0] * 9696381619417806478UL) + ((uint64_t)op[1] * 556469831145731517UL) + ((uint64_t)op[2] * 10173615455381375160UL) + ((uint64_t)op[3] * 7399373302649539965UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 7264615008749L) - ((-((int128)tmp_q[1] * 88939616629286L) - ((int128)tmp_q[2] * 66657285036677L) - ((int128)tmp_q[3] * 226989354655632L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 226989354655632L) + ((int128)tmp_q[1] * 7264615008749L) - ((-((int128)tmp_q[2] * 88939616629286L) - ((int128)tmp_q[3] * 66657285036677L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 66657285036677L) - ((int128)tmp_q[1] * 226989354655632L) + ((int128)tmp_q[2] * 7264615008749L) + ((int128)tmp_q[3] * 533637699775716L);
	tmp_zero[3] = -((int128)tmp_q[0] * 88939616629286L) - ((int128)tmp_q[1] * 66657285036677L) - ((int128)tmp_q[2] * 226989354655632L) + ((int128)tmp_q[3] * 7264615008749L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
}

