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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 16936283791067095213UL) + ((((uint64_t)op[1] * 17627454774789180721UL) + ((uint64_t)op[2] * 15707752331465963490UL) + ((uint64_t)op[3] * 844687849023206254UL) + ((uint64_t)op[4] * 2581914686047080725UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 2581914686047080725UL) + ((uint64_t)op[1] * 16936283791067095213UL) + ((((uint64_t)op[2] * 17627454774789180721UL) + ((uint64_t)op[3] * 15707752331465963490UL) + ((uint64_t)op[4] * 844687849023206254UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 844687849023206254UL) + ((uint64_t)op[1] * 2581914686047080725UL) + ((uint64_t)op[2] * 16936283791067095213UL) + ((((uint64_t)op[3] * 17627454774789180721UL) + ((uint64_t)op[4] * 15707752331465963490UL)) * 18446744073709551609);
	tmp_q[3] = ((uint64_t)op[0] * 15707752331465963490UL) + ((uint64_t)op[1] * 844687849023206254UL) + ((uint64_t)op[2] * 2581914686047080725UL) + ((uint64_t)op[3] * 16936283791067095213UL) + ((uint64_t)op[4] * 5735025092442596265UL);
	tmp_q[4] = ((uint64_t)op[0] * 17627454774789180721UL) + ((uint64_t)op[1] * 15707752331465963490UL) + ((uint64_t)op[2] * 844687849023206254UL) + ((uint64_t)op[3] * 2581914686047080725UL) + ((uint64_t)op[4] * 16936283791067095213UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 9515995970051L) - ((-((int128)tmp_q[1] * 10658669759320L) - ((int128)tmp_q[2] * 2519606656801L) - ((int128)tmp_q[3] * 9805933392781L) + ((int128)tmp_q[4] * 9257866839064L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 9257866839064L) + ((int128)tmp_q[1] * 9515995970051L) - ((-((int128)tmp_q[2] * 10658669759320L) - ((int128)tmp_q[3] * 2519606656801L) - ((int128)tmp_q[4] * 9805933392781L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 9805933392781L) + ((int128)tmp_q[1] * 9257866839064L) + ((int128)tmp_q[2] * 9515995970051L) - ((-((int128)tmp_q[3] * 10658669759320L) - ((int128)tmp_q[4] * 2519606656801L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 2519606656801L) - ((int128)tmp_q[1] * 9805933392781L) + ((int128)tmp_q[2] * 9257866839064L) + ((int128)tmp_q[3] * 9515995970051L) + ((int128)tmp_q[4] * 74610688315240L);
	tmp_zero[4] = -((int128)tmp_q[0] * 10658669759320L) - ((int128)tmp_q[1] * 2519606656801L) - ((int128)tmp_q[2] * 9805933392781L) + ((int128)tmp_q[3] * 9257866839064L) + ((int128)tmp_q[4] * 9515995970051L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

