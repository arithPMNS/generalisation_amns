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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8498372873996857538UL) + ((((uint64_t)op[1] * 3026624438781307316UL) + ((uint64_t)op[2] * 12851941307027432127UL) + ((uint64_t)op[3] * 5473122098184222973UL) + ((uint64_t)op[4] * 9845139156838266927UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 9845139156838266927UL) + ((uint64_t)op[1] * 8498372873996857538UL) + ((((uint64_t)op[2] * 3026624438781307316UL) + ((uint64_t)op[3] * 12851941307027432127UL) + ((uint64_t)op[4] * 5473122098184222973UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 5473122098184222973UL) + ((uint64_t)op[1] * 9845139156838266927UL) + ((uint64_t)op[2] * 8498372873996857538UL) + ((((uint64_t)op[3] * 3026624438781307316UL) + ((uint64_t)op[4] * 12851941307027432127UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 12851941307027432127UL) + ((uint64_t)op[1] * 5473122098184222973UL) + ((uint64_t)op[2] * 9845139156838266927UL) + ((uint64_t)op[3] * 8498372873996857538UL) + ((uint64_t)op[4] * 3313621879803015036UL);
	tmp_q[4] = ((uint64_t)op[0] * 3026624438781307316UL) + ((uint64_t)op[1] * 12851941307027432127UL) + ((uint64_t)op[2] * 5473122098184222973UL) + ((uint64_t)op[3] * 9845139156838266927UL) + ((uint64_t)op[4] * 8498372873996857538UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1120795650036437L) - ((-((int128)tmp_q[1] * 1666820078597186L) - ((int128)tmp_q[2] * 1055910103388135L) + ((int128)tmp_q[3] * 21215254189604L) - ((int128)tmp_q[4] * 274089830677205L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 274089830677205L) + ((int128)tmp_q[1] * 1120795650036437L) - ((-((int128)tmp_q[2] * 1666820078597186L) - ((int128)tmp_q[3] * 1055910103388135L) + ((int128)tmp_q[4] * 21215254189604L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 21215254189604L) - ((int128)tmp_q[1] * 274089830677205L) + ((int128)tmp_q[2] * 1120795650036437L) - ((-((int128)tmp_q[3] * 1666820078597186L) - ((int128)tmp_q[4] * 1055910103388135L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 1055910103388135L) + ((int128)tmp_q[1] * 21215254189604L) - ((int128)tmp_q[2] * 274089830677205L) + ((int128)tmp_q[3] * 1120795650036437L) + ((int128)tmp_q[4] * 8334100392985930L);
	tmp_zero[4] = -((int128)tmp_q[0] * 1666820078597186L) - ((int128)tmp_q[1] * 1055910103388135L) + ((int128)tmp_q[2] * 21215254189604L) - ((int128)tmp_q[3] * 274089830677205L) + ((int128)tmp_q[4] * 1120795650036437L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

