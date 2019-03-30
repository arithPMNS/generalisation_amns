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
	tmp_q[0] = ((uint64_t)op[0] * 242605730256542133UL) + ((((uint64_t)op[1] * 9011479981342725627UL) + ((uint64_t)op[2] * 847921297101743733UL) + ((uint64_t)op[3] * 18106276210236792569UL) + ((uint64_t)op[4] * 10766709172276137184UL) + ((uint64_t)op[5] * 6354290425853340169UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 6354290425853340169UL) + ((uint64_t)op[1] * 242605730256542133UL) + ((((uint64_t)op[2] * 9011479981342725627UL) + ((uint64_t)op[3] * 847921297101743733UL) + ((uint64_t)op[4] * 18106276210236792569UL) + ((uint64_t)op[5] * 10766709172276137184UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 10766709172276137184UL) + ((uint64_t)op[1] * 6354290425853340169UL) + ((uint64_t)op[2] * 242605730256542133UL) + ((((uint64_t)op[3] * 9011479981342725627UL) + ((uint64_t)op[4] * 847921297101743733UL) + ((uint64_t)op[5] * 18106276210236792569UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 18106276210236792569UL) + ((uint64_t)op[1] * 10766709172276137184UL) + ((uint64_t)op[2] * 6354290425853340169UL) + ((uint64_t)op[3] * 242605730256542133UL) + ((((uint64_t)op[4] * 9011479981342725627UL) + ((uint64_t)op[5] * 847921297101743733UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 847921297101743733UL) + ((uint64_t)op[1] * 18106276210236792569UL) + ((uint64_t)op[2] * 10766709172276137184UL) + ((uint64_t)op[3] * 6354290425853340169UL) + ((uint64_t)op[4] * 242605730256542133UL) + ((uint64_t)op[5] * 847568222048200724UL);
	tmp_q[5] = ((uint64_t)op[0] * 9011479981342725627UL) + ((uint64_t)op[1] * 847921297101743733UL) + ((uint64_t)op[2] * 18106276210236792569UL) + ((uint64_t)op[3] * 10766709172276137184UL) + ((uint64_t)op[4] * 6354290425853340169UL) + ((uint64_t)op[5] * 242605730256542133UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 77900188569L) - ((-((int128)tmp_q[1] * 5015350967L) - ((int128)tmp_q[2] * 99412835622L) - ((int128)tmp_q[3] * 11314313074L) - ((int128)tmp_q[4] * 69745948417L) - ((int128)tmp_q[5] * 7578245119L)) * 4);
	tmp_zero[1] = -((int128)tmp_q[0] * 7578245119L) - ((int128)tmp_q[1] * 77900188569L) - ((-((int128)tmp_q[2] * 5015350967L) - ((int128)tmp_q[3] * 99412835622L) - ((int128)tmp_q[4] * 11314313074L) - ((int128)tmp_q[5] * 69745948417L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 69745948417L) - ((int128)tmp_q[1] * 7578245119L) - ((int128)tmp_q[2] * 77900188569L) - ((-((int128)tmp_q[3] * 5015350967L) - ((int128)tmp_q[4] * 99412835622L) - ((int128)tmp_q[5] * 11314313074L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 11314313074L) - ((int128)tmp_q[1] * 69745948417L) - ((int128)tmp_q[2] * 7578245119L) - ((int128)tmp_q[3] * 77900188569L) - ((-((int128)tmp_q[4] * 5015350967L) - ((int128)tmp_q[5] * 99412835622L)) * 4);
	tmp_zero[4] = -((int128)tmp_q[0] * 99412835622L) - ((int128)tmp_q[1] * 11314313074L) - ((int128)tmp_q[2] * 69745948417L) - ((int128)tmp_q[3] * 7578245119L) - ((int128)tmp_q[4] * 77900188569L) + ((int128)tmp_q[5] * 20061403868L);
	tmp_zero[5] = -((int128)tmp_q[0] * 5015350967L) - ((int128)tmp_q[1] * 99412835622L) - ((int128)tmp_q[2] * 11314313074L) - ((int128)tmp_q[3] * 69745948417L) - ((int128)tmp_q[4] * 7578245119L) - ((int128)tmp_q[5] * 77900188569L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

