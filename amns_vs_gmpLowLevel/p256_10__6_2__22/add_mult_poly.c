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
	tmp_q[0] = ((uint64_t)op[0] * 14919950344199431827UL) + ((((uint64_t)op[1] * 10998056177567464328UL) + ((uint64_t)op[2] * 3738710252817854208UL) + ((uint64_t)op[3] * 8175673011684429558UL) + ((uint64_t)op[4] * 13655148318591634010UL) + ((uint64_t)op[5] * 2549890466626672535UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 2549890466626672535UL) + ((uint64_t)op[1] * 14919950344199431827UL) + ((((uint64_t)op[2] * 10998056177567464328UL) + ((uint64_t)op[3] * 3738710252817854208UL) + ((uint64_t)op[4] * 8175673011684429558UL) + ((uint64_t)op[5] * 13655148318591634010UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 13655148318591634010UL) + ((uint64_t)op[1] * 2549890466626672535UL) + ((uint64_t)op[2] * 14919950344199431827UL) + ((((uint64_t)op[3] * 10998056177567464328UL) + ((uint64_t)op[4] * 3738710252817854208UL) + ((uint64_t)op[5] * 8175673011684429558UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 8175673011684429558UL) + ((uint64_t)op[1] * 13655148318591634010UL) + ((uint64_t)op[2] * 2549890466626672535UL) + ((uint64_t)op[3] * 14919950344199431827UL) + ((((uint64_t)op[4] * 10998056177567464328UL) + ((uint64_t)op[5] * 3738710252817854208UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 3738710252817854208UL) + ((uint64_t)op[1] * 8175673011684429558UL) + ((uint64_t)op[2] * 13655148318591634010UL) + ((uint64_t)op[3] * 2549890466626672535UL) + ((uint64_t)op[4] * 14919950344199431827UL) + ((uint64_t)op[5] * 3549368281425377040UL);
	tmp_q[5] = ((uint64_t)op[0] * 10998056177567464328UL) + ((uint64_t)op[1] * 3738710252817854208UL) + ((uint64_t)op[2] * 8175673011684429558UL) + ((uint64_t)op[3] * 13655148318591634010UL) + ((uint64_t)op[4] * 2549890466626672535UL) + ((uint64_t)op[5] * 14919950344199431827UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 971833314881L) + ((-((int128)tmp_q[1] * 4313126417261L) + ((int128)tmp_q[2] * 1870305650869L) - ((int128)tmp_q[3] * 4061697788497L) + ((int128)tmp_q[4] * 145229857233L) - ((int128)tmp_q[5] * 3073360413675L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 3073360413675L) - ((int128)tmp_q[1] * 971833314881L) + ((-((int128)tmp_q[2] * 4313126417261L) + ((int128)tmp_q[3] * 1870305650869L) - ((int128)tmp_q[4] * 4061697788497L) + ((int128)tmp_q[5] * 145229857233L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 145229857233L) - ((int128)tmp_q[1] * 3073360413675L) - ((int128)tmp_q[2] * 971833314881L) + ((-((int128)tmp_q[3] * 4313126417261L) + ((int128)tmp_q[4] * 1870305650869L) - ((int128)tmp_q[5] * 4061697788497L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 4061697788497L) + ((int128)tmp_q[1] * 145229857233L) - ((int128)tmp_q[2] * 3073360413675L) - ((int128)tmp_q[3] * 971833314881L) + ((-((int128)tmp_q[4] * 4313126417261L) + ((int128)tmp_q[5] * 1870305650869L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 1870305650869L) - ((int128)tmp_q[1] * 4061697788497L) + ((int128)tmp_q[2] * 145229857233L) - ((int128)tmp_q[3] * 3073360413675L) - ((int128)tmp_q[4] * 971833314881L) - ((int128)tmp_q[5] * 8626252834522L);
	tmp_zero[5] = -((int128)tmp_q[0] * 4313126417261L) + ((int128)tmp_q[1] * 1870305650869L) - ((int128)tmp_q[2] * 4061697788497L) + ((int128)tmp_q[3] * 145229857233L) - ((int128)tmp_q[4] * 3073360413675L) - ((int128)tmp_q[5] * 971833314881L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

