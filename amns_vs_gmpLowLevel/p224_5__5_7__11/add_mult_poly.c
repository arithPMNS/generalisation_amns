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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9645404766833929966UL) + ((((uint64_t)op[1] * 15405583628885701079UL) + ((uint64_t)op[2] * 1910615545346763121UL) + ((uint64_t)op[3] * 8690852494130369280UL) + ((uint64_t)op[4] * 10188863780493511609UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 10188863780493511609UL) + ((uint64_t)op[1] * 9645404766833929966UL) + ((((uint64_t)op[2] * 15405583628885701079UL) + ((uint64_t)op[3] * 1910615545346763121UL) + ((uint64_t)op[4] * 8690852494130369280UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 8690852494130369280UL) + ((uint64_t)op[1] * 10188863780493511609UL) + ((uint64_t)op[2] * 9645404766833929966UL) + ((((uint64_t)op[3] * 15405583628885701079UL) + ((uint64_t)op[4] * 1910615545346763121UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 1910615545346763121UL) + ((uint64_t)op[1] * 8690852494130369280UL) + ((uint64_t)op[2] * 10188863780493511609UL) + ((uint64_t)op[3] * 9645404766833929966UL) + ((uint64_t)op[4] * 15605365033652149473UL);
	tmp_q[4] = ((uint64_t)op[0] * 15405583628885701079UL) + ((uint64_t)op[1] * 1910615545346763121UL) + ((uint64_t)op[2] * 8690852494130369280UL) + ((uint64_t)op[3] * 10188863780493511609UL) + ((uint64_t)op[4] * 9645404766833929966UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 6603455902843L) + ((((int128)tmp_q[1] * 15021357045881L) + ((int128)tmp_q[2] * 10999067504481L) - ((int128)tmp_q[3] * 9905048540514L) - ((int128)tmp_q[4] * 16061599767332L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 16061599767332L) + ((int128)tmp_q[1] * 6603455902843L) + ((((int128)tmp_q[2] * 15021357045881L) + ((int128)tmp_q[3] * 10999067504481L) - ((int128)tmp_q[4] * 9905048540514L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 9905048540514L) - ((int128)tmp_q[1] * 16061599767332L) + ((int128)tmp_q[2] * 6603455902843L) + ((((int128)tmp_q[3] * 15021357045881L) + ((int128)tmp_q[4] * 10999067504481L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 10999067504481L) - ((int128)tmp_q[1] * 9905048540514L) - ((int128)tmp_q[2] * 16061599767332L) + ((int128)tmp_q[3] * 6603455902843L) + ((int128)tmp_q[4] * 105149499321167L);
	tmp_zero[4] = ((int128)tmp_q[0] * 15021357045881L) + ((int128)tmp_q[1] * 10999067504481L) - ((int128)tmp_q[2] * 9905048540514L) - ((int128)tmp_q[3] * 16061599767332L) + ((int128)tmp_q[4] * 6603455902843L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

