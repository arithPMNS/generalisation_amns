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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17938432657587224280UL) + ((((uint64_t)op[1] * 3057272475926735445UL) + ((uint64_t)op[2] * 12472208059277854851UL) + ((uint64_t)op[3] * 15697300765962808534UL) + ((uint64_t)op[4] * 15422499593974599487UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 15422499593974599487UL) + ((uint64_t)op[1] * 17938432657587224280UL) + ((((uint64_t)op[2] * 3057272475926735445UL) + ((uint64_t)op[3] * 12472208059277854851UL) + ((uint64_t)op[4] * 15697300765962808534UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 15697300765962808534UL) + ((uint64_t)op[1] * 15422499593974599487UL) + ((uint64_t)op[2] * 17938432657587224280UL) + ((((uint64_t)op[3] * 3057272475926735445UL) + ((uint64_t)op[4] * 12472208059277854851UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 12472208059277854851UL) + ((uint64_t)op[1] * 15697300765962808534UL) + ((uint64_t)op[2] * 15422499593974599487UL) + ((uint64_t)op[3] * 17938432657587224280UL) + ((uint64_t)op[4] * 15286362379633677225UL);
	tmp_q[4] = ((uint64_t)op[0] * 3057272475926735445UL) + ((uint64_t)op[1] * 12472208059277854851UL) + ((uint64_t)op[2] * 15697300765962808534UL) + ((uint64_t)op[3] * 15422499593974599487UL) + ((uint64_t)op[4] * 17938432657587224280UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 149396806903L) + ((-((int128)tmp_q[1] * 48484076389L) + ((int128)tmp_q[2] * 88836261515L) + ((int128)tmp_q[3] * 33847226090L) - ((int128)tmp_q[4] * 240411783218L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 240411783218L) - ((int128)tmp_q[1] * 149396806903L) + ((-((int128)tmp_q[2] * 48484076389L) + ((int128)tmp_q[3] * 88836261515L) + ((int128)tmp_q[4] * 33847226090L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 33847226090L) - ((int128)tmp_q[1] * 240411783218L) - ((int128)tmp_q[2] * 149396806903L) + ((-((int128)tmp_q[3] * 48484076389L) + ((int128)tmp_q[4] * 88836261515L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 88836261515L) + ((int128)tmp_q[1] * 33847226090L) - ((int128)tmp_q[2] * 240411783218L) - ((int128)tmp_q[3] * 149396806903L) - ((int128)tmp_q[4] * 242420381945L);
	tmp_zero[4] = -((int128)tmp_q[0] * 48484076389L) + ((int128)tmp_q[1] * 88836261515L) + ((int128)tmp_q[2] * 33847226090L) - ((int128)tmp_q[3] * 240411783218L) - ((int128)tmp_q[4] * 149396806903L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

