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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9725043111099715293UL) + ((((uint64_t)op[1] * 6805857254792713640UL) + ((uint64_t)op[2] * 1649493958188036816UL) + ((uint64_t)op[3] * 1710334226614777016UL) + ((uint64_t)op[4] * 10659527242168318536UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 10659527242168318536UL) + ((uint64_t)op[1] * 9725043111099715293UL) + ((((uint64_t)op[2] * 6805857254792713640UL) + ((uint64_t)op[3] * 1649493958188036816UL) + ((uint64_t)op[4] * 1710334226614777016UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 1710334226614777016UL) + ((uint64_t)op[1] * 10659527242168318536UL) + ((uint64_t)op[2] * 9725043111099715293UL) + ((((uint64_t)op[3] * 6805857254792713640UL) + ((uint64_t)op[4] * 1649493958188036816UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 1649493958188036816UL) + ((uint64_t)op[1] * 1710334226614777016UL) + ((uint64_t)op[2] * 10659527242168318536UL) + ((uint64_t)op[3] * 9725043111099715293UL) + ((uint64_t)op[4] * 9670059128248248672UL);
	tmp_q[4] = ((uint64_t)op[0] * 6805857254792713640UL) + ((uint64_t)op[1] * 1649493958188036816UL) + ((uint64_t)op[2] * 1710334226614777016UL) + ((uint64_t)op[3] * 10659527242168318536UL) + ((uint64_t)op[4] * 9725043111099715293UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 279817651778421L) - ((-((int128)tmp_q[1] * 594089823650264L) - ((int128)tmp_q[2] * 404351826169392L) - ((int128)tmp_q[3] * 1727832042786888L) + ((int128)tmp_q[4] * 911234734642696L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 911234734642696L) - ((int128)tmp_q[1] * 279817651778421L) - ((-((int128)tmp_q[2] * 594089823650264L) - ((int128)tmp_q[3] * 404351826169392L) - ((int128)tmp_q[4] * 1727832042786888L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 1727832042786888L) + ((int128)tmp_q[1] * 911234734642696L) - ((int128)tmp_q[2] * 279817651778421L) - ((-((int128)tmp_q[3] * 594089823650264L) - ((int128)tmp_q[4] * 404351826169392L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 404351826169392L) - ((int128)tmp_q[1] * 1727832042786888L) + ((int128)tmp_q[2] * 911234734642696L) - ((int128)tmp_q[3] * 279817651778421L) + ((int128)tmp_q[4] * 2376359294601056L);
	tmp_zero[4] = -((int128)tmp_q[0] * 594089823650264L) - ((int128)tmp_q[1] * 404351826169392L) - ((int128)tmp_q[2] * 1727832042786888L) + ((int128)tmp_q[3] * 911234734642696L) - ((int128)tmp_q[4] * 279817651778421L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

