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

//~ assumes 'scalar' and/or coeffs of 'op' small enough to avoid an overflow.
void double_poly_coeffs(int64_t *rop, int64_t *op){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = op[j] << 1;
}

//~ assumes 'nb_pos' and/or coeffs of 'op' small enough to avoid an overflow.
void lshift_poly_coeffs(int64_t *rop, int64_t *op, int nb_pos){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = op[j] << nb_pos;
}

//~ Computes pa(X)*pb(X) mod(E)
void mult_mod_poly(int64_t *rop, int64_t *pa, int64_t *pb){

	int128 tmp_prod_result[NB_COEFF];
	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(E)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];
	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 5187472560635430935UL) + ((((uint64_t)op[1] * 2177329607927910459UL) + ((uint64_t)op[2] * 7147122787381547036UL) + ((uint64_t)op[3] * 1358633054727937561UL) + ((uint64_t)op[4] * 4072120099921067645UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 4072120099921067645UL) + ((uint64_t)op[1] * 5187472560635430935UL) + ((((uint64_t)op[2] * 2177329607927910459UL) + ((uint64_t)op[3] * 7147122787381547036UL) + ((uint64_t)op[4] * 1358633054727937561UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 1358633054727937561UL) + ((uint64_t)op[1] * 4072120099921067645UL) + ((uint64_t)op[2] * 5187472560635430935UL) + ((((uint64_t)op[3] * 2177329607927910459UL) + ((uint64_t)op[4] * 7147122787381547036UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 7147122787381547036UL) + ((uint64_t)op[1] * 1358633054727937561UL) + ((uint64_t)op[2] * 4072120099921067645UL) + ((uint64_t)op[3] * 5187472560635430935UL) + ((uint64_t)op[4] * 4354659215855820918UL);
	tmp_q[4] = ((uint64_t)op[0] * 2177329607927910459UL) + ((uint64_t)op[1] * 7147122787381547036UL) + ((uint64_t)op[2] * 1358633054727937561UL) + ((uint64_t)op[3] * 4072120099921067645UL) + ((uint64_t)op[4] * 5187472560635430935UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 474667993690623L) + ((((int128)tmp_q[1] * 1474731667115964L) + ((int128)tmp_q[2] * 148617409519721L) + ((int128)tmp_q[3] * 38605154947216L) + ((int128)tmp_q[4] * 1402623091101957L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 1402623091101957L) - ((int128)tmp_q[1] * 474667993690623L) + ((((int128)tmp_q[2] * 1474731667115964L) + ((int128)tmp_q[3] * 148617409519721L) + ((int128)tmp_q[4] * 38605154947216L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 38605154947216L) + ((int128)tmp_q[1] * 1402623091101957L) - ((int128)tmp_q[2] * 474667993690623L) + ((((int128)tmp_q[3] * 1474731667115964L) + ((int128)tmp_q[4] * 148617409519721L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 148617409519721L) + ((int128)tmp_q[1] * 38605154947216L) + ((int128)tmp_q[2] * 1402623091101957L) - ((int128)tmp_q[3] * 474667993690623L) + ((int128)tmp_q[4] * 2949463334231928L);
	tmp_zero[4] = ((int128)tmp_q[0] * 1474731667115964L) + ((int128)tmp_q[1] * 148617409519721L) + ((int128)tmp_q[2] * 38605154947216L) + ((int128)tmp_q[3] * 1402623091101957L) - ((int128)tmp_q[4] * 474667993690623L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

void exact_coeffs_reduction(int64_t *rop, int64_t *op){

	int i;
	int128 tmp[NB_COEFF];

	for(i=0; i<NB_COEFF; i++)
		tmp[i] = (int128) op[i];

	internal_reduction(rop, tmp);

	tmp[0] = (int128)rop[0] * poly_P0[0] + (((int128)rop[1] * poly_P0[4] + (int128)rop[2] * poly_P0[3] + (int128)rop[3] * poly_P0[2] + (int128)rop[4] * poly_P0[1]) << 1);
	tmp[1] = (int128)rop[0] * poly_P0[1] + (int128)rop[1] * poly_P0[0] + (((int128)rop[2] * poly_P0[4] + (int128)rop[3] * poly_P0[3] + (int128)rop[4] * poly_P0[2]) << 1);
	tmp[2] = (int128)rop[0] * poly_P0[2] + (int128)rop[1] * poly_P0[1] + (int128)rop[2] * poly_P0[0] + (((int128)rop[3] * poly_P0[4] + (int128)rop[4] * poly_P0[3]) << 1);
	tmp[3] = (int128)rop[0] * poly_P0[3] + (int128)rop[1] * poly_P0[2] + (int128)rop[2] * poly_P0[1] + (int128)rop[3] * poly_P0[0] + (((int128)rop[4] * poly_P0[4]) << 1);
	tmp[4] = (int128)rop[0] * poly_P0[4] + (int128)rop[1] * poly_P0[3] + (int128)rop[2] * poly_P0[2] + (int128)rop[3] * poly_P0[1] + (int128)rop[4] * poly_P0[0];

	internal_reduction(rop, tmp);
}

