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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 16992171857263082351UL) + ((((uint64_t)op[1] * 9015769154231507125UL) + ((uint64_t)op[2] * 5050184328706037716UL) + ((uint64_t)op[3] * 18266563821396381327UL) + ((uint64_t)op[4] * 18262060407116180834UL)) * 6);
	tmp_q[1] = ((uint64_t)op[0] * 18262060407116180834UL) + ((uint64_t)op[1] * 16992171857263082351UL) + ((((uint64_t)op[2] * 9015769154231507125UL) + ((uint64_t)op[3] * 5050184328706037716UL) + ((uint64_t)op[4] * 18266563821396381327UL)) * 6);
	tmp_q[2] = ((uint64_t)op[0] * 18266563821396381327UL) + ((uint64_t)op[1] * 18262060407116180834UL) + ((uint64_t)op[2] * 16992171857263082351UL) + ((((uint64_t)op[3] * 9015769154231507125UL) + ((uint64_t)op[4] * 5050184328706037716UL)) * 6);
	tmp_q[3] = ((uint64_t)op[0] * 5050184328706037716UL) + ((uint64_t)op[1] * 18266563821396381327UL) + ((uint64_t)op[2] * 18262060407116180834UL) + ((uint64_t)op[3] * 16992171857263082351UL) + ((uint64_t)op[4] * 17201126777969939518UL);
	tmp_q[4] = ((uint64_t)op[0] * 9015769154231507125UL) + ((uint64_t)op[1] * 5050184328706037716UL) + ((uint64_t)op[2] * 18266563821396381327UL) + ((uint64_t)op[3] * 18262060407116180834UL) + ((uint64_t)op[4] * 16992171857263082351UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 14162789399619L) + ((((int128)tmp_q[1] * 837331460134L) - ((int128)tmp_q[2] * 23270885655922L) + ((int128)tmp_q[3] * 13138246039239L) + ((int128)tmp_q[4] * 17312134742264L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 17312134742264L) - ((int128)tmp_q[1] * 14162789399619L) + ((((int128)tmp_q[2] * 837331460134L) - ((int128)tmp_q[3] * 23270885655922L) + ((int128)tmp_q[4] * 13138246039239L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 13138246039239L) + ((int128)tmp_q[1] * 17312134742264L) - ((int128)tmp_q[2] * 14162789399619L) + ((((int128)tmp_q[3] * 837331460134L) - ((int128)tmp_q[4] * 23270885655922L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 23270885655922L) + ((int128)tmp_q[1] * 13138246039239L) + ((int128)tmp_q[2] * 17312134742264L) - ((int128)tmp_q[3] * 14162789399619L) + ((int128)tmp_q[4] * 5023988760804L);
	tmp_zero[4] = ((int128)tmp_q[0] * 837331460134L) - ((int128)tmp_q[1] * 23270885655922L) + ((int128)tmp_q[2] * 13138246039239L) + ((int128)tmp_q[3] * 17312134742264L) - ((int128)tmp_q[4] * 14162789399619L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

