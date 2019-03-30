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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 6016805552455520194UL) + ((((uint64_t)op[1] * 370568231674460387UL) + ((uint64_t)op[2] * 1997124399621172359UL) + ((uint64_t)op[3] * 7928607506248410045UL) + ((uint64_t)op[4] * 5399591082696582497UL) + ((uint64_t)op[5] * 17747767300811649083UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 17747767300811649083UL) + ((uint64_t)op[1] * 6016805552455520194UL) + ((((uint64_t)op[2] * 370568231674460387UL) + ((uint64_t)op[3] * 1997124399621172359UL) + ((uint64_t)op[4] * 7928607506248410045UL) + ((uint64_t)op[5] * 5399591082696582497UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 5399591082696582497UL) + ((uint64_t)op[1] * 17747767300811649083UL) + ((uint64_t)op[2] * 6016805552455520194UL) + ((((uint64_t)op[3] * 370568231674460387UL) + ((uint64_t)op[4] * 1997124399621172359UL) + ((uint64_t)op[5] * 7928607506248410045UL)) * 18446744073709551609);
	tmp_q[3] = ((uint64_t)op[0] * 7928607506248410045UL) + ((uint64_t)op[1] * 5399591082696582497UL) + ((uint64_t)op[2] * 17747767300811649083UL) + ((uint64_t)op[3] * 6016805552455520194UL) + ((((uint64_t)op[4] * 370568231674460387UL) + ((uint64_t)op[5] * 1997124399621172359UL)) * 18446744073709551609);
	tmp_q[4] = ((uint64_t)op[0] * 1997124399621172359UL) + ((uint64_t)op[1] * 7928607506248410045UL) + ((uint64_t)op[2] * 5399591082696582497UL) + ((uint64_t)op[3] * 17747767300811649083UL) + ((uint64_t)op[4] * 6016805552455520194UL) + ((uint64_t)op[5] * 15852766451988328907UL);
	tmp_q[5] = ((uint64_t)op[0] * 370568231674460387UL) + ((uint64_t)op[1] * 1997124399621172359UL) + ((uint64_t)op[2] * 7928607506248410045UL) + ((uint64_t)op[3] * 5399591082696582497UL) + ((uint64_t)op[4] * 17747767300811649083UL) + ((uint64_t)op[5] * 6016805552455520194UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 2660905878L) - ((((int128)tmp_q[1] * 9721053L) + ((int128)tmp_q[2] * 3084607525L) + ((int128)tmp_q[3] * 233614967L) + ((int128)tmp_q[4] * 676953203L) + ((int128)tmp_q[5] * 2295179429L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 2295179429L) - ((int128)tmp_q[1] * 2660905878L) - ((((int128)tmp_q[2] * 9721053L) + ((int128)tmp_q[3] * 3084607525L) + ((int128)tmp_q[4] * 233614967L) + ((int128)tmp_q[5] * 676953203L)) * 7);
	tmp_zero[2] = ((int128)tmp_q[0] * 676953203L) + ((int128)tmp_q[1] * 2295179429L) - ((int128)tmp_q[2] * 2660905878L) - ((((int128)tmp_q[3] * 9721053L) + ((int128)tmp_q[4] * 3084607525L) + ((int128)tmp_q[5] * 233614967L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 233614967L) + ((int128)tmp_q[1] * 676953203L) + ((int128)tmp_q[2] * 2295179429L) - ((int128)tmp_q[3] * 2660905878L) - ((((int128)tmp_q[4] * 9721053L) + ((int128)tmp_q[5] * 3084607525L)) * 7);
	tmp_zero[4] = ((int128)tmp_q[0] * 3084607525L) + ((int128)tmp_q[1] * 233614967L) + ((int128)tmp_q[2] * 676953203L) + ((int128)tmp_q[3] * 2295179429L) - ((int128)tmp_q[4] * 2660905878L) - ((int128)tmp_q[5] * 68047371L);
	tmp_zero[5] = ((int128)tmp_q[0] * 9721053L) + ((int128)tmp_q[1] * 3084607525L) + ((int128)tmp_q[2] * 233614967L) + ((int128)tmp_q[3] * 676953203L) + ((int128)tmp_q[4] * 2295179429L) - ((int128)tmp_q[5] * 2660905878L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

