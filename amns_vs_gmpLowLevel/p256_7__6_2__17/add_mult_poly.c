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
	tmp_q[0] = ((uint64_t)op[0] * 5535990506466656481UL) + ((((uint64_t)op[1] * 3005017451729258402UL) + ((uint64_t)op[2] * 4947971490383875886UL) + ((uint64_t)op[3] * 11920798710963246841UL) + ((uint64_t)op[4] * 16147546380131977318UL) + ((uint64_t)op[5] * 11422361791120664570UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 11422361791120664570UL) + ((uint64_t)op[1] * 5535990506466656481UL) + ((((uint64_t)op[2] * 3005017451729258402UL) + ((uint64_t)op[3] * 4947971490383875886UL) + ((uint64_t)op[4] * 11920798710963246841UL) + ((uint64_t)op[5] * 16147546380131977318UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 16147546380131977318UL) + ((uint64_t)op[1] * 11422361791120664570UL) + ((uint64_t)op[2] * 5535990506466656481UL) + ((((uint64_t)op[3] * 3005017451729258402UL) + ((uint64_t)op[4] * 4947971490383875886UL) + ((uint64_t)op[5] * 11920798710963246841UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 11920798710963246841UL) + ((uint64_t)op[1] * 16147546380131977318UL) + ((uint64_t)op[2] * 11422361791120664570UL) + ((uint64_t)op[3] * 5535990506466656481UL) + ((((uint64_t)op[4] * 3005017451729258402UL) + ((uint64_t)op[5] * 4947971490383875886UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 4947971490383875886UL) + ((uint64_t)op[1] * 11920798710963246841UL) + ((uint64_t)op[2] * 16147546380131977318UL) + ((uint64_t)op[3] * 11422361791120664570UL) + ((uint64_t)op[4] * 5535990506466656481UL) + ((uint64_t)op[5] * 6010034903458516804UL);
	tmp_q[5] = ((uint64_t)op[0] * 3005017451729258402UL) + ((uint64_t)op[1] * 4947971490383875886UL) + ((uint64_t)op[2] * 11920798710963246841UL) + ((uint64_t)op[3] * 16147546380131977318UL) + ((uint64_t)op[4] * 11422361791120664570UL) + ((uint64_t)op[5] * 5535990506466656481UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 3933396324049L) + ((-((int128)tmp_q[1] * 441119170386L) - ((int128)tmp_q[2] * 1983648169854L) - ((int128)tmp_q[3] * 1478366753465L) - ((int128)tmp_q[4] * 733858059490L) - ((int128)tmp_q[5] * 2513286432354L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 2513286432354L) + ((int128)tmp_q[1] * 3933396324049L) + ((-((int128)tmp_q[2] * 441119170386L) - ((int128)tmp_q[3] * 1983648169854L) - ((int128)tmp_q[4] * 1478366753465L) - ((int128)tmp_q[5] * 733858059490L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 733858059490L) - ((int128)tmp_q[1] * 2513286432354L) + ((int128)tmp_q[2] * 3933396324049L) + ((-((int128)tmp_q[3] * 441119170386L) - ((int128)tmp_q[4] * 1983648169854L) - ((int128)tmp_q[5] * 1478366753465L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 1478366753465L) - ((int128)tmp_q[1] * 733858059490L) - ((int128)tmp_q[2] * 2513286432354L) + ((int128)tmp_q[3] * 3933396324049L) + ((-((int128)tmp_q[4] * 441119170386L) - ((int128)tmp_q[5] * 1983648169854L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 1983648169854L) - ((int128)tmp_q[1] * 1478366753465L) - ((int128)tmp_q[2] * 733858059490L) - ((int128)tmp_q[3] * 2513286432354L) + ((int128)tmp_q[4] * 3933396324049L) - ((int128)tmp_q[5] * 882238340772L);
	tmp_zero[5] = -((int128)tmp_q[0] * 441119170386L) - ((int128)tmp_q[1] * 1983648169854L) - ((int128)tmp_q[2] * 1478366753465L) - ((int128)tmp_q[3] * 733858059490L) - ((int128)tmp_q[4] * 2513286432354L) + ((int128)tmp_q[5] * 3933396324049L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

