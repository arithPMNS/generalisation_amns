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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 12923174411101449726UL) + ((((uint64_t)op[1] * 6779119785774550500UL) + ((uint64_t)op[2] * 4049263870938315725UL) + ((uint64_t)op[3] * 11909685136488313388UL) + ((uint64_t)op[4] * 9142678203312377969UL) + ((uint64_t)op[5] * 11651497520573810251UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 11651497520573810251UL) + ((uint64_t)op[1] * 12923174411101449726UL) + ((((uint64_t)op[2] * 6779119785774550500UL) + ((uint64_t)op[3] * 4049263870938315725UL) + ((uint64_t)op[4] * 11909685136488313388UL) + ((uint64_t)op[5] * 9142678203312377969UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 9142678203312377969UL) + ((uint64_t)op[1] * 11651497520573810251UL) + ((uint64_t)op[2] * 12923174411101449726UL) + ((((uint64_t)op[3] * 6779119785774550500UL) + ((uint64_t)op[4] * 4049263870938315725UL) + ((uint64_t)op[5] * 11909685136488313388UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 11909685136488313388UL) + ((uint64_t)op[1] * 9142678203312377969UL) + ((uint64_t)op[2] * 11651497520573810251UL) + ((uint64_t)op[3] * 12923174411101449726UL) + ((((uint64_t)op[4] * 6779119785774550500UL) + ((uint64_t)op[5] * 4049263870938315725UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 4049263870938315725UL) + ((uint64_t)op[1] * 11909685136488313388UL) + ((uint64_t)op[2] * 9142678203312377969UL) + ((uint64_t)op[3] * 11651497520573810251UL) + ((uint64_t)op[4] * 12923174411101449726UL) + ((uint64_t)op[5] * 16556128790095451732UL);
	tmp_q[5] = ((uint64_t)op[0] * 6779119785774550500UL) + ((uint64_t)op[1] * 4049263870938315725UL) + ((uint64_t)op[2] * 11909685136488313388UL) + ((uint64_t)op[3] * 9142678203312377969UL) + ((uint64_t)op[4] * 11651497520573810251UL) + ((uint64_t)op[5] * 12923174411101449726UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 114114755L) - ((((int128)tmp_q[1] * 1061862834L) + ((int128)tmp_q[2] * 826926767L) + ((int128)tmp_q[3] * 184828891L) + ((int128)tmp_q[4] * 2093728500L) + ((int128)tmp_q[5] * 1136337904L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 1136337904L) - ((int128)tmp_q[1] * 114114755L) - ((((int128)tmp_q[2] * 1061862834L) + ((int128)tmp_q[3] * 826926767L) + ((int128)tmp_q[4] * 184828891L) + ((int128)tmp_q[5] * 2093728500L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 2093728500L) + ((int128)tmp_q[1] * 1136337904L) - ((int128)tmp_q[2] * 114114755L) - ((((int128)tmp_q[3] * 1061862834L) + ((int128)tmp_q[4] * 826926767L) + ((int128)tmp_q[5] * 184828891L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 184828891L) + ((int128)tmp_q[1] * 2093728500L) + ((int128)tmp_q[2] * 1136337904L) - ((int128)tmp_q[3] * 114114755L) - ((((int128)tmp_q[4] * 1061862834L) + ((int128)tmp_q[5] * 826926767L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 826926767L) + ((int128)tmp_q[1] * 184828891L) + ((int128)tmp_q[2] * 2093728500L) + ((int128)tmp_q[3] * 1136337904L) - ((int128)tmp_q[4] * 114114755L) - ((int128)tmp_q[5] * 3185588502L);
	tmp_zero[5] = ((int128)tmp_q[0] * 1061862834L) + ((int128)tmp_q[1] * 826926767L) + ((int128)tmp_q[2] * 184828891L) + ((int128)tmp_q[3] * 2093728500L) + ((int128)tmp_q[4] * 1136337904L) - ((int128)tmp_q[5] * 114114755L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

