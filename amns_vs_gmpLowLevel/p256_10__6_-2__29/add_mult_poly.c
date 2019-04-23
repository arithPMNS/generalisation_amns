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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 15987870344290775861UL) + ((((uint64_t)op[1] * 9877489577903502938UL) + ((uint64_t)op[2] * 5190425846985348347UL) + ((uint64_t)op[3] * 6519639995484981198UL) + ((uint64_t)op[4] * 1904513922317272457UL) + ((uint64_t)op[5] * 6963193297314109121UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 6963193297314109121UL) + ((uint64_t)op[1] * 15987870344290775861UL) + ((((uint64_t)op[2] * 9877489577903502938UL) + ((uint64_t)op[3] * 5190425846985348347UL) + ((uint64_t)op[4] * 6519639995484981198UL) + ((uint64_t)op[5] * 1904513922317272457UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 1904513922317272457UL) + ((uint64_t)op[1] * 6963193297314109121UL) + ((uint64_t)op[2] * 15987870344290775861UL) + ((((uint64_t)op[3] * 9877489577903502938UL) + ((uint64_t)op[4] * 5190425846985348347UL) + ((uint64_t)op[5] * 6519639995484981198UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 6519639995484981198UL) + ((uint64_t)op[1] * 1904513922317272457UL) + ((uint64_t)op[2] * 6963193297314109121UL) + ((uint64_t)op[3] * 15987870344290775861UL) + ((((uint64_t)op[4] * 9877489577903502938UL) + ((uint64_t)op[5] * 5190425846985348347UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 5190425846985348347UL) + ((uint64_t)op[1] * 6519639995484981198UL) + ((uint64_t)op[2] * 1904513922317272457UL) + ((uint64_t)op[3] * 6963193297314109121UL) + ((uint64_t)op[4] * 15987870344290775861UL) + ((uint64_t)op[5] * 17138508991612097356UL);
	tmp_q[5] = ((uint64_t)op[0] * 9877489577903502938UL) + ((uint64_t)op[1] * 5190425846985348347UL) + ((uint64_t)op[2] * 6519639995484981198UL) + ((uint64_t)op[3] * 1904513922317272457UL) + ((uint64_t)op[4] * 6963193297314109121UL) + ((uint64_t)op[5] * 15987870344290775861UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 6013862895201L) - ((-((int128)tmp_q[1] * 903711916838L) + ((int128)tmp_q[2] * 339359158106L) - ((int128)tmp_q[3] * 51888874755L) + ((int128)tmp_q[4] * 2307439660562L) - ((int128)tmp_q[5] * 1400551118041L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 1400551118041L) - ((int128)tmp_q[1] * 6013862895201L) - ((-((int128)tmp_q[2] * 903711916838L) + ((int128)tmp_q[3] * 339359158106L) - ((int128)tmp_q[4] * 51888874755L) + ((int128)tmp_q[5] * 2307439660562L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 2307439660562L) - ((int128)tmp_q[1] * 1400551118041L) - ((int128)tmp_q[2] * 6013862895201L) - ((-((int128)tmp_q[3] * 903711916838L) + ((int128)tmp_q[4] * 339359158106L) - ((int128)tmp_q[5] * 51888874755L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 51888874755L) + ((int128)tmp_q[1] * 2307439660562L) - ((int128)tmp_q[2] * 1400551118041L) - ((int128)tmp_q[3] * 6013862895201L) - ((-((int128)tmp_q[4] * 903711916838L) + ((int128)tmp_q[5] * 339359158106L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 339359158106L) - ((int128)tmp_q[1] * 51888874755L) + ((int128)tmp_q[2] * 2307439660562L) - ((int128)tmp_q[3] * 1400551118041L) - ((int128)tmp_q[4] * 6013862895201L) + ((int128)tmp_q[5] * 1807423833676L);
	tmp_zero[5] = -((int128)tmp_q[0] * 903711916838L) + ((int128)tmp_q[1] * 339359158106L) - ((int128)tmp_q[2] * 51888874755L) + ((int128)tmp_q[3] * 2307439660562L) - ((int128)tmp_q[4] * 1400551118041L) - ((int128)tmp_q[5] * 6013862895201L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

