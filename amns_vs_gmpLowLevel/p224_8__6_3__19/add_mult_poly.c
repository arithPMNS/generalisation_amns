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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 15154032348299423323UL) + ((((uint64_t)op[1] * 10871583149702586215UL) + ((uint64_t)op[2] * 15351441352570507189UL) + ((uint64_t)op[3] * 12699254344351552873UL) + ((uint64_t)op[4] * 11790641955545644144UL) + ((uint64_t)op[5] * 16282528343281499181UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 16282528343281499181UL) + ((uint64_t)op[1] * 15154032348299423323UL) + ((((uint64_t)op[2] * 10871583149702586215UL) + ((uint64_t)op[3] * 15351441352570507189UL) + ((uint64_t)op[4] * 12699254344351552873UL) + ((uint64_t)op[5] * 11790641955545644144UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 11790641955545644144UL) + ((uint64_t)op[1] * 16282528343281499181UL) + ((uint64_t)op[2] * 15154032348299423323UL) + ((((uint64_t)op[3] * 10871583149702586215UL) + ((uint64_t)op[4] * 15351441352570507189UL) + ((uint64_t)op[5] * 12699254344351552873UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 12699254344351552873UL) + ((uint64_t)op[1] * 11790641955545644144UL) + ((uint64_t)op[2] * 16282528343281499181UL) + ((uint64_t)op[3] * 15154032348299423323UL) + ((((uint64_t)op[4] * 10871583149702586215UL) + ((uint64_t)op[5] * 15351441352570507189UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 15351441352570507189UL) + ((uint64_t)op[1] * 12699254344351552873UL) + ((uint64_t)op[2] * 11790641955545644144UL) + ((uint64_t)op[3] * 16282528343281499181UL) + ((uint64_t)op[4] * 15154032348299423323UL) + ((uint64_t)op[5] * 14168005375398207029UL);
	tmp_q[5] = ((uint64_t)op[0] * 10871583149702586215UL) + ((uint64_t)op[1] * 15351441352570507189UL) + ((uint64_t)op[2] * 12699254344351552873UL) + ((uint64_t)op[3] * 11790641955545644144UL) + ((uint64_t)op[4] * 16282528343281499181UL) + ((uint64_t)op[5] * 15154032348299423323UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 11945463949L) + ((((int128)tmp_q[1] * 19126841159L) + ((int128)tmp_q[2] * 17952757412L) - ((int128)tmp_q[3] * 44512270717L) + ((int128)tmp_q[4] * 26619744611L) - ((int128)tmp_q[5] * 104776402653L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 104776402653L) - ((int128)tmp_q[1] * 11945463949L) + ((((int128)tmp_q[2] * 19126841159L) + ((int128)tmp_q[3] * 17952757412L) - ((int128)tmp_q[4] * 44512270717L) + ((int128)tmp_q[5] * 26619744611L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 26619744611L) - ((int128)tmp_q[1] * 104776402653L) - ((int128)tmp_q[2] * 11945463949L) + ((((int128)tmp_q[3] * 19126841159L) + ((int128)tmp_q[4] * 17952757412L) - ((int128)tmp_q[5] * 44512270717L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 44512270717L) + ((int128)tmp_q[1] * 26619744611L) - ((int128)tmp_q[2] * 104776402653L) - ((int128)tmp_q[3] * 11945463949L) + ((((int128)tmp_q[4] * 19126841159L) + ((int128)tmp_q[5] * 17952757412L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 17952757412L) - ((int128)tmp_q[1] * 44512270717L) + ((int128)tmp_q[2] * 26619744611L) - ((int128)tmp_q[3] * 104776402653L) - ((int128)tmp_q[4] * 11945463949L) + ((int128)tmp_q[5] * 57380523477L);
	tmp_zero[5] = ((int128)tmp_q[0] * 19126841159L) + ((int128)tmp_q[1] * 17952757412L) - ((int128)tmp_q[2] * 44512270717L) + ((int128)tmp_q[3] * 26619744611L) - ((int128)tmp_q[4] * 104776402653L) - ((int128)tmp_q[5] * 11945463949L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

