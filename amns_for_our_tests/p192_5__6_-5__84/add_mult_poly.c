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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 1130872062525751682UL) + ((((uint64_t)op[1] * 18308846456164130666UL) + ((uint64_t)op[2] * 4470426268912178369UL) + ((uint64_t)op[3] * 13636862696412875622UL) + ((uint64_t)op[4] * 1564642001959410790UL) + ((uint64_t)op[5] * 4531518877655207550UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 4531518877655207550UL) + ((uint64_t)op[1] * 1130872062525751682UL) + ((((uint64_t)op[2] * 18308846456164130666UL) + ((uint64_t)op[3] * 4470426268912178369UL) + ((uint64_t)op[4] * 13636862696412875622UL) + ((uint64_t)op[5] * 1564642001959410790UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 1564642001959410790UL) + ((uint64_t)op[1] * 4531518877655207550UL) + ((uint64_t)op[2] * 1130872062525751682UL) + ((((uint64_t)op[3] * 18308846456164130666UL) + ((uint64_t)op[4] * 4470426268912178369UL) + ((uint64_t)op[5] * 13636862696412875622UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 13636862696412875622UL) + ((uint64_t)op[1] * 1564642001959410790UL) + ((uint64_t)op[2] * 4531518877655207550UL) + ((uint64_t)op[3] * 1130872062525751682UL) + ((((uint64_t)op[4] * 18308846456164130666UL) + ((uint64_t)op[5] * 4470426268912178369UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 4470426268912178369UL) + ((uint64_t)op[1] * 13636862696412875622UL) + ((uint64_t)op[2] * 1564642001959410790UL) + ((uint64_t)op[3] * 4531518877655207550UL) + ((uint64_t)op[4] * 1130872062525751682UL) + ((uint64_t)op[5] * 689488087727104750UL);
	tmp_q[5] = ((uint64_t)op[0] * 18308846456164130666UL) + ((uint64_t)op[1] * 4470426268912178369UL) + ((uint64_t)op[2] * 13636862696412875622UL) + ((uint64_t)op[3] * 1564642001959410790UL) + ((uint64_t)op[4] * 4531518877655207550UL) + ((uint64_t)op[5] * 1130872062525751682UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 125975678L) - ((-((int128)tmp_q[1] * 2333892954L) + ((int128)tmp_q[2] * 539899394L) - ((int128)tmp_q[3] * 227746490L) - ((int128)tmp_q[4] * 2361973447L) - ((int128)tmp_q[5] * 1291276854L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 1291276854L) - ((int128)tmp_q[1] * 125975678L) - ((-((int128)tmp_q[2] * 2333892954L) + ((int128)tmp_q[3] * 539899394L) - ((int128)tmp_q[4] * 227746490L) - ((int128)tmp_q[5] * 2361973447L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 2361973447L) - ((int128)tmp_q[1] * 1291276854L) - ((int128)tmp_q[2] * 125975678L) - ((-((int128)tmp_q[3] * 2333892954L) + ((int128)tmp_q[4] * 539899394L) - ((int128)tmp_q[5] * 227746490L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 227746490L) - ((int128)tmp_q[1] * 2361973447L) - ((int128)tmp_q[2] * 1291276854L) - ((int128)tmp_q[3] * 125975678L) - ((-((int128)tmp_q[4] * 2333892954L) + ((int128)tmp_q[5] * 539899394L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 539899394L) - ((int128)tmp_q[1] * 227746490L) - ((int128)tmp_q[2] * 2361973447L) - ((int128)tmp_q[3] * 1291276854L) - ((int128)tmp_q[4] * 125975678L) + ((int128)tmp_q[5] * 11669464770L);
	tmp_zero[5] = -((int128)tmp_q[0] * 2333892954L) + ((int128)tmp_q[1] * 539899394L) - ((int128)tmp_q[2] * 227746490L) - ((int128)tmp_q[3] * 2361973447L) - ((int128)tmp_q[4] * 1291276854L) - ((int128)tmp_q[5] * 125975678L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

