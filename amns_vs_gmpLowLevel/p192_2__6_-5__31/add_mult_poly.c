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
	tmp_q[0] = ((uint64_t)op[0] * 13509358164794378647UL) + ((((uint64_t)op[1] * 18191453500437988773UL) + ((uint64_t)op[2] * 10131509960313030048UL) + ((uint64_t)op[3] * 14303843615732932567UL) + ((uint64_t)op[4] * 14954915072424562808UL) + ((uint64_t)op[5] * 6789643007526163606UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 6789643007526163606UL) + ((uint64_t)op[1] * 13509358164794378647UL) + ((((uint64_t)op[2] * 18191453500437988773UL) + ((uint64_t)op[3] * 10131509960313030048UL) + ((uint64_t)op[4] * 14303843615732932567UL) + ((uint64_t)op[5] * 14954915072424562808UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 14954915072424562808UL) + ((uint64_t)op[1] * 6789643007526163606UL) + ((uint64_t)op[2] * 13509358164794378647UL) + ((((uint64_t)op[3] * 18191453500437988773UL) + ((uint64_t)op[4] * 10131509960313030048UL) + ((uint64_t)op[5] * 14303843615732932567UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 14303843615732932567UL) + ((uint64_t)op[1] * 14954915072424562808UL) + ((uint64_t)op[2] * 6789643007526163606UL) + ((uint64_t)op[3] * 13509358164794378647UL) + ((((uint64_t)op[4] * 18191453500437988773UL) + ((uint64_t)op[5] * 10131509960313030048UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 10131509960313030048UL) + ((uint64_t)op[1] * 14303843615732932567UL) + ((uint64_t)op[2] * 14954915072424562808UL) + ((uint64_t)op[3] * 6789643007526163606UL) + ((uint64_t)op[4] * 13509358164794378647UL) + ((uint64_t)op[5] * 1276452866357814215UL);
	tmp_q[5] = ((uint64_t)op[0] * 18191453500437988773UL) + ((uint64_t)op[1] * 10131509960313030048UL) + ((uint64_t)op[2] * 14303843615732932567UL) + ((uint64_t)op[3] * 14954915072424562808UL) + ((uint64_t)op[4] * 6789643007526163606UL) + ((uint64_t)op[5] * 13509358164794378647UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 2170975586L) - ((-((int128)tmp_q[1] * 765884585L) + ((int128)tmp_q[2] * 3049338610L) + ((int128)tmp_q[3] * 1461675386L) + ((int128)tmp_q[4] * 2829073995L) + ((int128)tmp_q[5] * 201333771L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 201333771L) - ((int128)tmp_q[1] * 2170975586L) - ((-((int128)tmp_q[2] * 765884585L) + ((int128)tmp_q[3] * 3049338610L) + ((int128)tmp_q[4] * 1461675386L) + ((int128)tmp_q[5] * 2829073995L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 2829073995L) + ((int128)tmp_q[1] * 201333771L) - ((int128)tmp_q[2] * 2170975586L) - ((-((int128)tmp_q[3] * 765884585L) + ((int128)tmp_q[4] * 3049338610L) + ((int128)tmp_q[5] * 1461675386L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 1461675386L) + ((int128)tmp_q[1] * 2829073995L) + ((int128)tmp_q[2] * 201333771L) - ((int128)tmp_q[3] * 2170975586L) - ((-((int128)tmp_q[4] * 765884585L) + ((int128)tmp_q[5] * 3049338610L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 3049338610L) + ((int128)tmp_q[1] * 1461675386L) + ((int128)tmp_q[2] * 2829073995L) + ((int128)tmp_q[3] * 201333771L) - ((int128)tmp_q[4] * 2170975586L) + ((int128)tmp_q[5] * 3829422925L);
	tmp_zero[5] = -((int128)tmp_q[0] * 765884585L) + ((int128)tmp_q[1] * 3049338610L) + ((int128)tmp_q[2] * 1461675386L) + ((int128)tmp_q[3] * 2829073995L) + ((int128)tmp_q[4] * 201333771L) - ((int128)tmp_q[5] * 2170975586L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

