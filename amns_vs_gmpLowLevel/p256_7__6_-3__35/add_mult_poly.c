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
	tmp_q[0] = ((uint64_t)op[0] * 7694687945627488685UL) + ((((uint64_t)op[1] * 6576743591567670640UL) + ((uint64_t)op[2] * 1621745733546245521UL) + ((uint64_t)op[3] * 10873980292445264822UL) + ((uint64_t)op[4] * 18202941241653918330UL) + ((uint64_t)op[5] * 10859035981776509861UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 10859035981776509861UL) + ((uint64_t)op[1] * 7694687945627488685UL) + ((((uint64_t)op[2] * 6576743591567670640UL) + ((uint64_t)op[3] * 1621745733546245521UL) + ((uint64_t)op[4] * 10873980292445264822UL) + ((uint64_t)op[5] * 18202941241653918330UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 18202941241653918330UL) + ((uint64_t)op[1] * 10859035981776509861UL) + ((uint64_t)op[2] * 7694687945627488685UL) + ((((uint64_t)op[3] * 6576743591567670640UL) + ((uint64_t)op[4] * 1621745733546245521UL) + ((uint64_t)op[5] * 10873980292445264822UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 10873980292445264822UL) + ((uint64_t)op[1] * 18202941241653918330UL) + ((uint64_t)op[2] * 10859035981776509861UL) + ((uint64_t)op[3] * 7694687945627488685UL) + ((((uint64_t)op[4] * 6576743591567670640UL) + ((uint64_t)op[5] * 1621745733546245521UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 1621745733546245521UL) + ((uint64_t)op[1] * 10873980292445264822UL) + ((uint64_t)op[2] * 18202941241653918330UL) + ((uint64_t)op[3] * 10859035981776509861UL) + ((uint64_t)op[4] * 7694687945627488685UL) + ((uint64_t)op[5] * 17163257372716091312UL);
	tmp_q[5] = ((uint64_t)op[0] * 6576743591567670640UL) + ((uint64_t)op[1] * 1621745733546245521UL) + ((uint64_t)op[2] * 10873980292445264822UL) + ((uint64_t)op[3] * 18202941241653918330UL) + ((uint64_t)op[4] * 10859035981776509861UL) + ((uint64_t)op[5] * 7694687945627488685UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 905234880197L) - ((-((int128)tmp_q[1] * 2054045311938L) - ((int128)tmp_q[2] * 220802796095L) + ((int128)tmp_q[3] * 1232196310098L) + ((int128)tmp_q[4] * 962507224628L) - ((int128)tmp_q[5] * 5725812058935L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 5725812058935L) - ((int128)tmp_q[1] * 905234880197L) - ((-((int128)tmp_q[2] * 2054045311938L) - ((int128)tmp_q[3] * 220802796095L) + ((int128)tmp_q[4] * 1232196310098L) + ((int128)tmp_q[5] * 962507224628L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 962507224628L) - ((int128)tmp_q[1] * 5725812058935L) - ((int128)tmp_q[2] * 905234880197L) - ((-((int128)tmp_q[3] * 2054045311938L) - ((int128)tmp_q[4] * 220802796095L) + ((int128)tmp_q[5] * 1232196310098L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 1232196310098L) + ((int128)tmp_q[1] * 962507224628L) - ((int128)tmp_q[2] * 5725812058935L) - ((int128)tmp_q[3] * 905234880197L) - ((-((int128)tmp_q[4] * 2054045311938L) - ((int128)tmp_q[5] * 220802796095L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 220802796095L) + ((int128)tmp_q[1] * 1232196310098L) + ((int128)tmp_q[2] * 962507224628L) - ((int128)tmp_q[3] * 5725812058935L) - ((int128)tmp_q[4] * 905234880197L) + ((int128)tmp_q[5] * 6162135935814L);
	tmp_zero[5] = -((int128)tmp_q[0] * 2054045311938L) - ((int128)tmp_q[1] * 220802796095L) + ((int128)tmp_q[2] * 1232196310098L) + ((int128)tmp_q[3] * 962507224628L) - ((int128)tmp_q[4] * 5725812058935L) - ((int128)tmp_q[5] * 905234880197L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

