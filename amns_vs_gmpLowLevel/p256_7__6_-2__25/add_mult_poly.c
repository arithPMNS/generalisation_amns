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
	tmp_q[0] = ((uint64_t)op[0] * 11730257001970320779UL) + ((((uint64_t)op[1] * 9027706382056148621UL) + ((uint64_t)op[2] * 6446316833248371443UL) + ((uint64_t)op[3] * 5564293098301397251UL) + ((uint64_t)op[4] * 7397601331250865344UL) + ((uint64_t)op[5] * 3228773011272606025UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 3228773011272606025UL) + ((uint64_t)op[1] * 11730257001970320779UL) + ((((uint64_t)op[2] * 9027706382056148621UL) + ((uint64_t)op[3] * 6446316833248371443UL) + ((uint64_t)op[4] * 5564293098301397251UL) + ((uint64_t)op[5] * 7397601331250865344UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 7397601331250865344UL) + ((uint64_t)op[1] * 3228773011272606025UL) + ((uint64_t)op[2] * 11730257001970320779UL) + ((((uint64_t)op[3] * 9027706382056148621UL) + ((uint64_t)op[4] * 6446316833248371443UL) + ((uint64_t)op[5] * 5564293098301397251UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 5564293098301397251UL) + ((uint64_t)op[1] * 7397601331250865344UL) + ((uint64_t)op[2] * 3228773011272606025UL) + ((uint64_t)op[3] * 11730257001970320779UL) + ((((uint64_t)op[4] * 9027706382056148621UL) + ((uint64_t)op[5] * 6446316833248371443UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 6446316833248371443UL) + ((uint64_t)op[1] * 5564293098301397251UL) + ((uint64_t)op[2] * 7397601331250865344UL) + ((uint64_t)op[3] * 3228773011272606025UL) + ((uint64_t)op[4] * 11730257001970320779UL) + ((uint64_t)op[5] * 391331309597254374UL);
	tmp_q[5] = ((uint64_t)op[0] * 9027706382056148621UL) + ((uint64_t)op[1] * 6446316833248371443UL) + ((uint64_t)op[2] * 5564293098301397251UL) + ((uint64_t)op[3] * 7397601331250865344UL) + ((uint64_t)op[4] * 3228773011272606025UL) + ((uint64_t)op[5] * 11730257001970320779UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 655361941023L) - ((-((int128)tmp_q[1] * 212902645429L) - ((int128)tmp_q[2] * 3578164725036L) - ((int128)tmp_q[3] * 2619096791750L) + ((int128)tmp_q[4] * 4853619446231L) - ((int128)tmp_q[5] * 1548574765031L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 1548574765031L) + ((int128)tmp_q[1] * 655361941023L) - ((-((int128)tmp_q[2] * 212902645429L) - ((int128)tmp_q[3] * 3578164725036L) - ((int128)tmp_q[4] * 2619096791750L) + ((int128)tmp_q[5] * 4853619446231L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 4853619446231L) - ((int128)tmp_q[1] * 1548574765031L) + ((int128)tmp_q[2] * 655361941023L) - ((-((int128)tmp_q[3] * 212902645429L) - ((int128)tmp_q[4] * 3578164725036L) - ((int128)tmp_q[5] * 2619096791750L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 2619096791750L) + ((int128)tmp_q[1] * 4853619446231L) - ((int128)tmp_q[2] * 1548574765031L) + ((int128)tmp_q[3] * 655361941023L) - ((-((int128)tmp_q[4] * 212902645429L) - ((int128)tmp_q[5] * 3578164725036L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 3578164725036L) - ((int128)tmp_q[1] * 2619096791750L) + ((int128)tmp_q[2] * 4853619446231L) - ((int128)tmp_q[3] * 1548574765031L) + ((int128)tmp_q[4] * 655361941023L) + ((int128)tmp_q[5] * 425805290858L);
	tmp_zero[5] = -((int128)tmp_q[0] * 212902645429L) - ((int128)tmp_q[1] * 3578164725036L) - ((int128)tmp_q[2] * 2619096791750L) + ((int128)tmp_q[3] * 4853619446231L) - ((int128)tmp_q[4] * 1548574765031L) + ((int128)tmp_q[5] * 655361941023L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

