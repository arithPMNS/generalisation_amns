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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 4);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) << 4);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8756201372580914963UL) + ((((uint64_t)op[1] * 18339909272683396649UL) + ((uint64_t)op[2] * 10227301331451651194UL) + ((uint64_t)op[3] * 6597304741699617838UL) + ((uint64_t)op[4] * 9572919216684215471UL) + ((uint64_t)op[5] * 12777676616383265288UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 12777676616383265288UL) + ((uint64_t)op[1] * 8756201372580914963UL) + ((((uint64_t)op[2] * 18339909272683396649UL) + ((uint64_t)op[3] * 10227301331451651194UL) + ((uint64_t)op[4] * 6597304741699617838UL) + ((uint64_t)op[5] * 9572919216684215471UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 9572919216684215471UL) + ((uint64_t)op[1] * 12777676616383265288UL) + ((uint64_t)op[2] * 8756201372580914963UL) + ((((uint64_t)op[3] * 18339909272683396649UL) + ((uint64_t)op[4] * 10227301331451651194UL) + ((uint64_t)op[5] * 6597304741699617838UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 6597304741699617838UL) + ((uint64_t)op[1] * 9572919216684215471UL) + ((uint64_t)op[2] * 12777676616383265288UL) + ((uint64_t)op[3] * 8756201372580914963UL) + ((((uint64_t)op[4] * 18339909272683396649UL) + ((uint64_t)op[5] * 10227301331451651194UL)) * 8);
	tmp_q[4] = ((uint64_t)op[0] * 10227301331451651194UL) + ((uint64_t)op[1] * 6597304741699617838UL) + ((uint64_t)op[2] * 9572919216684215471UL) + ((uint64_t)op[3] * 12777676616383265288UL) + ((uint64_t)op[4] * 8756201372580914963UL) + ((uint64_t)op[5] * 17592065665500311880UL);
	tmp_q[5] = ((uint64_t)op[0] * 18339909272683396649UL) + ((uint64_t)op[1] * 10227301331451651194UL) + ((uint64_t)op[2] * 6597304741699617838UL) + ((uint64_t)op[3] * 9572919216684215471UL) + ((uint64_t)op[4] * 12777676616383265288UL) + ((uint64_t)op[5] * 8756201372580914963UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 3687154682147L) + ((((int128)tmp_q[1] * 53244085578013L) + ((int128)tmp_q[2] * 47451342928743L) + ((int128)tmp_q[3] * 75444878711014L) - ((int128)tmp_q[4] * 16351362100081L) - ((int128)tmp_q[5] * 51033460750008L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 51033460750008L) - ((int128)tmp_q[1] * 3687154682147L) + ((((int128)tmp_q[2] * 53244085578013L) + ((int128)tmp_q[3] * 47451342928743L) + ((int128)tmp_q[4] * 75444878711014L) - ((int128)tmp_q[5] * 16351362100081L)) * 8);
	tmp_zero[2] = -((int128)tmp_q[0] * 16351362100081L) - ((int128)tmp_q[1] * 51033460750008L) - ((int128)tmp_q[2] * 3687154682147L) + ((((int128)tmp_q[3] * 53244085578013L) + ((int128)tmp_q[4] * 47451342928743L) + ((int128)tmp_q[5] * 75444878711014L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 75444878711014L) - ((int128)tmp_q[1] * 16351362100081L) - ((int128)tmp_q[2] * 51033460750008L) - ((int128)tmp_q[3] * 3687154682147L) + ((((int128)tmp_q[4] * 53244085578013L) + ((int128)tmp_q[5] * 47451342928743L)) * 8);
	tmp_zero[4] = ((int128)tmp_q[0] * 47451342928743L) + ((int128)tmp_q[1] * 75444878711014L) - ((int128)tmp_q[2] * 16351362100081L) - ((int128)tmp_q[3] * 51033460750008L) - ((int128)tmp_q[4] * 3687154682147L) + ((int128)tmp_q[5] * 425952684624104L);
	tmp_zero[5] = ((int128)tmp_q[0] * 53244085578013L) + ((int128)tmp_q[1] * 47451342928743L) + ((int128)tmp_q[2] * 75444878711014L) - ((int128)tmp_q[3] * 16351362100081L) - ((int128)tmp_q[4] * 51033460750008L) - ((int128)tmp_q[5] * 3687154682147L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

