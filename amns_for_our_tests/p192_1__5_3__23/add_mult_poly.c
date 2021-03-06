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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9649612614571870073UL) + ((((uint64_t)op[1] * 13533496615038886284UL) + ((uint64_t)op[2] * 6230560592922411660UL) + ((uint64_t)op[3] * 3429885167738783203UL) + ((uint64_t)op[4] * 12142041143975485207UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 12142041143975485207UL) + ((uint64_t)op[1] * 9649612614571870073UL) + ((((uint64_t)op[2] * 13533496615038886284UL) + ((uint64_t)op[3] * 6230560592922411660UL) + ((uint64_t)op[4] * 3429885167738783203UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 3429885167738783203UL) + ((uint64_t)op[1] * 12142041143975485207UL) + ((uint64_t)op[2] * 9649612614571870073UL) + ((((uint64_t)op[3] * 13533496615038886284UL) + ((uint64_t)op[4] * 6230560592922411660UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 6230560592922411660UL) + ((uint64_t)op[1] * 3429885167738783203UL) + ((uint64_t)op[2] * 12142041143975485207UL) + ((uint64_t)op[3] * 9649612614571870073UL) + ((uint64_t)op[4] * 3707001697697555620UL);
	tmp_q[4] = ((uint64_t)op[0] * 13533496615038886284UL) + ((uint64_t)op[1] * 6230560592922411660UL) + ((uint64_t)op[2] * 3429885167738783203UL) + ((uint64_t)op[3] * 12142041143975485207UL) + ((uint64_t)op[4] * 9649612614571870073UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 45180157208L) + ((-((int128)tmp_q[1] * 30389362799L) + ((int128)tmp_q[2] * 176267378346L) + ((int128)tmp_q[3] * 71600914735L) - ((int128)tmp_q[4] * 34979488253L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 34979488253L) - ((int128)tmp_q[1] * 45180157208L) + ((-((int128)tmp_q[2] * 30389362799L) + ((int128)tmp_q[3] * 176267378346L) + ((int128)tmp_q[4] * 71600914735L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 71600914735L) - ((int128)tmp_q[1] * 34979488253L) - ((int128)tmp_q[2] * 45180157208L) + ((-((int128)tmp_q[3] * 30389362799L) + ((int128)tmp_q[4] * 176267378346L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 176267378346L) + ((int128)tmp_q[1] * 71600914735L) - ((int128)tmp_q[2] * 34979488253L) - ((int128)tmp_q[3] * 45180157208L) - ((int128)tmp_q[4] * 91168088397L);
	tmp_zero[4] = -((int128)tmp_q[0] * 30389362799L) + ((int128)tmp_q[1] * 176267378346L) + ((int128)tmp_q[2] * 71600914735L) - ((int128)tmp_q[3] * 34979488253L) - ((int128)tmp_q[4] * 45180157208L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

