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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) * 7);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) * 14);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) * 7);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 85169761754306884UL) + ((((uint64_t)op[1] * 6274753140656878189UL) + ((uint64_t)op[2] * 957644756629625206UL) + ((uint64_t)op[3] * 17765993184838414525UL) + ((uint64_t)op[4] * 8380338273045324582UL) + ((uint64_t)op[5] * 852992413085808503UL) + ((uint64_t)op[6] * 18187849982775220396UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 18187849982775220396UL) + ((uint64_t)op[1] * 85169761754306884UL) + ((((uint64_t)op[2] * 6274753140656878189UL) + ((uint64_t)op[3] * 957644756629625206UL) + ((uint64_t)op[4] * 17765993184838414525UL) + ((uint64_t)op[5] * 8380338273045324582UL) + ((uint64_t)op[6] * 852992413085808503UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 852992413085808503UL) + ((uint64_t)op[1] * 18187849982775220396UL) + ((uint64_t)op[2] * 85169761754306884UL) + ((((uint64_t)op[3] * 6274753140656878189UL) + ((uint64_t)op[4] * 957644756629625206UL) + ((uint64_t)op[5] * 17765993184838414525UL) + ((uint64_t)op[6] * 8380338273045324582UL)) * 18446744073709551609);
	tmp_q[3] = ((uint64_t)op[0] * 8380338273045324582UL) + ((uint64_t)op[1] * 852992413085808503UL) + ((uint64_t)op[2] * 18187849982775220396UL) + ((uint64_t)op[3] * 85169761754306884UL) + ((((uint64_t)op[4] * 6274753140656878189UL) + ((uint64_t)op[5] * 957644756629625206UL) + ((uint64_t)op[6] * 17765993184838414525UL)) * 18446744073709551609);
	tmp_q[4] = ((uint64_t)op[0] * 17765993184838414525UL) + ((uint64_t)op[1] * 8380338273045324582UL) + ((uint64_t)op[2] * 852992413085808503UL) + ((uint64_t)op[3] * 18187849982775220396UL) + ((uint64_t)op[4] * 85169761754306884UL) + ((((uint64_t)op[5] * 6274753140656878189UL) + ((uint64_t)op[6] * 957644756629625206UL)) * 18446744073709551609);
	tmp_q[5] = ((uint64_t)op[0] * 957644756629625206UL) + ((uint64_t)op[1] * 17765993184838414525UL) + ((uint64_t)op[2] * 8380338273045324582UL) + ((uint64_t)op[3] * 852992413085808503UL) + ((uint64_t)op[4] * 18187849982775220396UL) + ((uint64_t)op[5] * 85169761754306884UL) + ((uint64_t)op[6] * 11416960236530507525UL);
	tmp_q[6] = ((uint64_t)op[0] * 6274753140656878189UL) + ((uint64_t)op[1] * 957644756629625206UL) + ((uint64_t)op[2] * 17765993184838414525UL) + ((uint64_t)op[3] * 8380338273045324582UL) + ((uint64_t)op[4] * 852992413085808503UL) + ((uint64_t)op[5] * 18187849982775220396UL) + ((uint64_t)op[6] * 85169761754306884UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 33279420926L) - ((-((int128)tmp_q[1] * 25196387348L) - ((int128)tmp_q[2] * 13546913283L) - ((int128)tmp_q[3] * 33146317837L) - ((int128)tmp_q[4] * 19425544569L) - ((int128)tmp_q[5] * 15118856955L) - ((int128)tmp_q[6] * 28411115059L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 28411115059L) + ((int128)tmp_q[1] * 33279420926L) - ((-((int128)tmp_q[2] * 25196387348L) - ((int128)tmp_q[3] * 13546913283L) - ((int128)tmp_q[4] * 33146317837L) - ((int128)tmp_q[5] * 19425544569L) - ((int128)tmp_q[6] * 15118856955L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 15118856955L) - ((int128)tmp_q[1] * 28411115059L) + ((int128)tmp_q[2] * 33279420926L) - ((-((int128)tmp_q[3] * 25196387348L) - ((int128)tmp_q[4] * 13546913283L) - ((int128)tmp_q[5] * 33146317837L) - ((int128)tmp_q[6] * 19425544569L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 19425544569L) - ((int128)tmp_q[1] * 15118856955L) - ((int128)tmp_q[2] * 28411115059L) + ((int128)tmp_q[3] * 33279420926L) - ((-((int128)tmp_q[4] * 25196387348L) - ((int128)tmp_q[5] * 13546913283L) - ((int128)tmp_q[6] * 33146317837L)) * 7);
	tmp_zero[4] = -((int128)tmp_q[0] * 33146317837L) - ((int128)tmp_q[1] * 19425544569L) - ((int128)tmp_q[2] * 15118856955L) - ((int128)tmp_q[3] * 28411115059L) + ((int128)tmp_q[4] * 33279420926L) - ((-((int128)tmp_q[5] * 25196387348L) - ((int128)tmp_q[6] * 13546913283L)) * 7);
	tmp_zero[5] = -((int128)tmp_q[0] * 13546913283L) - ((int128)tmp_q[1] * 33146317837L) - ((int128)tmp_q[2] * 19425544569L) - ((int128)tmp_q[3] * 15118856955L) - ((int128)tmp_q[4] * 28411115059L) + ((int128)tmp_q[5] * 33279420926L) + ((int128)tmp_q[6] * 176374711436L);
	tmp_zero[6] = -((int128)tmp_q[0] * 25196387348L) - ((int128)tmp_q[1] * 13546913283L) - ((int128)tmp_q[2] * 33146317837L) - ((int128)tmp_q[3] * 19425544569L) - ((int128)tmp_q[4] * 15118856955L) - ((int128)tmp_q[5] * 28411115059L) + ((int128)tmp_q[6] * 33279420926L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

