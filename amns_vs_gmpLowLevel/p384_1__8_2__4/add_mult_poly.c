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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[7]) << 1);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[7] * pa[6]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[7] * pa[7]) << 1);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3986862024429470917UL) + ((((uint64_t)op[1] * 12545084410724456261UL) + ((uint64_t)op[2] * 15396245711359731548UL) + ((uint64_t)op[3] * 1673052076892358590UL) + ((uint64_t)op[4] * 8781318811651109556UL) + ((uint64_t)op[5] * 17824143636338083919UL) + ((uint64_t)op[6] * 10824372135841103491UL) + ((uint64_t)op[7] * 8218314441052219590UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 8218314441052219590UL) + ((uint64_t)op[1] * 3986862024429470917UL) + ((((uint64_t)op[2] * 12545084410724456261UL) + ((uint64_t)op[3] * 15396245711359731548UL) + ((uint64_t)op[4] * 1673052076892358590UL) + ((uint64_t)op[5] * 8781318811651109556UL) + ((uint64_t)op[6] * 17824143636338083919UL) + ((uint64_t)op[7] * 10824372135841103491UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 10824372135841103491UL) + ((uint64_t)op[1] * 8218314441052219590UL) + ((uint64_t)op[2] * 3986862024429470917UL) + ((((uint64_t)op[3] * 12545084410724456261UL) + ((uint64_t)op[4] * 15396245711359731548UL) + ((uint64_t)op[5] * 1673052076892358590UL) + ((uint64_t)op[6] * 8781318811651109556UL) + ((uint64_t)op[7] * 17824143636338083919UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 17824143636338083919UL) + ((uint64_t)op[1] * 10824372135841103491UL) + ((uint64_t)op[2] * 8218314441052219590UL) + ((uint64_t)op[3] * 3986862024429470917UL) + ((((uint64_t)op[4] * 12545084410724456261UL) + ((uint64_t)op[5] * 15396245711359731548UL) + ((uint64_t)op[6] * 1673052076892358590UL) + ((uint64_t)op[7] * 8781318811651109556UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 8781318811651109556UL) + ((uint64_t)op[1] * 17824143636338083919UL) + ((uint64_t)op[2] * 10824372135841103491UL) + ((uint64_t)op[3] * 8218314441052219590UL) + ((uint64_t)op[4] * 3986862024429470917UL) + ((((uint64_t)op[5] * 12545084410724456261UL) + ((uint64_t)op[6] * 15396245711359731548UL) + ((uint64_t)op[7] * 1673052076892358590UL)) * 2);
	tmp_q[5] = ((uint64_t)op[0] * 1673052076892358590UL) + ((uint64_t)op[1] * 8781318811651109556UL) + ((uint64_t)op[2] * 17824143636338083919UL) + ((uint64_t)op[3] * 10824372135841103491UL) + ((uint64_t)op[4] * 8218314441052219590UL) + ((uint64_t)op[5] * 3986862024429470917UL) + ((((uint64_t)op[6] * 12545084410724456261UL) + ((uint64_t)op[7] * 15396245711359731548UL)) * 2);
	tmp_q[6] = ((uint64_t)op[0] * 15396245711359731548UL) + ((uint64_t)op[1] * 1673052076892358590UL) + ((uint64_t)op[2] * 8781318811651109556UL) + ((uint64_t)op[3] * 17824143636338083919UL) + ((uint64_t)op[4] * 10824372135841103491UL) + ((uint64_t)op[5] * 8218314441052219590UL) + ((uint64_t)op[6] * 3986862024429470917UL) + ((uint64_t)op[7] * 6643424747739360906UL);
	tmp_q[7] = ((uint64_t)op[0] * 12545084410724456261UL) + ((uint64_t)op[1] * 15396245711359731548UL) + ((uint64_t)op[2] * 1673052076892358590UL) + ((uint64_t)op[3] * 8781318811651109556UL) + ((uint64_t)op[4] * 17824143636338083919UL) + ((uint64_t)op[5] * 10824372135841103491UL) + ((uint64_t)op[6] * 8218314441052219590UL) + ((uint64_t)op[7] * 3986862024429470917UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 69116946384611L) + ((((int128)tmp_q[1] * 159037028416942L) - ((int128)tmp_q[2] * 7040316108362L) + ((int128)tmp_q[3] * 16008927558432L) + ((int128)tmp_q[4] * 33312180538131L) + ((int128)tmp_q[5] * 108357106243595L) - ((int128)tmp_q[6] * 13723089772375L) - ((int128)tmp_q[7] * 3831151373936L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 3831151373936L) + ((int128)tmp_q[1] * 69116946384611L) + ((((int128)tmp_q[2] * 159037028416942L) - ((int128)tmp_q[3] * 7040316108362L) + ((int128)tmp_q[4] * 16008927558432L) + ((int128)tmp_q[5] * 33312180538131L) + ((int128)tmp_q[6] * 108357106243595L) - ((int128)tmp_q[7] * 13723089772375L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 13723089772375L) - ((int128)tmp_q[1] * 3831151373936L) + ((int128)tmp_q[2] * 69116946384611L) + ((((int128)tmp_q[3] * 159037028416942L) - ((int128)tmp_q[4] * 7040316108362L) + ((int128)tmp_q[5] * 16008927558432L) + ((int128)tmp_q[6] * 33312180538131L) + ((int128)tmp_q[7] * 108357106243595L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 108357106243595L) - ((int128)tmp_q[1] * 13723089772375L) - ((int128)tmp_q[2] * 3831151373936L) + ((int128)tmp_q[3] * 69116946384611L) + ((((int128)tmp_q[4] * 159037028416942L) - ((int128)tmp_q[5] * 7040316108362L) + ((int128)tmp_q[6] * 16008927558432L) + ((int128)tmp_q[7] * 33312180538131L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 33312180538131L) + ((int128)tmp_q[1] * 108357106243595L) - ((int128)tmp_q[2] * 13723089772375L) - ((int128)tmp_q[3] * 3831151373936L) + ((int128)tmp_q[4] * 69116946384611L) + ((((int128)tmp_q[5] * 159037028416942L) - ((int128)tmp_q[6] * 7040316108362L) + ((int128)tmp_q[7] * 16008927558432L)) * 2);
	tmp_zero[5] = ((int128)tmp_q[0] * 16008927558432L) + ((int128)tmp_q[1] * 33312180538131L) + ((int128)tmp_q[2] * 108357106243595L) - ((int128)tmp_q[3] * 13723089772375L) - ((int128)tmp_q[4] * 3831151373936L) + ((int128)tmp_q[5] * 69116946384611L) + ((((int128)tmp_q[6] * 159037028416942L) - ((int128)tmp_q[7] * 7040316108362L)) * 2);
	tmp_zero[6] = -((int128)tmp_q[0] * 7040316108362L) + ((int128)tmp_q[1] * 16008927558432L) + ((int128)tmp_q[2] * 33312180538131L) + ((int128)tmp_q[3] * 108357106243595L) - ((int128)tmp_q[4] * 13723089772375L) - ((int128)tmp_q[5] * 3831151373936L) + ((int128)tmp_q[6] * 69116946384611L) + ((int128)tmp_q[7] * 318074056833884L);
	tmp_zero[7] = ((int128)tmp_q[0] * 159037028416942L) - ((int128)tmp_q[1] * 7040316108362L) + ((int128)tmp_q[2] * 16008927558432L) + ((int128)tmp_q[3] * 33312180538131L) + ((int128)tmp_q[4] * 108357106243595L) - ((int128)tmp_q[5] * 13723089772375L) - ((int128)tmp_q[6] * 3831151373936L) + ((int128)tmp_q[7] * 69116946384611L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
}

