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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 15928871211480478635UL) + ((((uint64_t)op[1] * 2330661720382525752UL) + ((uint64_t)op[2] * 7601590303682541292UL) + ((uint64_t)op[3] * 2345750929951063538UL) + ((uint64_t)op[4] * 5884611759837138603UL) + ((uint64_t)op[5] * 6352150653096248311UL) + ((uint64_t)op[6] * 3612464145882324109UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 3612464145882324109UL) + ((uint64_t)op[1] * 15928871211480478635UL) + ((((uint64_t)op[2] * 2330661720382525752UL) + ((uint64_t)op[3] * 7601590303682541292UL) + ((uint64_t)op[4] * 2345750929951063538UL) + ((uint64_t)op[5] * 5884611759837138603UL) + ((uint64_t)op[6] * 6352150653096248311UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 6352150653096248311UL) + ((uint64_t)op[1] * 3612464145882324109UL) + ((uint64_t)op[2] * 15928871211480478635UL) + ((((uint64_t)op[3] * 2330661720382525752UL) + ((uint64_t)op[4] * 7601590303682541292UL) + ((uint64_t)op[5] * 2345750929951063538UL) + ((uint64_t)op[6] * 5884611759837138603UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 5884611759837138603UL) + ((uint64_t)op[1] * 6352150653096248311UL) + ((uint64_t)op[2] * 3612464145882324109UL) + ((uint64_t)op[3] * 15928871211480478635UL) + ((((uint64_t)op[4] * 2330661720382525752UL) + ((uint64_t)op[5] * 7601590303682541292UL) + ((uint64_t)op[6] * 2345750929951063538UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 2345750929951063538UL) + ((uint64_t)op[1] * 5884611759837138603UL) + ((uint64_t)op[2] * 6352150653096248311UL) + ((uint64_t)op[3] * 3612464145882324109UL) + ((uint64_t)op[4] * 15928871211480478635UL) + ((((uint64_t)op[5] * 2330661720382525752UL) + ((uint64_t)op[6] * 7601590303682541292UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 7601590303682541292UL) + ((uint64_t)op[1] * 2345750929951063538UL) + ((uint64_t)op[2] * 5884611759837138603UL) + ((uint64_t)op[3] * 6352150653096248311UL) + ((uint64_t)op[4] * 3612464145882324109UL) + ((uint64_t)op[5] * 15928871211480478635UL) + ((uint64_t)op[6] * 13785420632944500112UL);
	tmp_q[6] = ((uint64_t)op[0] * 2330661720382525752UL) + ((uint64_t)op[1] * 7601590303682541292UL) + ((uint64_t)op[2] * 2345750929951063538UL) + ((uint64_t)op[3] * 5884611759837138603UL) + ((uint64_t)op[4] * 6352150653096248311UL) + ((uint64_t)op[5] * 3612464145882324109UL) + ((uint64_t)op[6] * 15928871211480478635UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 8927860460542149L) - ((((int128)tmp_q[1] * 11019096048996156L) + ((int128)tmp_q[2] * 13532039508611349L) - ((int128)tmp_q[3] * 4491823099680165L) - ((int128)tmp_q[4] * 10570517433743390L) - ((int128)tmp_q[5] * 5840413509264274L) - ((int128)tmp_q[6] * 5371230250152557L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 5371230250152557L) + ((int128)tmp_q[1] * 8927860460542149L) - ((((int128)tmp_q[2] * 11019096048996156L) + ((int128)tmp_q[3] * 13532039508611349L) - ((int128)tmp_q[4] * 4491823099680165L) - ((int128)tmp_q[5] * 10570517433743390L) - ((int128)tmp_q[6] * 5840413509264274L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 5840413509264274L) - ((int128)tmp_q[1] * 5371230250152557L) + ((int128)tmp_q[2] * 8927860460542149L) - ((((int128)tmp_q[3] * 11019096048996156L) + ((int128)tmp_q[4] * 13532039508611349L) - ((int128)tmp_q[5] * 4491823099680165L) - ((int128)tmp_q[6] * 10570517433743390L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 10570517433743390L) - ((int128)tmp_q[1] * 5840413509264274L) - ((int128)tmp_q[2] * 5371230250152557L) + ((int128)tmp_q[3] * 8927860460542149L) - ((((int128)tmp_q[4] * 11019096048996156L) + ((int128)tmp_q[5] * 13532039508611349L) - ((int128)tmp_q[6] * 4491823099680165L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 4491823099680165L) - ((int128)tmp_q[1] * 10570517433743390L) - ((int128)tmp_q[2] * 5840413509264274L) - ((int128)tmp_q[3] * 5371230250152557L) + ((int128)tmp_q[4] * 8927860460542149L) - ((((int128)tmp_q[5] * 11019096048996156L) + ((int128)tmp_q[6] * 13532039508611349L)) * 2);
	tmp_zero[5] = ((int128)tmp_q[0] * 13532039508611349L) - ((int128)tmp_q[1] * 4491823099680165L) - ((int128)tmp_q[2] * 10570517433743390L) - ((int128)tmp_q[3] * 5840413509264274L) - ((int128)tmp_q[4] * 5371230250152557L) + ((int128)tmp_q[5] * 8927860460542149L) - ((int128)tmp_q[6] * 22038192097992312L);
	tmp_zero[6] = ((int128)tmp_q[0] * 11019096048996156L) + ((int128)tmp_q[1] * 13532039508611349L) - ((int128)tmp_q[2] * 4491823099680165L) - ((int128)tmp_q[3] * 10570517433743390L) - ((int128)tmp_q[4] * 5840413509264274L) - ((int128)tmp_q[5] * 5371230250152557L) + ((int128)tmp_q[6] * 8927860460542149L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

