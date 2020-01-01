#include "add_mult_poly.h"


void add_poly(int *rop, int *pa, int *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] + pb[j];
}

void sub_poly(int *rop, int *pa, int *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] - pb[j];
}

void neg_poly(int *rop, int *op){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = -op[j];
}

//~ assumes 'scalar' and/or coeffs of 'op' small enough to avoid an overflow.
void scalar_mult_poly(int *rop, int *op, int scalar){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = scalar * op[j];
}

//~ assumes 'scalar' and/or coeffs of 'op' small enough to avoid an overflow.
void double_poly_coeffs(int *rop, int *op){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = op[j] << 1;
}

//~ assumes 'nb_pos' and/or coeffs of 'op' small enough to avoid an overflow.
void lshift_poly_coeffs(int *rop, int *op, int nb_pos){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = op[j] << nb_pos;
}

//~ Computes pa(X)*pb(X) mod(E)
void mult_mod_poly(int *rop, int *pa, int *pb){

	llong tmp_prod_result[NB_COEFF];
	tmp_prod_result[0] = (llong)pa[0] * pb[0] + (((llong)pa[1] * pb[12] + (llong)pa[2] * pb[11] + (llong)pa[3] * pb[10] + (llong)pa[4] * pb[9] + (llong)pa[5] * pb[8] + (llong)pa[6] * pb[7] + (llong)pa[7] * pb[6] + (llong)pa[8] * pb[5] + (llong)pa[9] * pb[4] + (llong)pa[10] * pb[3] + (llong)pa[11] * pb[2] + (llong)pa[12] * pb[1]) << 1);
	tmp_prod_result[1] = (llong)pa[0] * pb[1] + (llong)pa[1] * pb[0] + (((llong)pa[2] * pb[12] + (llong)pa[3] * pb[11] + (llong)pa[4] * pb[10] + (llong)pa[5] * pb[9] + (llong)pa[6] * pb[8] + (llong)pa[7] * pb[7] + (llong)pa[8] * pb[6] + (llong)pa[9] * pb[5] + (llong)pa[10] * pb[4] + (llong)pa[11] * pb[3] + (llong)pa[12] * pb[2]) << 1);
	tmp_prod_result[2] = (llong)pa[0] * pb[2] + (llong)pa[1] * pb[1] + (llong)pa[2] * pb[0] + (((llong)pa[3] * pb[12] + (llong)pa[4] * pb[11] + (llong)pa[5] * pb[10] + (llong)pa[6] * pb[9] + (llong)pa[7] * pb[8] + (llong)pa[8] * pb[7] + (llong)pa[9] * pb[6] + (llong)pa[10] * pb[5] + (llong)pa[11] * pb[4] + (llong)pa[12] * pb[3]) << 1);
	tmp_prod_result[3] = (llong)pa[0] * pb[3] + (llong)pa[1] * pb[2] + (llong)pa[2] * pb[1] + (llong)pa[3] * pb[0] + (((llong)pa[4] * pb[12] + (llong)pa[5] * pb[11] + (llong)pa[6] * pb[10] + (llong)pa[7] * pb[9] + (llong)pa[8] * pb[8] + (llong)pa[9] * pb[7] + (llong)pa[10] * pb[6] + (llong)pa[11] * pb[5] + (llong)pa[12] * pb[4]) << 1);
	tmp_prod_result[4] = (llong)pa[0] * pb[4] + (llong)pa[1] * pb[3] + (llong)pa[2] * pb[2] + (llong)pa[3] * pb[1] + (llong)pa[4] * pb[0] + (((llong)pa[5] * pb[12] + (llong)pa[6] * pb[11] + (llong)pa[7] * pb[10] + (llong)pa[8] * pb[9] + (llong)pa[9] * pb[8] + (llong)pa[10] * pb[7] + (llong)pa[11] * pb[6] + (llong)pa[12] * pb[5]) << 1);
	tmp_prod_result[5] = (llong)pa[0] * pb[5] + (llong)pa[1] * pb[4] + (llong)pa[2] * pb[3] + (llong)pa[3] * pb[2] + (llong)pa[4] * pb[1] + (llong)pa[5] * pb[0] + (((llong)pa[6] * pb[12] + (llong)pa[7] * pb[11] + (llong)pa[8] * pb[10] + (llong)pa[9] * pb[9] + (llong)pa[10] * pb[8] + (llong)pa[11] * pb[7] + (llong)pa[12] * pb[6]) << 1);
	tmp_prod_result[6] = (llong)pa[0] * pb[6] + (llong)pa[1] * pb[5] + (llong)pa[2] * pb[4] + (llong)pa[3] * pb[3] + (llong)pa[4] * pb[2] + (llong)pa[5] * pb[1] + (llong)pa[6] * pb[0] + (((llong)pa[7] * pb[12] + (llong)pa[8] * pb[11] + (llong)pa[9] * pb[10] + (llong)pa[10] * pb[9] + (llong)pa[11] * pb[8] + (llong)pa[12] * pb[7]) << 1);
	tmp_prod_result[7] = (llong)pa[0] * pb[7] + (llong)pa[1] * pb[6] + (llong)pa[2] * pb[5] + (llong)pa[3] * pb[4] + (llong)pa[4] * pb[3] + (llong)pa[5] * pb[2] + (llong)pa[6] * pb[1] + (llong)pa[7] * pb[0] + (((llong)pa[8] * pb[12] + (llong)pa[9] * pb[11] + (llong)pa[10] * pb[10] + (llong)pa[11] * pb[9] + (llong)pa[12] * pb[8]) << 1);
	tmp_prod_result[8] = (llong)pa[0] * pb[8] + (llong)pa[1] * pb[7] + (llong)pa[2] * pb[6] + (llong)pa[3] * pb[5] + (llong)pa[4] * pb[4] + (llong)pa[5] * pb[3] + (llong)pa[6] * pb[2] + (llong)pa[7] * pb[1] + (llong)pa[8] * pb[0] + (((llong)pa[9] * pb[12] + (llong)pa[10] * pb[11] + (llong)pa[11] * pb[10] + (llong)pa[12] * pb[9]) << 1);
	tmp_prod_result[9] = (llong)pa[0] * pb[9] + (llong)pa[1] * pb[8] + (llong)pa[2] * pb[7] + (llong)pa[3] * pb[6] + (llong)pa[4] * pb[5] + (llong)pa[5] * pb[4] + (llong)pa[6] * pb[3] + (llong)pa[7] * pb[2] + (llong)pa[8] * pb[1] + (llong)pa[9] * pb[0] + (((llong)pa[10] * pb[12] + (llong)pa[11] * pb[11] + (llong)pa[12] * pb[10]) << 1);
	tmp_prod_result[10] = (llong)pa[0] * pb[10] + (llong)pa[1] * pb[9] + (llong)pa[2] * pb[8] + (llong)pa[3] * pb[7] + (llong)pa[4] * pb[6] + (llong)pa[5] * pb[5] + (llong)pa[6] * pb[4] + (llong)pa[7] * pb[3] + (llong)pa[8] * pb[2] + (llong)pa[9] * pb[1] + (llong)pa[10] * pb[0] + (((llong)pa[11] * pb[12] + (llong)pa[12] * pb[11]) << 1);
	tmp_prod_result[11] = (llong)pa[0] * pb[11] + (llong)pa[1] * pb[10] + (llong)pa[2] * pb[9] + (llong)pa[3] * pb[8] + (llong)pa[4] * pb[7] + (llong)pa[5] * pb[6] + (llong)pa[6] * pb[5] + (llong)pa[7] * pb[4] + (llong)pa[8] * pb[3] + (llong)pa[9] * pb[2] + (llong)pa[10] * pb[1] + (llong)pa[11] * pb[0] + (((llong)pa[12] * pb[12]) << 1);
	tmp_prod_result[12] = (llong)pa[0] * pb[12] + (llong)pa[1] * pb[11] + (llong)pa[2] * pb[10] + (llong)pa[3] * pb[9] + (llong)pa[4] * pb[8] + (llong)pa[5] * pb[7] + (llong)pa[6] * pb[6] + (llong)pa[7] * pb[5] + (llong)pa[8] * pb[4] + (llong)pa[9] * pb[3] + (llong)pa[10] * pb[2] + (llong)pa[11] * pb[1] + (llong)pa[12] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(E)
void square_mod_poly(int *rop, int *pa){

	llong tmp_prod_result[NB_COEFF];
	tmp_prod_result[0] = (llong)pa[0] * pa[0] + (((llong)pa[7] * pa[6] + (llong)pa[8] * pa[5] + (llong)pa[9] * pa[4] + (llong)pa[10] * pa[3] + (llong)pa[11] * pa[2] + (llong)pa[12] * pa[1]) << 2);
	tmp_prod_result[1] = (((llong)pa[1] * pa[0]) << 1) + (((((llong)pa[8] * pa[6] + (llong)pa[9] * pa[5] + (llong)pa[10] * pa[4] + (llong)pa[11] * pa[3] + (llong)pa[12] * pa[2]) << 1) + (llong)pa[7] * pa[7]) << 1);
	tmp_prod_result[2] = (((llong)pa[2] * pa[0]) << 1) + (llong)pa[1] * pa[1] + (((llong)pa[8] * pa[7] + (llong)pa[9] * pa[6] + (llong)pa[10] * pa[5] + (llong)pa[11] * pa[4] + (llong)pa[12] * pa[3]) << 2);
	tmp_prod_result[3] = (((llong)pa[2] * pa[1] + (llong)pa[3] * pa[0]) << 1) + (((((llong)pa[9] * pa[7] + (llong)pa[10] * pa[6] + (llong)pa[11] * pa[5] + (llong)pa[12] * pa[4]) << 1) + (llong)pa[8] * pa[8]) << 1);
	tmp_prod_result[4] = (((llong)pa[3] * pa[1] + (llong)pa[4] * pa[0]) << 1) + (llong)pa[2] * pa[2] + (((llong)pa[9] * pa[8] + (llong)pa[10] * pa[7] + (llong)pa[11] * pa[6] + (llong)pa[12] * pa[5]) << 2);
	tmp_prod_result[5] = (((llong)pa[3] * pa[2] + (llong)pa[4] * pa[1] + (llong)pa[5] * pa[0]) << 1) + (((((llong)pa[10] * pa[8] + (llong)pa[11] * pa[7] + (llong)pa[12] * pa[6]) << 1) + (llong)pa[9] * pa[9]) << 1);
	tmp_prod_result[6] = (((llong)pa[4] * pa[2] + (llong)pa[5] * pa[1] + (llong)pa[6] * pa[0]) << 1) + (llong)pa[3] * pa[3] + (((llong)pa[10] * pa[9] + (llong)pa[11] * pa[8] + (llong)pa[12] * pa[7]) << 2);
	tmp_prod_result[7] = (((llong)pa[4] * pa[3] + (llong)pa[5] * pa[2] + (llong)pa[6] * pa[1] + (llong)pa[7] * pa[0]) << 1) + (((((llong)pa[11] * pa[9] + (llong)pa[12] * pa[8]) << 1) + (llong)pa[10] * pa[10]) << 1);
	tmp_prod_result[8] = (((llong)pa[5] * pa[3] + (llong)pa[6] * pa[2] + (llong)pa[7] * pa[1] + (llong)pa[8] * pa[0]) << 1) + (llong)pa[4] * pa[4] + (((llong)pa[11] * pa[10] + (llong)pa[12] * pa[9]) << 2);
	tmp_prod_result[9] = (((llong)pa[5] * pa[4] + (llong)pa[6] * pa[3] + (llong)pa[7] * pa[2] + (llong)pa[8] * pa[1] + (llong)pa[9] * pa[0]) << 1) + (((((llong)pa[12] * pa[10]) << 1) + (llong)pa[11] * pa[11]) << 1);
	tmp_prod_result[10] = (((llong)pa[6] * pa[4] + (llong)pa[7] * pa[3] + (llong)pa[8] * pa[2] + (llong)pa[9] * pa[1] + (llong)pa[10] * pa[0]) << 1) + (llong)pa[5] * pa[5] + (((llong)pa[12] * pa[11]) << 2);
	tmp_prod_result[11] = (((llong)pa[6] * pa[5] + (llong)pa[7] * pa[4] + (llong)pa[8] * pa[3] + (llong)pa[9] * pa[2] + (llong)pa[10] * pa[1] + (llong)pa[11] * pa[0]) << 1) + (((llong)pa[12] * pa[12]) << 1);
	tmp_prod_result[12] = (((llong)pa[7] * pa[5] + (llong)pa[8] * pa[4] + (llong)pa[9] * pa[3] + (llong)pa[10] * pa[2] + (llong)pa[11] * pa[1] + (llong)pa[12] * pa[0]) << 1) + (llong)pa[6] * pa[6];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int *rop, llong *op){

	uint tmp_q[NB_COEFF];
	llong tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint)op[0] * 2419461793UL) + ((((uint)op[1] * 1999110111UL) + ((uint)op[2] * 2862384428UL) + ((uint)op[3] * 444564538UL) + ((uint)op[4] * 3211163653UL) + ((uint)op[5] * 1749535181UL) + ((uint)op[6] * 2317451020UL) + ((uint)op[7] * 2810577607UL) + ((uint)op[8] * 2325814501UL) + ((uint)op[9] * 913268822UL) + ((uint)op[10] * 852435122UL) + ((uint)op[11] * 919492381UL) + ((uint)op[12] * 3573709715UL)) * 2);
	tmp_q[1] = ((uint)op[0] * 3573709715UL) + ((uint)op[1] * 2419461793UL) + ((((uint)op[2] * 1999110111UL) + ((uint)op[3] * 2862384428UL) + ((uint)op[4] * 444564538UL) + ((uint)op[5] * 3211163653UL) + ((uint)op[6] * 1749535181UL) + ((uint)op[7] * 2317451020UL) + ((uint)op[8] * 2810577607UL) + ((uint)op[9] * 2325814501UL) + ((uint)op[10] * 913268822UL) + ((uint)op[11] * 852435122UL) + ((uint)op[12] * 919492381UL)) * 2);
	tmp_q[2] = ((uint)op[0] * 919492381UL) + ((uint)op[1] * 3573709715UL) + ((uint)op[2] * 2419461793UL) + ((((uint)op[3] * 1999110111UL) + ((uint)op[4] * 2862384428UL) + ((uint)op[5] * 444564538UL) + ((uint)op[6] * 3211163653UL) + ((uint)op[7] * 1749535181UL) + ((uint)op[8] * 2317451020UL) + ((uint)op[9] * 2810577607UL) + ((uint)op[10] * 2325814501UL) + ((uint)op[11] * 913268822UL) + ((uint)op[12] * 852435122UL)) * 2);
	tmp_q[3] = ((uint)op[0] * 852435122UL) + ((uint)op[1] * 919492381UL) + ((uint)op[2] * 3573709715UL) + ((uint)op[3] * 2419461793UL) + ((((uint)op[4] * 1999110111UL) + ((uint)op[5] * 2862384428UL) + ((uint)op[6] * 444564538UL) + ((uint)op[7] * 3211163653UL) + ((uint)op[8] * 1749535181UL) + ((uint)op[9] * 2317451020UL) + ((uint)op[10] * 2810577607UL) + ((uint)op[11] * 2325814501UL) + ((uint)op[12] * 913268822UL)) * 2);
	tmp_q[4] = ((uint)op[0] * 913268822UL) + ((uint)op[1] * 852435122UL) + ((uint)op[2] * 919492381UL) + ((uint)op[3] * 3573709715UL) + ((uint)op[4] * 2419461793UL) + ((((uint)op[5] * 1999110111UL) + ((uint)op[6] * 2862384428UL) + ((uint)op[7] * 444564538UL) + ((uint)op[8] * 3211163653UL) + ((uint)op[9] * 1749535181UL) + ((uint)op[10] * 2317451020UL) + ((uint)op[11] * 2810577607UL) + ((uint)op[12] * 2325814501UL)) * 2);
	tmp_q[5] = ((uint)op[0] * 2325814501UL) + ((uint)op[1] * 913268822UL) + ((uint)op[2] * 852435122UL) + ((uint)op[3] * 919492381UL) + ((uint)op[4] * 3573709715UL) + ((uint)op[5] * 2419461793UL) + ((((uint)op[6] * 1999110111UL) + ((uint)op[7] * 2862384428UL) + ((uint)op[8] * 444564538UL) + ((uint)op[9] * 3211163653UL) + ((uint)op[10] * 1749535181UL) + ((uint)op[11] * 2317451020UL) + ((uint)op[12] * 2810577607UL)) * 2);
	tmp_q[6] = ((uint)op[0] * 2810577607UL) + ((uint)op[1] * 2325814501UL) + ((uint)op[2] * 913268822UL) + ((uint)op[3] * 852435122UL) + ((uint)op[4] * 919492381UL) + ((uint)op[5] * 3573709715UL) + ((uint)op[6] * 2419461793UL) + ((((uint)op[7] * 1999110111UL) + ((uint)op[8] * 2862384428UL) + ((uint)op[9] * 444564538UL) + ((uint)op[10] * 3211163653UL) + ((uint)op[11] * 1749535181UL) + ((uint)op[12] * 2317451020UL)) * 2);
	tmp_q[7] = ((uint)op[0] * 2317451020UL) + ((uint)op[1] * 2810577607UL) + ((uint)op[2] * 2325814501UL) + ((uint)op[3] * 913268822UL) + ((uint)op[4] * 852435122UL) + ((uint)op[5] * 919492381UL) + ((uint)op[6] * 3573709715UL) + ((uint)op[7] * 2419461793UL) + ((((uint)op[8] * 1999110111UL) + ((uint)op[9] * 2862384428UL) + ((uint)op[10] * 444564538UL) + ((uint)op[11] * 3211163653UL) + ((uint)op[12] * 1749535181UL)) * 2);
	tmp_q[8] = ((uint)op[0] * 1749535181UL) + ((uint)op[1] * 2317451020UL) + ((uint)op[2] * 2810577607UL) + ((uint)op[3] * 2325814501UL) + ((uint)op[4] * 913268822UL) + ((uint)op[5] * 852435122UL) + ((uint)op[6] * 919492381UL) + ((uint)op[7] * 3573709715UL) + ((uint)op[8] * 2419461793UL) + ((((uint)op[9] * 1999110111UL) + ((uint)op[10] * 2862384428UL) + ((uint)op[11] * 444564538UL) + ((uint)op[12] * 3211163653UL)) * 2);
	tmp_q[9] = ((uint)op[0] * 3211163653UL) + ((uint)op[1] * 1749535181UL) + ((uint)op[2] * 2317451020UL) + ((uint)op[3] * 2810577607UL) + ((uint)op[4] * 2325814501UL) + ((uint)op[5] * 913268822UL) + ((uint)op[6] * 852435122UL) + ((uint)op[7] * 919492381UL) + ((uint)op[8] * 3573709715UL) + ((uint)op[9] * 2419461793UL) + ((((uint)op[10] * 1999110111UL) + ((uint)op[11] * 2862384428UL) + ((uint)op[12] * 444564538UL)) * 2);
	tmp_q[10] = ((uint)op[0] * 444564538UL) + ((uint)op[1] * 3211163653UL) + ((uint)op[2] * 1749535181UL) + ((uint)op[3] * 2317451020UL) + ((uint)op[4] * 2810577607UL) + ((uint)op[5] * 2325814501UL) + ((uint)op[6] * 913268822UL) + ((uint)op[7] * 852435122UL) + ((uint)op[8] * 919492381UL) + ((uint)op[9] * 3573709715UL) + ((uint)op[10] * 2419461793UL) + ((((uint)op[11] * 1999110111UL) + ((uint)op[12] * 2862384428UL)) * 2);
	tmp_q[11] = ((uint)op[0] * 2862384428UL) + ((uint)op[1] * 444564538UL) + ((uint)op[2] * 3211163653UL) + ((uint)op[3] * 1749535181UL) + ((uint)op[4] * 2317451020UL) + ((uint)op[5] * 2810577607UL) + ((uint)op[6] * 2325814501UL) + ((uint)op[7] * 913268822UL) + ((uint)op[8] * 852435122UL) + ((uint)op[9] * 919492381UL) + ((uint)op[10] * 3573709715UL) + ((uint)op[11] * 2419461793UL) + ((uint)op[12] * 3998220222UL);
	tmp_q[12] = ((uint)op[0] * 1999110111UL) + ((uint)op[1] * 2862384428UL) + ((uint)op[2] * 444564538UL) + ((uint)op[3] * 3211163653UL) + ((uint)op[4] * 1749535181UL) + ((uint)op[5] * 2317451020UL) + ((uint)op[6] * 2810577607UL) + ((uint)op[7] * 2325814501UL) + ((uint)op[8] * 913268822UL) + ((uint)op[9] * 852435122UL) + ((uint)op[10] * 919492381UL) + ((uint)op[11] * 3573709715UL) + ((uint)op[12] * 2419461793UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((llong)tmp_q[0] * 339445L) + ((((llong)tmp_q[1] * 109713L) + ((llong)tmp_q[2] * 139921L) + ((llong)tmp_q[3] * 60007L) + ((llong)tmp_q[4] * 16418L) + ((llong)tmp_q[5] * 392928L) + ((llong)tmp_q[6] * 372734L) - ((llong)tmp_q[7] * 90544L) + ((llong)tmp_q[8] * 146431L) + ((llong)tmp_q[9] * 64533L) + ((llong)tmp_q[10] * 54047L) - ((llong)tmp_q[11] * 195776L) - ((llong)tmp_q[12] * 89767L)) * 2);
	tmp_zero[1] = -((llong)tmp_q[0] * 89767L) + ((llong)tmp_q[1] * 339445L) + ((((llong)tmp_q[2] * 109713L) + ((llong)tmp_q[3] * 139921L) + ((llong)tmp_q[4] * 60007L) + ((llong)tmp_q[5] * 16418L) + ((llong)tmp_q[6] * 392928L) + ((llong)tmp_q[7] * 372734L) - ((llong)tmp_q[8] * 90544L) + ((llong)tmp_q[9] * 146431L) + ((llong)tmp_q[10] * 64533L) + ((llong)tmp_q[11] * 54047L) - ((llong)tmp_q[12] * 195776L)) * 2);
	tmp_zero[2] = -((llong)tmp_q[0] * 195776L) - ((llong)tmp_q[1] * 89767L) + ((llong)tmp_q[2] * 339445L) + ((((llong)tmp_q[3] * 109713L) + ((llong)tmp_q[4] * 139921L) + ((llong)tmp_q[5] * 60007L) + ((llong)tmp_q[6] * 16418L) + ((llong)tmp_q[7] * 392928L) + ((llong)tmp_q[8] * 372734L) - ((llong)tmp_q[9] * 90544L) + ((llong)tmp_q[10] * 146431L) + ((llong)tmp_q[11] * 64533L) + ((llong)tmp_q[12] * 54047L)) * 2);
	tmp_zero[3] = ((llong)tmp_q[0] * 54047L) - ((llong)tmp_q[1] * 195776L) - ((llong)tmp_q[2] * 89767L) + ((llong)tmp_q[3] * 339445L) + ((((llong)tmp_q[4] * 109713L) + ((llong)tmp_q[5] * 139921L) + ((llong)tmp_q[6] * 60007L) + ((llong)tmp_q[7] * 16418L) + ((llong)tmp_q[8] * 392928L) + ((llong)tmp_q[9] * 372734L) - ((llong)tmp_q[10] * 90544L) + ((llong)tmp_q[11] * 146431L) + ((llong)tmp_q[12] * 64533L)) * 2);
	tmp_zero[4] = ((llong)tmp_q[0] * 64533L) + ((llong)tmp_q[1] * 54047L) - ((llong)tmp_q[2] * 195776L) - ((llong)tmp_q[3] * 89767L) + ((llong)tmp_q[4] * 339445L) + ((((llong)tmp_q[5] * 109713L) + ((llong)tmp_q[6] * 139921L) + ((llong)tmp_q[7] * 60007L) + ((llong)tmp_q[8] * 16418L) + ((llong)tmp_q[9] * 392928L) + ((llong)tmp_q[10] * 372734L) - ((llong)tmp_q[11] * 90544L) + ((llong)tmp_q[12] * 146431L)) * 2);
	tmp_zero[5] = ((llong)tmp_q[0] * 146431L) + ((llong)tmp_q[1] * 64533L) + ((llong)tmp_q[2] * 54047L) - ((llong)tmp_q[3] * 195776L) - ((llong)tmp_q[4] * 89767L) + ((llong)tmp_q[5] * 339445L) + ((((llong)tmp_q[6] * 109713L) + ((llong)tmp_q[7] * 139921L) + ((llong)tmp_q[8] * 60007L) + ((llong)tmp_q[9] * 16418L) + ((llong)tmp_q[10] * 392928L) + ((llong)tmp_q[11] * 372734L) - ((llong)tmp_q[12] * 90544L)) * 2);
	tmp_zero[6] = -((llong)tmp_q[0] * 90544L) + ((llong)tmp_q[1] * 146431L) + ((llong)tmp_q[2] * 64533L) + ((llong)tmp_q[3] * 54047L) - ((llong)tmp_q[4] * 195776L) - ((llong)tmp_q[5] * 89767L) + ((llong)tmp_q[6] * 339445L) + ((((llong)tmp_q[7] * 109713L) + ((llong)tmp_q[8] * 139921L) + ((llong)tmp_q[9] * 60007L) + ((llong)tmp_q[10] * 16418L) + ((llong)tmp_q[11] * 392928L) + ((llong)tmp_q[12] * 372734L)) * 2);
	tmp_zero[7] = ((llong)tmp_q[0] * 372734L) - ((llong)tmp_q[1] * 90544L) + ((llong)tmp_q[2] * 146431L) + ((llong)tmp_q[3] * 64533L) + ((llong)tmp_q[4] * 54047L) - ((llong)tmp_q[5] * 195776L) - ((llong)tmp_q[6] * 89767L) + ((llong)tmp_q[7] * 339445L) + ((((llong)tmp_q[8] * 109713L) + ((llong)tmp_q[9] * 139921L) + ((llong)tmp_q[10] * 60007L) + ((llong)tmp_q[11] * 16418L) + ((llong)tmp_q[12] * 392928L)) * 2);
	tmp_zero[8] = ((llong)tmp_q[0] * 392928L) + ((llong)tmp_q[1] * 372734L) - ((llong)tmp_q[2] * 90544L) + ((llong)tmp_q[3] * 146431L) + ((llong)tmp_q[4] * 64533L) + ((llong)tmp_q[5] * 54047L) - ((llong)tmp_q[6] * 195776L) - ((llong)tmp_q[7] * 89767L) + ((llong)tmp_q[8] * 339445L) + ((((llong)tmp_q[9] * 109713L) + ((llong)tmp_q[10] * 139921L) + ((llong)tmp_q[11] * 60007L) + ((llong)tmp_q[12] * 16418L)) * 2);
	tmp_zero[9] = ((llong)tmp_q[0] * 16418L) + ((llong)tmp_q[1] * 392928L) + ((llong)tmp_q[2] * 372734L) - ((llong)tmp_q[3] * 90544L) + ((llong)tmp_q[4] * 146431L) + ((llong)tmp_q[5] * 64533L) + ((llong)tmp_q[6] * 54047L) - ((llong)tmp_q[7] * 195776L) - ((llong)tmp_q[8] * 89767L) + ((llong)tmp_q[9] * 339445L) + ((((llong)tmp_q[10] * 109713L) + ((llong)tmp_q[11] * 139921L) + ((llong)tmp_q[12] * 60007L)) * 2);
	tmp_zero[10] = ((llong)tmp_q[0] * 60007L) + ((llong)tmp_q[1] * 16418L) + ((llong)tmp_q[2] * 392928L) + ((llong)tmp_q[3] * 372734L) - ((llong)tmp_q[4] * 90544L) + ((llong)tmp_q[5] * 146431L) + ((llong)tmp_q[6] * 64533L) + ((llong)tmp_q[7] * 54047L) - ((llong)tmp_q[8] * 195776L) - ((llong)tmp_q[9] * 89767L) + ((llong)tmp_q[10] * 339445L) + ((((llong)tmp_q[11] * 109713L) + ((llong)tmp_q[12] * 139921L)) * 2);
	tmp_zero[11] = ((llong)tmp_q[0] * 139921L) + ((llong)tmp_q[1] * 60007L) + ((llong)tmp_q[2] * 16418L) + ((llong)tmp_q[3] * 392928L) + ((llong)tmp_q[4] * 372734L) - ((llong)tmp_q[5] * 90544L) + ((llong)tmp_q[6] * 146431L) + ((llong)tmp_q[7] * 64533L) + ((llong)tmp_q[8] * 54047L) - ((llong)tmp_q[9] * 195776L) - ((llong)tmp_q[10] * 89767L) + ((llong)tmp_q[11] * 339445L) + ((llong)tmp_q[12] * 219426L);
	tmp_zero[12] = ((llong)tmp_q[0] * 109713L) + ((llong)tmp_q[1] * 139921L) + ((llong)tmp_q[2] * 60007L) + ((llong)tmp_q[3] * 16418L) + ((llong)tmp_q[4] * 392928L) + ((llong)tmp_q[5] * 372734L) - ((llong)tmp_q[6] * 90544L) + ((llong)tmp_q[7] * 146431L) + ((llong)tmp_q[8] * 64533L) + ((llong)tmp_q[9] * 54047L) - ((llong)tmp_q[10] * 195776L) - ((llong)tmp_q[11] * 89767L) + ((llong)tmp_q[12] * 339445L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
	rop[8] = (op[8] + tmp_zero[8]) >> WORD_SIZE;
	rop[9] = (op[9] + tmp_zero[9]) >> WORD_SIZE;
	rop[10] = (op[10] + tmp_zero[10]) >> WORD_SIZE;
	rop[11] = (op[11] + tmp_zero[11]) >> WORD_SIZE;
	rop[12] = (op[12] + tmp_zero[12]) >> WORD_SIZE;
}

void exact_coeffs_reduction(int *rop, int *op){

	int i;
	llong tmp[NB_COEFF];

	for(i=0; i<NB_COEFF; i++)
		tmp[i] = (llong) op[i];

	internal_reduction(rop, tmp);

	tmp[0] = (llong)rop[0] * poly_P0[0] + (((llong)rop[1] * poly_P0[12] + (llong)rop[2] * poly_P0[11] + (llong)rop[3] * poly_P0[10] + (llong)rop[4] * poly_P0[9] + (llong)rop[5] * poly_P0[8] + (llong)rop[6] * poly_P0[7] + (llong)rop[7] * poly_P0[6] + (llong)rop[8] * poly_P0[5] + (llong)rop[9] * poly_P0[4] + (llong)rop[10] * poly_P0[3] + (llong)rop[11] * poly_P0[2] + (llong)rop[12] * poly_P0[1]) << 1);
	tmp[1] = (llong)rop[0] * poly_P0[1] + (llong)rop[1] * poly_P0[0] + (((llong)rop[2] * poly_P0[12] + (llong)rop[3] * poly_P0[11] + (llong)rop[4] * poly_P0[10] + (llong)rop[5] * poly_P0[9] + (llong)rop[6] * poly_P0[8] + (llong)rop[7] * poly_P0[7] + (llong)rop[8] * poly_P0[6] + (llong)rop[9] * poly_P0[5] + (llong)rop[10] * poly_P0[4] + (llong)rop[11] * poly_P0[3] + (llong)rop[12] * poly_P0[2]) << 1);
	tmp[2] = (llong)rop[0] * poly_P0[2] + (llong)rop[1] * poly_P0[1] + (llong)rop[2] * poly_P0[0] + (((llong)rop[3] * poly_P0[12] + (llong)rop[4] * poly_P0[11] + (llong)rop[5] * poly_P0[10] + (llong)rop[6] * poly_P0[9] + (llong)rop[7] * poly_P0[8] + (llong)rop[8] * poly_P0[7] + (llong)rop[9] * poly_P0[6] + (llong)rop[10] * poly_P0[5] + (llong)rop[11] * poly_P0[4] + (llong)rop[12] * poly_P0[3]) << 1);
	tmp[3] = (llong)rop[0] * poly_P0[3] + (llong)rop[1] * poly_P0[2] + (llong)rop[2] * poly_P0[1] + (llong)rop[3] * poly_P0[0] + (((llong)rop[4] * poly_P0[12] + (llong)rop[5] * poly_P0[11] + (llong)rop[6] * poly_P0[10] + (llong)rop[7] * poly_P0[9] + (llong)rop[8] * poly_P0[8] + (llong)rop[9] * poly_P0[7] + (llong)rop[10] * poly_P0[6] + (llong)rop[11] * poly_P0[5] + (llong)rop[12] * poly_P0[4]) << 1);
	tmp[4] = (llong)rop[0] * poly_P0[4] + (llong)rop[1] * poly_P0[3] + (llong)rop[2] * poly_P0[2] + (llong)rop[3] * poly_P0[1] + (llong)rop[4] * poly_P0[0] + (((llong)rop[5] * poly_P0[12] + (llong)rop[6] * poly_P0[11] + (llong)rop[7] * poly_P0[10] + (llong)rop[8] * poly_P0[9] + (llong)rop[9] * poly_P0[8] + (llong)rop[10] * poly_P0[7] + (llong)rop[11] * poly_P0[6] + (llong)rop[12] * poly_P0[5]) << 1);
	tmp[5] = (llong)rop[0] * poly_P0[5] + (llong)rop[1] * poly_P0[4] + (llong)rop[2] * poly_P0[3] + (llong)rop[3] * poly_P0[2] + (llong)rop[4] * poly_P0[1] + (llong)rop[5] * poly_P0[0] + (((llong)rop[6] * poly_P0[12] + (llong)rop[7] * poly_P0[11] + (llong)rop[8] * poly_P0[10] + (llong)rop[9] * poly_P0[9] + (llong)rop[10] * poly_P0[8] + (llong)rop[11] * poly_P0[7] + (llong)rop[12] * poly_P0[6]) << 1);
	tmp[6] = (llong)rop[0] * poly_P0[6] + (llong)rop[1] * poly_P0[5] + (llong)rop[2] * poly_P0[4] + (llong)rop[3] * poly_P0[3] + (llong)rop[4] * poly_P0[2] + (llong)rop[5] * poly_P0[1] + (llong)rop[6] * poly_P0[0] + (((llong)rop[7] * poly_P0[12] + (llong)rop[8] * poly_P0[11] + (llong)rop[9] * poly_P0[10] + (llong)rop[10] * poly_P0[9] + (llong)rop[11] * poly_P0[8] + (llong)rop[12] * poly_P0[7]) << 1);
	tmp[7] = (llong)rop[0] * poly_P0[7] + (llong)rop[1] * poly_P0[6] + (llong)rop[2] * poly_P0[5] + (llong)rop[3] * poly_P0[4] + (llong)rop[4] * poly_P0[3] + (llong)rop[5] * poly_P0[2] + (llong)rop[6] * poly_P0[1] + (llong)rop[7] * poly_P0[0] + (((llong)rop[8] * poly_P0[12] + (llong)rop[9] * poly_P0[11] + (llong)rop[10] * poly_P0[10] + (llong)rop[11] * poly_P0[9] + (llong)rop[12] * poly_P0[8]) << 1);
	tmp[8] = (llong)rop[0] * poly_P0[8] + (llong)rop[1] * poly_P0[7] + (llong)rop[2] * poly_P0[6] + (llong)rop[3] * poly_P0[5] + (llong)rop[4] * poly_P0[4] + (llong)rop[5] * poly_P0[3] + (llong)rop[6] * poly_P0[2] + (llong)rop[7] * poly_P0[1] + (llong)rop[8] * poly_P0[0] + (((llong)rop[9] * poly_P0[12] + (llong)rop[10] * poly_P0[11] + (llong)rop[11] * poly_P0[10] + (llong)rop[12] * poly_P0[9]) << 1);
	tmp[9] = (llong)rop[0] * poly_P0[9] + (llong)rop[1] * poly_P0[8] + (llong)rop[2] * poly_P0[7] + (llong)rop[3] * poly_P0[6] + (llong)rop[4] * poly_P0[5] + (llong)rop[5] * poly_P0[4] + (llong)rop[6] * poly_P0[3] + (llong)rop[7] * poly_P0[2] + (llong)rop[8] * poly_P0[1] + (llong)rop[9] * poly_P0[0] + (((llong)rop[10] * poly_P0[12] + (llong)rop[11] * poly_P0[11] + (llong)rop[12] * poly_P0[10]) << 1);
	tmp[10] = (llong)rop[0] * poly_P0[10] + (llong)rop[1] * poly_P0[9] + (llong)rop[2] * poly_P0[8] + (llong)rop[3] * poly_P0[7] + (llong)rop[4] * poly_P0[6] + (llong)rop[5] * poly_P0[5] + (llong)rop[6] * poly_P0[4] + (llong)rop[7] * poly_P0[3] + (llong)rop[8] * poly_P0[2] + (llong)rop[9] * poly_P0[1] + (llong)rop[10] * poly_P0[0] + (((llong)rop[11] * poly_P0[12] + (llong)rop[12] * poly_P0[11]) << 1);
	tmp[11] = (llong)rop[0] * poly_P0[11] + (llong)rop[1] * poly_P0[10] + (llong)rop[2] * poly_P0[9] + (llong)rop[3] * poly_P0[8] + (llong)rop[4] * poly_P0[7] + (llong)rop[5] * poly_P0[6] + (llong)rop[6] * poly_P0[5] + (llong)rop[7] * poly_P0[4] + (llong)rop[8] * poly_P0[3] + (llong)rop[9] * poly_P0[2] + (llong)rop[10] * poly_P0[1] + (llong)rop[11] * poly_P0[0] + (((llong)rop[12] * poly_P0[12]) << 1);
	tmp[12] = (llong)rop[0] * poly_P0[12] + (llong)rop[1] * poly_P0[11] + (llong)rop[2] * poly_P0[10] + (llong)rop[3] * poly_P0[9] + (llong)rop[4] * poly_P0[8] + (llong)rop[5] * poly_P0[7] + (llong)rop[6] * poly_P0[6] + (llong)rop[7] * poly_P0[5] + (llong)rop[8] * poly_P0[4] + (llong)rop[9] * poly_P0[3] + (llong)rop[10] * poly_P0[2] + (llong)rop[11] * poly_P0[1] + (llong)rop[12] * poly_P0[0];

	internal_reduction(rop, tmp);
}

