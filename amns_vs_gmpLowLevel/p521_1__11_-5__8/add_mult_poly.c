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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 5);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 5);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 5);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[10]) * 5);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 10);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 5);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[10] * pa[9]) * 10);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[10] * pa[10]) * 5);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 4748985873728691325UL) + ((((uint64_t)op[1] * 10498079930589236068UL) + ((uint64_t)op[2] * 15970145461866830186UL) + ((uint64_t)op[3] * 10828751350178080950UL) + ((uint64_t)op[4] * 15597870121506163345UL) + ((uint64_t)op[5] * 9028625076367531239UL) + ((uint64_t)op[6] * 2406721272165513848UL) + ((uint64_t)op[7] * 16262024427402024307UL) + ((uint64_t)op[8] * 8623119534880492327UL) + ((uint64_t)op[9] * 11377748919695545531UL) + ((uint64_t)op[10] * 3864368310349030593UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 3864368310349030593UL) + ((uint64_t)op[1] * 4748985873728691325UL) + ((((uint64_t)op[2] * 10498079930589236068UL) + ((uint64_t)op[3] * 15970145461866830186UL) + ((uint64_t)op[4] * 10828751350178080950UL) + ((uint64_t)op[5] * 15597870121506163345UL) + ((uint64_t)op[6] * 9028625076367531239UL) + ((uint64_t)op[7] * 2406721272165513848UL) + ((uint64_t)op[8] * 16262024427402024307UL) + ((uint64_t)op[9] * 8623119534880492327UL) + ((uint64_t)op[10] * 11377748919695545531UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 11377748919695545531UL) + ((uint64_t)op[1] * 3864368310349030593UL) + ((uint64_t)op[2] * 4748985873728691325UL) + ((((uint64_t)op[3] * 10498079930589236068UL) + ((uint64_t)op[4] * 15970145461866830186UL) + ((uint64_t)op[5] * 10828751350178080950UL) + ((uint64_t)op[6] * 15597870121506163345UL) + ((uint64_t)op[7] * 9028625076367531239UL) + ((uint64_t)op[8] * 2406721272165513848UL) + ((uint64_t)op[9] * 16262024427402024307UL) + ((uint64_t)op[10] * 8623119534880492327UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 8623119534880492327UL) + ((uint64_t)op[1] * 11377748919695545531UL) + ((uint64_t)op[2] * 3864368310349030593UL) + ((uint64_t)op[3] * 4748985873728691325UL) + ((((uint64_t)op[4] * 10498079930589236068UL) + ((uint64_t)op[5] * 15970145461866830186UL) + ((uint64_t)op[6] * 10828751350178080950UL) + ((uint64_t)op[7] * 15597870121506163345UL) + ((uint64_t)op[8] * 9028625076367531239UL) + ((uint64_t)op[9] * 2406721272165513848UL) + ((uint64_t)op[10] * 16262024427402024307UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 16262024427402024307UL) + ((uint64_t)op[1] * 8623119534880492327UL) + ((uint64_t)op[2] * 11377748919695545531UL) + ((uint64_t)op[3] * 3864368310349030593UL) + ((uint64_t)op[4] * 4748985873728691325UL) + ((((uint64_t)op[5] * 10498079930589236068UL) + ((uint64_t)op[6] * 15970145461866830186UL) + ((uint64_t)op[7] * 10828751350178080950UL) + ((uint64_t)op[8] * 15597870121506163345UL) + ((uint64_t)op[9] * 9028625076367531239UL) + ((uint64_t)op[10] * 2406721272165513848UL)) * 18446744073709551611);
	tmp_q[5] = ((uint64_t)op[0] * 2406721272165513848UL) + ((uint64_t)op[1] * 16262024427402024307UL) + ((uint64_t)op[2] * 8623119534880492327UL) + ((uint64_t)op[3] * 11377748919695545531UL) + ((uint64_t)op[4] * 3864368310349030593UL) + ((uint64_t)op[5] * 4748985873728691325UL) + ((((uint64_t)op[6] * 10498079930589236068UL) + ((uint64_t)op[7] * 15970145461866830186UL) + ((uint64_t)op[8] * 10828751350178080950UL) + ((uint64_t)op[9] * 15597870121506163345UL) + ((uint64_t)op[10] * 9028625076367531239UL)) * 18446744073709551611);
	tmp_q[6] = ((uint64_t)op[0] * 9028625076367531239UL) + ((uint64_t)op[1] * 2406721272165513848UL) + ((uint64_t)op[2] * 16262024427402024307UL) + ((uint64_t)op[3] * 8623119534880492327UL) + ((uint64_t)op[4] * 11377748919695545531UL) + ((uint64_t)op[5] * 3864368310349030593UL) + ((uint64_t)op[6] * 4748985873728691325UL) + ((((uint64_t)op[7] * 10498079930589236068UL) + ((uint64_t)op[8] * 15970145461866830186UL) + ((uint64_t)op[9] * 10828751350178080950UL) + ((uint64_t)op[10] * 15597870121506163345UL)) * 18446744073709551611);
	tmp_q[7] = ((uint64_t)op[0] * 15597870121506163345UL) + ((uint64_t)op[1] * 9028625076367531239UL) + ((uint64_t)op[2] * 2406721272165513848UL) + ((uint64_t)op[3] * 16262024427402024307UL) + ((uint64_t)op[4] * 8623119534880492327UL) + ((uint64_t)op[5] * 11377748919695545531UL) + ((uint64_t)op[6] * 3864368310349030593UL) + ((uint64_t)op[7] * 4748985873728691325UL) + ((((uint64_t)op[8] * 10498079930589236068UL) + ((uint64_t)op[9] * 15970145461866830186UL) + ((uint64_t)op[10] * 10828751350178080950UL)) * 18446744073709551611);
	tmp_q[8] = ((uint64_t)op[0] * 10828751350178080950UL) + ((uint64_t)op[1] * 15597870121506163345UL) + ((uint64_t)op[2] * 9028625076367531239UL) + ((uint64_t)op[3] * 2406721272165513848UL) + ((uint64_t)op[4] * 16262024427402024307UL) + ((uint64_t)op[5] * 8623119534880492327UL) + ((uint64_t)op[6] * 11377748919695545531UL) + ((uint64_t)op[7] * 3864368310349030593UL) + ((uint64_t)op[8] * 4748985873728691325UL) + ((((uint64_t)op[9] * 10498079930589236068UL) + ((uint64_t)op[10] * 15970145461866830186UL)) * 18446744073709551611);
	tmp_q[9] = ((uint64_t)op[0] * 15970145461866830186UL) + ((uint64_t)op[1] * 10828751350178080950UL) + ((uint64_t)op[2] * 15597870121506163345UL) + ((uint64_t)op[3] * 9028625076367531239UL) + ((uint64_t)op[4] * 2406721272165513848UL) + ((uint64_t)op[5] * 16262024427402024307UL) + ((uint64_t)op[6] * 8623119534880492327UL) + ((uint64_t)op[7] * 11377748919695545531UL) + ((uint64_t)op[8] * 3864368310349030593UL) + ((uint64_t)op[9] * 4748985873728691325UL) + ((uint64_t)op[10] * 2849832568182474508UL);
	tmp_q[10] = ((uint64_t)op[0] * 10498079930589236068UL) + ((uint64_t)op[1] * 15970145461866830186UL) + ((uint64_t)op[2] * 10828751350178080950UL) + ((uint64_t)op[3] * 15597870121506163345UL) + ((uint64_t)op[4] * 9028625076367531239UL) + ((uint64_t)op[5] * 2406721272165513848UL) + ((uint64_t)op[6] * 16262024427402024307UL) + ((uint64_t)op[7] * 8623119534880492327UL) + ((uint64_t)op[8] * 11377748919695545531UL) + ((uint64_t)op[9] * 3864368310349030593UL) + ((uint64_t)op[10] * 4748985873728691325UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 5769364062286L) - ((-((int128)tmp_q[1] * 39885879277319L) - ((int128)tmp_q[2] * 30479918027186L) + ((int128)tmp_q[3] * 50019844059288L) - ((int128)tmp_q[4] * 5837994208687L) + ((int128)tmp_q[5] * 66994899244386L) - ((int128)tmp_q[6] * 14571984830463L) + ((int128)tmp_q[7] * 95095792590062L) + ((int128)tmp_q[8] * 108829942252234L) - ((int128)tmp_q[9] * 109926844467260L) - ((int128)tmp_q[10] * 64875712781902L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 64875712781902L) - ((int128)tmp_q[1] * 5769364062286L) - ((-((int128)tmp_q[2] * 39885879277319L) - ((int128)tmp_q[3] * 30479918027186L) + ((int128)tmp_q[4] * 50019844059288L) - ((int128)tmp_q[5] * 5837994208687L) + ((int128)tmp_q[6] * 66994899244386L) - ((int128)tmp_q[7] * 14571984830463L) + ((int128)tmp_q[8] * 95095792590062L) + ((int128)tmp_q[9] * 108829942252234L) - ((int128)tmp_q[10] * 109926844467260L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 109926844467260L) - ((int128)tmp_q[1] * 64875712781902L) - ((int128)tmp_q[2] * 5769364062286L) - ((-((int128)tmp_q[3] * 39885879277319L) - ((int128)tmp_q[4] * 30479918027186L) + ((int128)tmp_q[5] * 50019844059288L) - ((int128)tmp_q[6] * 5837994208687L) + ((int128)tmp_q[7] * 66994899244386L) - ((int128)tmp_q[8] * 14571984830463L) + ((int128)tmp_q[9] * 95095792590062L) + ((int128)tmp_q[10] * 108829942252234L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 108829942252234L) - ((int128)tmp_q[1] * 109926844467260L) - ((int128)tmp_q[2] * 64875712781902L) - ((int128)tmp_q[3] * 5769364062286L) - ((-((int128)tmp_q[4] * 39885879277319L) - ((int128)tmp_q[5] * 30479918027186L) + ((int128)tmp_q[6] * 50019844059288L) - ((int128)tmp_q[7] * 5837994208687L) + ((int128)tmp_q[8] * 66994899244386L) - ((int128)tmp_q[9] * 14571984830463L) + ((int128)tmp_q[10] * 95095792590062L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 95095792590062L) + ((int128)tmp_q[1] * 108829942252234L) - ((int128)tmp_q[2] * 109926844467260L) - ((int128)tmp_q[3] * 64875712781902L) - ((int128)tmp_q[4] * 5769364062286L) - ((-((int128)tmp_q[5] * 39885879277319L) - ((int128)tmp_q[6] * 30479918027186L) + ((int128)tmp_q[7] * 50019844059288L) - ((int128)tmp_q[8] * 5837994208687L) + ((int128)tmp_q[9] * 66994899244386L) - ((int128)tmp_q[10] * 14571984830463L)) * 5);
	tmp_zero[5] = -((int128)tmp_q[0] * 14571984830463L) + ((int128)tmp_q[1] * 95095792590062L) + ((int128)tmp_q[2] * 108829942252234L) - ((int128)tmp_q[3] * 109926844467260L) - ((int128)tmp_q[4] * 64875712781902L) - ((int128)tmp_q[5] * 5769364062286L) - ((-((int128)tmp_q[6] * 39885879277319L) - ((int128)tmp_q[7] * 30479918027186L) + ((int128)tmp_q[8] * 50019844059288L) - ((int128)tmp_q[9] * 5837994208687L) + ((int128)tmp_q[10] * 66994899244386L)) * 5);
	tmp_zero[6] = ((int128)tmp_q[0] * 66994899244386L) - ((int128)tmp_q[1] * 14571984830463L) + ((int128)tmp_q[2] * 95095792590062L) + ((int128)tmp_q[3] * 108829942252234L) - ((int128)tmp_q[4] * 109926844467260L) - ((int128)tmp_q[5] * 64875712781902L) - ((int128)tmp_q[6] * 5769364062286L) - ((-((int128)tmp_q[7] * 39885879277319L) - ((int128)tmp_q[8] * 30479918027186L) + ((int128)tmp_q[9] * 50019844059288L) - ((int128)tmp_q[10] * 5837994208687L)) * 5);
	tmp_zero[7] = -((int128)tmp_q[0] * 5837994208687L) + ((int128)tmp_q[1] * 66994899244386L) - ((int128)tmp_q[2] * 14571984830463L) + ((int128)tmp_q[3] * 95095792590062L) + ((int128)tmp_q[4] * 108829942252234L) - ((int128)tmp_q[5] * 109926844467260L) - ((int128)tmp_q[6] * 64875712781902L) - ((int128)tmp_q[7] * 5769364062286L) - ((-((int128)tmp_q[8] * 39885879277319L) - ((int128)tmp_q[9] * 30479918027186L) + ((int128)tmp_q[10] * 50019844059288L)) * 5);
	tmp_zero[8] = ((int128)tmp_q[0] * 50019844059288L) - ((int128)tmp_q[1] * 5837994208687L) + ((int128)tmp_q[2] * 66994899244386L) - ((int128)tmp_q[3] * 14571984830463L) + ((int128)tmp_q[4] * 95095792590062L) + ((int128)tmp_q[5] * 108829942252234L) - ((int128)tmp_q[6] * 109926844467260L) - ((int128)tmp_q[7] * 64875712781902L) - ((int128)tmp_q[8] * 5769364062286L) - ((-((int128)tmp_q[9] * 39885879277319L) - ((int128)tmp_q[10] * 30479918027186L)) * 5);
	tmp_zero[9] = -((int128)tmp_q[0] * 30479918027186L) + ((int128)tmp_q[1] * 50019844059288L) - ((int128)tmp_q[2] * 5837994208687L) + ((int128)tmp_q[3] * 66994899244386L) - ((int128)tmp_q[4] * 14571984830463L) + ((int128)tmp_q[5] * 95095792590062L) + ((int128)tmp_q[6] * 108829942252234L) - ((int128)tmp_q[7] * 109926844467260L) - ((int128)tmp_q[8] * 64875712781902L) - ((int128)tmp_q[9] * 5769364062286L) + ((int128)tmp_q[10] * 199429396386595L);
	tmp_zero[10] = -((int128)tmp_q[0] * 39885879277319L) - ((int128)tmp_q[1] * 30479918027186L) + ((int128)tmp_q[2] * 50019844059288L) - ((int128)tmp_q[3] * 5837994208687L) + ((int128)tmp_q[4] * 66994899244386L) - ((int128)tmp_q[5] * 14571984830463L) + ((int128)tmp_q[6] * 95095792590062L) + ((int128)tmp_q[7] * 108829942252234L) - ((int128)tmp_q[8] * 109926844467260L) - ((int128)tmp_q[9] * 64875712781902L) - ((int128)tmp_q[10] * 5769364062286L);

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
}

