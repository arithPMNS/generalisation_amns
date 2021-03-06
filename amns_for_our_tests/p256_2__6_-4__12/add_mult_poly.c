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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8424404232240434215UL) + ((((uint64_t)op[1] * 5367489611988368752UL) + ((uint64_t)op[2] * 1525196241188178378UL) + ((uint64_t)op[3] * 4996612649194714638UL) + ((uint64_t)op[4] * 6890741552638037896UL) + ((uint64_t)op[5] * 15045104381345770818UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 15045104381345770818UL) + ((uint64_t)op[1] * 8424404232240434215UL) + ((((uint64_t)op[2] * 5367489611988368752UL) + ((uint64_t)op[3] * 1525196241188178378UL) + ((uint64_t)op[4] * 4996612649194714638UL) + ((uint64_t)op[5] * 6890741552638037896UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 6890741552638037896UL) + ((uint64_t)op[1] * 15045104381345770818UL) + ((uint64_t)op[2] * 8424404232240434215UL) + ((((uint64_t)op[3] * 5367489611988368752UL) + ((uint64_t)op[4] * 1525196241188178378UL) + ((uint64_t)op[5] * 4996612649194714638UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 4996612649194714638UL) + ((uint64_t)op[1] * 6890741552638037896UL) + ((uint64_t)op[2] * 15045104381345770818UL) + ((uint64_t)op[3] * 8424404232240434215UL) + ((((uint64_t)op[4] * 5367489611988368752UL) + ((uint64_t)op[5] * 1525196241188178378UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 1525196241188178378UL) + ((uint64_t)op[1] * 4996612649194714638UL) + ((uint64_t)op[2] * 6890741552638037896UL) + ((uint64_t)op[3] * 15045104381345770818UL) + ((uint64_t)op[4] * 8424404232240434215UL) + ((uint64_t)op[5] * 15423529699465628224UL);
	tmp_q[5] = ((uint64_t)op[0] * 5367489611988368752UL) + ((uint64_t)op[1] * 1525196241188178378UL) + ((uint64_t)op[2] * 4996612649194714638UL) + ((uint64_t)op[3] * 6890741552638037896UL) + ((uint64_t)op[4] * 15045104381345770818UL) + ((uint64_t)op[5] * 8424404232240434215UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 3103029602567L) - ((-((int128)tmp_q[1] * 1490264390688L) - ((int128)tmp_q[2] * 663537588430L) - ((int128)tmp_q[3] * 3110623879210L) - ((int128)tmp_q[4] * 2347545691364L) + ((int128)tmp_q[5] * 2669310841634L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 2669310841634L) - ((int128)tmp_q[1] * 3103029602567L) - ((-((int128)tmp_q[2] * 1490264390688L) - ((int128)tmp_q[3] * 663537588430L) - ((int128)tmp_q[4] * 3110623879210L) - ((int128)tmp_q[5] * 2347545691364L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 2347545691364L) + ((int128)tmp_q[1] * 2669310841634L) - ((int128)tmp_q[2] * 3103029602567L) - ((-((int128)tmp_q[3] * 1490264390688L) - ((int128)tmp_q[4] * 663537588430L) - ((int128)tmp_q[5] * 3110623879210L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 3110623879210L) - ((int128)tmp_q[1] * 2347545691364L) + ((int128)tmp_q[2] * 2669310841634L) - ((int128)tmp_q[3] * 3103029602567L) - ((-((int128)tmp_q[4] * 1490264390688L) - ((int128)tmp_q[5] * 663537588430L)) * 4);
	tmp_zero[4] = -((int128)tmp_q[0] * 663537588430L) - ((int128)tmp_q[1] * 3110623879210L) - ((int128)tmp_q[2] * 2347545691364L) + ((int128)tmp_q[3] * 2669310841634L) - ((int128)tmp_q[4] * 3103029602567L) + ((int128)tmp_q[5] * 5961057562752L);
	tmp_zero[5] = -((int128)tmp_q[0] * 1490264390688L) - ((int128)tmp_q[1] * 663537588430L) - ((int128)tmp_q[2] * 3110623879210L) - ((int128)tmp_q[3] * 2347545691364L) + ((int128)tmp_q[4] * 2669310841634L) - ((int128)tmp_q[5] * 3103029602567L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

