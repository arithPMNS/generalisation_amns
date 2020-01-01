
def build_add_mult_poly_h_file(dir_path, small_int, big_int):
	with open(dir_path+"/add_mult_poly.h", "w") as f:
		
		f.write("#ifndef POLY_MULT_ADD\n")
		f.write("#define POLY_MULT_ADD\n\n\n")
		
		f.write("void sub_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb);\n")
		f.write("void add_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb);\n")
		f.write("void neg_poly(" + small_int + " *rop, " + small_int + " *op);\n")
		f.write("void scalar_mult_poly(" + small_int + " *rop, " + small_int + " *op, " + small_int + " scalar);\n")
		f.write("void double_poly_coeffs(" + small_int + " *rop, " + small_int + " *op);\n")
		f.write("void lshift_poly_coeffs(" + small_int + " *rop, " + small_int + " *op, int nb_pos);\n\n")
		f.write("void mult_mod_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb);\n\n")
		f.write("void square_mod_poly(" + small_int + " *rop, " + small_int + " *pa);\n\n")	
		f.write("void internal_reduction(" + small_int + " *rop, " + big_int + " *op);\n\n")
		f.write("void exact_coeffs_reduction(" + small_int + " *rop, " + small_int + " *op);\n\n")
		
		f.write("#endif\n\n")


def build_add_mult_poly_c_file(dir_path, n, mont_phi, lambd, small_int, unsigned_small_int, big_int, red_int_coeff, neg_inv_ri_rep_coeff):
	with open(dir_path+"/add_mult_poly.c", "w") as f:
		
		f.write("#include \"add_mult_poly.h\"\n\n\n")
		
		f.write("void add_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb){\n")
		f.write("	int j;\n")
		f.write("	for (j=0; j<NB_COEFF; j++)\n")
		f.write("		rop[j] = pa[j] + pb[j];\n")
		f.write("}\n\n")

		f.write("void sub_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb){\n")
		f.write("	int j;\n")
		f.write("	for (j=0; j<NB_COEFF; j++)\n")
		f.write("		rop[j] = pa[j] - pb[j];\n")
		f.write("}\n\n")

		f.write("void neg_poly(" + small_int + " *rop, " + small_int + " *op){\n")
		f.write("	int j;\n")
		f.write("	for (j=0; j<NB_COEFF; j++)\n")
		f.write("		rop[j] = -op[j];\n")
		f.write("}\n\n")

		f.write("//~ assumes 'scalar' and/or coeffs of 'op' small enough to avoid an overflow.\n")
		f.write("void scalar_mult_poly(" + small_int + " *rop, " + small_int + " *op, " + small_int + " scalar){\n")
		f.write("	int j;\n")
		f.write("	for (j=0; j<NB_COEFF; j++)\n")
		f.write("		rop[j] = scalar * op[j];\n")
		f.write("}\n\n")
		
		f.write("//~ assumes 'scalar' and/or coeffs of 'op' small enough to avoid an overflow.\n")
		f.write("void double_poly_coeffs(" + small_int + " *rop, " + small_int + " *op){\n")
		f.write("	int j;\n")
		f.write("	for (j=0; j<NB_COEFF; j++)\n")
		f.write("		rop[j] = op[j] << 1;\n")
		f.write("}\n\n")
		
		f.write("//~ assumes 'nb_pos' and/or coeffs of 'op' small enough to avoid an overflow.\n")
		f.write("void lshift_poly_coeffs(" + small_int + " *rop, " + small_int + " *op, int nb_pos){\n")
		f.write("	int j;\n")
		f.write("	for (j=0; j<NB_COEFF; j++)\n")
		f.write("		rop[j] = op[j] << nb_pos;\n")
		f.write("}\n\n")
		
		f.write("//~ Computes pa(X)*pb(X) mod(E)\n")
		f.write("void mult_mod_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb){\n\n")
		f.write("	" + big_int + " tmp_prod_result[NB_COEFF];\n")
		f.write(build_prod_code(n, lambd, big_int))
		f.write("\n	internal_reduction(rop, tmp_prod_result);\n")
		f.write("}\n\n")
		
		f.write("//~ Computes pa(X)^2 mod(E)\n")
		f.write("void square_mod_poly(" + small_int + " *rop, " + small_int + " *pa){\n\n")
		f.write("	" + big_int + " tmp_prod_result[NB_COEFF];\n")
		f.write(build_square_code(n, lambd, big_int))
		f.write("\n	internal_reduction(rop, tmp_prod_result);\n")
		f.write("}\n\n")
		
		f.write("//~ performs the internal reduction on 'op' and puts the result in 'rop'\n")
		f.write("//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.\n")
		f.write("void internal_reduction(" + small_int + " *rop, " + big_int + " *op){\n\n")
		f.write("	" + unsigned_small_int + " tmp_q[NB_COEFF];\n")
		f.write("	" + big_int + " tmp_zero[NB_COEFF];\n\n")
		f.write(build_red_int_code(unsigned_small_int, big_int, n, mont_phi, lambd, (lambd%mont_phi), red_int_coeff, neg_inv_ri_rep_coeff))
		f.write("}\n\n")
		
		f.write("void exact_coeffs_reduction(" + small_int + " *rop, " + small_int + " *op){\n\n")
		f.write("	int i;\n")
		f.write("	" + big_int + " tmp[NB_COEFF];\n\n")
		f.write("	for(i=0; i<NB_COEFF; i++)\n")
		f.write("		tmp[i] = (" + big_int + ") op[i];\n")
		f.write("\n")
		f.write("	internal_reduction(rop, tmp);\n")
		f.write(build_exactRedCoeff_interProd_code(n, lambd, big_int))
		f.write("\n")
		f.write("	internal_reduction(rop, tmp);\n")
		f.write("}\n\n")







