
def build_structs_data_file(dir_path, word_size, poly_deg, rho_log2, nb_max_add, big_int_name, big_int, small_int, conv_P0, conv_P1):
	with open(dir_path+"/structs_data.h", "w") as f:
		
		f.write("#ifndef STRUCTS_DATA\n")
		f.write("#define STRUCTS_DATA\n\n\n")
		
		f.write("//~ IMPORTANT : We take 'phi = 1 << WORD_SIZE'\n")
		f.write("#define WORD_SIZE " + str(word_size) + "\n")
		f.write("#define POLY_DEG " + str(poly_deg) + "\n")
		f.write("#define NB_COEFF " + str(poly_deg+1) + "\n")
		f.write("#define NB_ADD_MAX " + str(nb_max_add) + "\n\n")
		
		f.write("#define RHO_LOG2 " + str(rho_log2) + "  // rho = 1 << RHO_LOG2.\n\n")
		
		f.write("typedef " + big_int_name + " " + big_int + ";\n\n")
		
		f.write("//~ representations of the polynomials P0 and P1, used for conversion into the AMNS\n")
		f.write(small_int + " poly_P0[NB_COEFF] = {" + str(conv_P0)[1:-1] + "};\n")
		f.write(small_int + " poly_P1[NB_COEFF] = {" + str(conv_P1)[1:-1] + "};\n\n")
		
		f.write("//~ representations of polynomials Pi, for i=2,...,n-1\n")
		f.write(small_int + " polys_P[(NB_COEFF - 2)][NB_COEFF];\n\n")
		
		f.write("mpz_t modul_p;\n")
		f.write("mpz_t gama_pow[POLY_DEG];\n\n")
		
		f.write("#endif\n\n")
	





