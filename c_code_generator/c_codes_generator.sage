from os import makedirs as create_dir
from shutil import rmtree as delete_code
load("c_codes_gen_utils.sage")
load("gen_file_main.sage")
load("gen_file_init_amns.sage")
load("gen_file_structs_data.sage")
load("gen_file_add_mult_poly.sage")
load("gen_file_useful_functs.sage")



#~ assumes that polynomial representation form is P(X) = a_0 + ... + a_n.X^n = [a_0, ..., a_n].
#~ assumes 'amns_data' is a list containning (in this order): [delta, n, E, rho_log2, gamma, M, M', conv_P0, conv_P1]
#~ assumes 'target_archi_info' is a list containning (in this order): [word_size, small_int_name, unsigned_small_int_name, big_int_name, big_int_new_name]
#~ WARNING : 'amns_data' is supposed generated for 'target_archi_info'
def build_amns_c_codes(p, amns_data, target_archi_info, p_num=0, amns_num=0):
	
	amnsData = list(amns_data)
	
	redExtPol_coeffs = amnsData[2]
	
	amnsData[2] = -redExtPol_coeffs[0]   # lambda
	
	gen_c_codes(p, amnsData, target_archi_info, p_num, amns_num)
	
	return;




def gen_c_codes(p, amns_data, target_archi_info, p_num, amns_num):
	
	word_size = target_archi_info[0]
	small_int = target_archi_info[1]
	unsigned_small_int = target_archi_info[2]
	big_int_name = target_archi_info[3]
	big_int = target_archi_info[4]
	
	nb_max_add = amns_data[0]
	n = amns_data[1]
	lambd = amns_data[2]
	rho_log2 = amns_data[3]
	gmm = amns_data[4]
	red_int_pol = amns_data[5]
	neg_inv_ri = amns_data[6]
	conv_P0 = amns_data[7]
	conv_P1 = amns_data[8]
	
	p = Integer(p)
	
	dir_name = 'p' + str(p.nbits()) + '_'  + str(p_num) + '__' + str(n) + '_' + str(lambd) + '__' + str(amns_num)
	dir_path = "c_codes/" + dir_name
	try:
		create_dir(dir_path)
	except OSError:  # if this directory already exist
		delete_code(dir_path)
		create_dir(dir_path)
	
	
	mont_phi = 1 << word_size
	
	
	build_main_file(dir_path, small_int)
	
	
	build_structs_data_file(dir_path, word_size, (n-1), rho_log2, nb_max_add, big_int_name, big_int, small_int, conv_P0, conv_P1)

	build_amns_init_h_file(dir_path)
	build_amns_init_c_file(dir_path, p, gmm, small_int)

	build_add_mult_poly_h_file(dir_path, small_int, big_int)
	build_add_mult_poly_c_file(dir_path, n, mont_phi, lambd, small_int, unsigned_small_int, big_int, red_int_pol, neg_inv_ri)

	build_useful_functs_h_file(dir_path, small_int, unsigned_small_int, big_int)
	build_useful_functs_c_file(dir_path, small_int, unsigned_small_int, big_int)
	




















