load("amns_generator.sage")


#~ NOTE: Here, AMNS generation is done according to the implementation strategy mentioned at the end of section 5.1 of the article.

#~ ---------------------------------------------------------------------------------------

word_size = 64

p_size = 256

n_min = (p_size//word_size) + 1 
n_max = n_min + 2

abs_lamb_max = 1 << 3

nth_root_max_duration_checks = 60  # seconds

find_all_nthroot = True


#~ random prime generation
p = random_prime(2**p_size, lbound=2**(p_size-1))
print (p)
print (p.is_prime())
print (p.nbits())


#~ amns generation
flt_amns = generate_amns_candidates_with_n_min_max(word_size, p, n_min, n_max, abs_lamb_max, nth_root_max_duration_checks, find_all_nthroot)



#~ Data structure for each AMNS generated: [delta, n, E, rho_log2, gamma, M, M', conv_P0, conv_P1]


nb_amns = len(flt_amns)

print nb_amns

if nb_amns != 0:
	print
	print flt_amns[0]
	print



































