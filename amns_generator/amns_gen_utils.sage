#~ return the infinite norm of 'vect'
def get_inf_norm(vect):
	tmp = [abs(a) for a in vect]
	return max(tmp)

#~ return an element of 'm' with the smallest infinite norm
def get_smallest_vect(m): 
	sml_vect = m[0]
	sv_norm = get_inf_norm(m[0])
	for v in m[1:] :
		tmp = get_inf_norm(v)
		if sv_norm > tmp :
			sv_norm = tmp
			sml_vect = v
	return list(sml_vect)

#~ -------------- Case where : E(X) = X^n - lambda, with abs(lambda) even ------------------------------

#~ note : each line corresponds to a poly of the form : a0 + a1.X + ... + a_{n-1}.X^{n-1}
def build_lattice_base_even_lambda(p, n, gmm):
	b = []
	l = [p] + [0]*(n-1)
	b.append(l)
	m = identity_matrix(n)
	for i in range(1,n):
		l = [(-gmm.powermod(i, p))%p] + list(m[i][1:])
		b.append(l)
	bb = matrix(b)
	return bb.LLL(delta=1, algorithm='NTL:LLL')

#~ an element 'pol' is valid if 'pol[0] % 2 == 1' 
def get_valid_elements_list_even_lambda(m):
	res = []
	for l in m:
		if l[0]%2 == 1 : # it must be odd for an odd resultant
			res.append(list(l))
	return res

#~ return a small valid vect
def find_valid_vect_even_lambda(p, n, gmm):
	
	m = build_lattice_base_even_lambda(p, n, gmm)
	
	vld_elmts = get_valid_elements_list_even_lambda(m)
	
	if vld_elmts == [] :
		print("WARNING : even lambda, no valid poly found; this shouldn't happen !!!")
	
	return get_smallest_vect(vld_elmts)

#~ -------------- Case where : E(X) = X^n - lambda, with abs(lambda) odd ------------------------------

#~ note : each line corresponds to a poly of the form : a0 + a1.X + ... + a_{n-1}.X^{n-1}
#~ note : 'p' is supposed odd
def build_lattice_base_odd_lambda(p, n, gmm):
	b = []
	l = [p] + [0]*(n-1)
	b.append(l)
	m = identity_matrix(n)
	for i in range(1,n):
		t = (-gmm.powermod(i, p))%p
		if t%2 == 1 :
			t += p
		l = [t] + list(m[i][1:])
		b.append(l)
	bb = matrix(b)
	return bb.LLL(delta=1, algorithm='NTL:LLL')

# assumes : k lower than 1<<n
def build_weight_vect(k, n):
	vt = list(bin(k)[2:])
	vt = [int(x) for x in vt]
	vt = [0]*(n-len(vt)) + vt
	return vt

#~ here, 'vect' is good if ... (see the article) 
def check_vect_odd_lambda(vect, e2, R):
	nb_odds = 0
	vect2 = []
	for x in vect :
		x2 = x % 2
		vect2.append(x2)
		if x2 == 1:
			nb_odds += 1
	if (nb_odds % 2) == 0 :
		#~ no need to continue if there is an even number of odd elements
		return False
	m = R(vect2)
	return (e2.gcd(m) == 1)

#~ return a small valid vect
def find_valid_vect_odd_lambda(p, n, gmm):
	
	R.<x> = GF(2)[]
	
	e2 = R(x**n - 1)
	
	mm = build_lattice_base_odd_lambda(p, n, gmm)
	
	vmax = 1 << n
	vld_vects = []
	for k in range(1, vmax): 
		k_vt = build_weight_vect(k, n)
		vsum = 0*mm[0] # to obtain the good type of vector
		for i in range(n):
			vsum += k_vt[i]*mm[i]
		if check_vect_odd_lambda(vsum, e2, R) :
			vld_vects.append(list(vsum))
	
	if vld_vects != [] :
		return get_smallest_vect(vld_vects)
	
	print("WARNING : odd lambda, no valid poly found; this shouldn't happen !!!")
	
	return []


#~ ---------------------------------------------------------------------------------

#~ note: we assume that elements coeffs can be negatives in the AMNS, so one bit will be allocated for that.
#~ NOTE : 'rho' is taken as a power of two
def compute_rhoUp_log2_and_nb_add_max(word_size, n, lambd, mont_phi, redIntPol_coeffs):
	
	nm_inf = get_inf_norm(redIntPol_coeffs)
	
	if nm_inf == 0: # should not happen, but ...
		return (0,-1)
	
	w = 1 + (n-1)*abs(lambd)
	
	rho_init = 2 * w * nm_inf
	rhoUp_log2 = ceil(log(rho_init, 2))
	rho = 1 << rhoUp_log2
	tmp_phi = 2 * w * rho
	
	coeff_ubound = 1 << (word_size - 1) # one bit is allocated for sign
	prod_acc_ubound = 1 << (2*word_size - 1) # one bit is also allocated for sign
	
	if (rho > coeff_ubound) or  (tmp_phi > mont_phi):
		return (0,-1) # no need to continue
	
	nb_max_add = -1
	while True :  
		
		# note: this loop is ended when 'nb_max_add' becomes too big
		# note: returns always 'nb_max_add-1' because the loop ended when 'nb_max_add' is too big.
		
		nb_max_add += 1
		
		c_maxVal = (nb_max_add + 1)*(rho - 1)
		if coeff_ubound <= c_maxVal :
			return (rhoUp_log2, (nb_max_add-1)) # here, 'rho' (and/or 'nb_max_add') is too big for 'word_size', (one bit is used for sign)
		
		add_effect_on_phi = (nb_max_add + 1)**2
		if mont_phi < (add_effect_on_phi * tmp_phi) :
			return (rhoUp_log2, (nb_max_add-1)) # here, 'rho' (and/or 'nb_max_add') is too big for the montgomery-like reduction to be correct
		
		prod_infNorm_max = w * (nm_inf*(mont_phi-1) + c_maxVal)
		if prod_acc_ubound <= prod_infNorm_max :
			return (rhoUp_log2, (nb_max_add-1)) 
	
	return;
















