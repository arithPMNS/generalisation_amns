# AMNS
This repository contains codes to generate Adapted Modular Number Systems (AMNS) for a given prime.
It also contains a C code generator for arithmetic (and conversion) operations using AMNS.


The subdirectory 'amns_generator' contains codes to generate AMNS given a prime and some parameters; see an example in file 'gen_example.py' of this subdirectory. NOTE: Here, AMNS generation is done according to the implementation strategy mentioned at the end of section 5.1 of the article.


The subdirectory 'c_code_generator' contains codes to generate a C code from an AMNS generated using codes in the subdirectory 'amns_generator'. This C code allows to perform arithmetic (and conversion) operations efficiently; see an example in file 'EXAMPLE.py' of this subdirectory.          
The C code generated will be in the subdirectory 'c_codes' of this subdirectory and its entry point is the file 'main.c' which contains the commands to compile and execute it. Information about the AMNS can be found in the file 'amns_init.c' (the value of the prime and the parameter 'gamma' are written in the function 'init_data'), the file 'structs_data.h' (contains, among others, the value of 'n' and 'rho') and the file 'add_mult_poly.c' (contains functions that directly use the polynomials M and M' for operations).


The subdirectory 'amns_for_our_tests' contains the C codes of the AMNS we used to build the tables of performances in the article.
In this subdirectory, files are named as follow : 'p'+pBitSize_num1__n_lambda__num2; n and lambda are some parameters of the AMNS. num1 and num2 are irrelevant numbers we used to distinguish AMNS in the generation process.
Notice that in each subdirectory of this subdirectory, there is a file 'timing.txt' which contains results we used to build the tables of performances in the article. This file shows computation time of modular multiplications using AMNS, GNU MP mpz_t and OPENSSL BIGNUM.
The subdirectory 'amns_vs_gmpLowLevel' is identical to the subdirectory 'amns_for_our_tests' except that it compares AMNS to GNU MP low level functions for modular multiplications. GNU MP low level functions are faster than the high level functions provided when using the mpz_t type.
Also, these subdirectories contain more than 1600 AMNS, so GitHub will not show all of them. Therefore, you should download these subdirectories to see all these AMNS.



Note : To use the AMNS generator and the C code generator, you will need SageMath library which can be found here: http://www.sagemath.org/

