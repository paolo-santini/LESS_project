reset();

import numpy as np



##Lee and Brickell ISD algorithm to find subcodes
##Input:
# - Fq: finite field;
# - G: generator matrix
# - n: code length;
# - k: code dimension;
# - w: searched weight
# - Perm_set: set of length-n permutations
# - op_mode: 0 if you want only subcodes with support size w, 1 if you want support size <= w

#The algorithm outputs three values:
# - first value: 0 if RREF is not possible, 1 if it has been computed
# - second value: 0 if no subcode has been found, 1 if ISD has found a subcode with the desired support size
# - third value: found subcode, or 0 (in case ISD has failed)

def subcodes_lee_brickell_isd(Fq,G,n,k,w,Perm_set,op_mode):

   
    ##Pick a random permutation
    P = Perm_set.random_element();
    Gp = G[:,P];
    Gp_0 = Gp[:,0:k];
   
    #Verify rank of left submatrix
    if rank(Gp_0) < k:
        del Gp;
        del Gp_0;
        return 0,0,0;
    else:
        #Go on and put matrix in systematic form
        Gp_0_inv = Gp_0^-1;
        new_G = Gp_0_inv*Gp;
        flag_weight = 0;    
        for i0 in range(0,k-1):
            g0 = new_G[i0,:];
            supp0 = vector(g0).support();
            for i1 in range(i0+1,k):
                g1 = new_G[i1,:];
                supp1 = vector(g1).support()
                this_w = len(Set(supp0).union(Set(supp1)));
               
                if this_w <= w:
                    if op_mode == 0:
                        if this_w == w:
                            full_cw = matrix(Fq,2,n);

                            for j in range(n):                                    
                                full_cw[0,j] = g0[0,j];
                                full_cw[1,j] = g1[0,j];
                            output_cw = matrix(Fq,2,n);
                            for j in range(n):
                                output_cw[0,P[j]] = full_cw[0,j];
                                output_cw[1,P[j]] = full_cw[1,j];
                            ok = 1;
                            del Gp;
                            del Gp_0;
                            del Gp_0_inv;
                            del new_G;
                            return 1,1,output_cw;
                        
                    else:
                        
                        full_cw = matrix(Fq,2,n);

                        for j in range(n):                                    
                            full_cw[0,j] = g0[0,j];
                            full_cw[1,j] = g1[0,j];
                        output_cw = matrix(Fq,2,n);
                        for j in range(n):
                            output_cw[0,P[j]] = full_cw[0,j];
                            output_cw[1,P[j]] = full_cw[1,j];

                        del Gp;
                        del Gp_0;
                        del Gp_0_inv;
                        del new_G;
                        return 1,1,output_cw;    
                    
        del Gp;
        del Gp_0;
        return 1,0,0;
    
########################################################################   

#Generate random monomial matrix
def random_monomial(Fq,Fq_star,Sn,n):

    P = Sn.random_element();
    Q = matrix(Fq,n,n);
    coeffs = [];
    P_list = [];
    for i in range(n):
        a = Fq_star.random_element();
        Q[i,P[i]] = a;
        coeffs.append([a]);
        P_list.append([P[i]]);
        
    return P_list, coeffs, Q;



#########################################################

import csv
load('list_sorting.sage');
load('subcodes_lex.sage');

#Code and algorithm parameters
n = 30; #code length
k = 20; #code dimension
q = 13; #finite field size
w = 10; #codewords weight
L = 100; #number of codewords from each code

test_id_max = 150; # max number of generated codes; each is employed as a seed to generate the code
max_num_calls = 100; #max number of successful ISD calls

op_mode = 0; #op_mode = 0 means that we only want codewords with weight w
    
##Set finite field and group of permutations 
Fq = GF(q);
Fq_star = Set(Set(Fq)[1:q]);
Sn = Permutations(range(n));
    
print("Expected number of subcodes is "+str((q^2-1)^w*binomial(n,w)/( (q^2-q) *(q^2-1) )*q^(-2.*(n-k))));    
#Start generating codes and launching ISD to find codewords    
for test_id in range(test_id_max):
    
    #Use test_id as the seed for the simulation
    set_random_seed(test_id);

    #Generate random code with dimension k
    r = 0;
    while r<k:
        G1 = random_matrix(Fq,k,n);
        r = rank(G1);
    

    #Extract random monomial and transform the code
    P_list,coeffs,Q = random_monomial(Fq,Fq_star,Sn,n);    
    G2 = G1*Q;

    
    #Write monomial into csv file
    perm_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w) + "_" + str(test_id) + "_monomial_perm.csv";

    with open(perm_csv_name, 'w') as f1:
        writefile = csv.writer(f1)
        writefile.writerows(P_list)
    
    coeffs_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w) + "_" + str(test_id) + "_monomial_coeffs.csv";

    with open(coeffs_csv_name, 'w') as f1:
        writefile = csv.writer(f1)
        writefile.writerows(coeffs)

    #Find subcodes from first code. Proceed until i) the number of distinct found subcodes reaches L, or ii) the total number of successful calls 
    #reaches max_num_calls (without considering different representations of the same subcode)
    num_found = 0; #number of successful ISD calls
    X = []; #set with found subcodes
    lex_X = []; #Lex values of the found subcodes
    while (len(X)<L)&(num_found<max_num_calls):
    
        ok_found = 0;
        while ok_found == 0:
            ok1,ok_found, c = subcodes_lee_brickell_isd(Fq,G1,n,k,w,Sn,op_mode);
        num_found += 1;
    
        print("Test id = "+str(test_id)+", Code 1: Num ISD successful calls = "+str(num_found)+", Num subcodes = "+str(len(X)));
    
        #Use echelon form to avoid duplicates 
        c = c.echelon_form();
        
      
        #See if the codeword is a new one, and eventually put it into X
        collisions = colliding_elements(X,[c]);
        if len(collisions) == 0:
            X.append(c);
            lex_c = compute_lex(Fq,Fq_star,n,vector(c[0,:]),vector(c[1,:]));
            lex_X.append(lex_c);


    #Write found subcodes into csv file
    
    subcodes_1_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w)+"_" + str(test_id) + "_list_subcodes_1.csv";
    list_for_csv = [];
    for i in range(len(X)):
        a = X[i][0,:];
        b = X[i][1,:];
        list_for_csv.append(vector(a));
        list_for_csv.append(vector(b));
        
    with open(subcodes_1_csv_name, 'w') as f1:
        writefile = csv.writer(f1)
        writefile.writerows(list_for_csv)
    
    #Write lex values into csv file
    lex_1_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w)+"_" + str(test_id) + "_list_lex_subcodes_1.csv";
    list_for_csv = [];
    for i in range(len(lex_X)):
        a = lex_X[i][0,:];
        b = lex_X[i][1,:];
        list_for_csv.append(vector(a));
        list_for_csv.append(vector(b));
    
    with open(lex_1_csv_name, 'w') as f1:
        writefile = csv.writer(f1)
        writefile.writerows(list_for_csv)    
    
    
    #Repeat the whole thing for second code
    
    
    num_found = 0; #number of successful ISD calls
    Y = []; #set with found subcodes
    lex_Y = []; #Lex values of the found subcodes
    while (len(Y)<L)&(num_found<max_num_calls):
    
        ok_found = 0;
        while ok_found == 0:
            ok1,ok_found, c = subcodes_lee_brickell_isd(Fq,G2,n,k,w,Sn,op_mode);
        num_found += 1;
    
        print("Test id = "+str(test_id)+", Code 2: Num ISD successful calls = "+str(num_found)+", Num subcodes = "+str(len(Y)));
    
        #Use echelon form to avoid duplicates 
        c = c.echelon_form();
        
      
        #See if the codeword is a new one, and eventually put it into X
        collisions = colliding_elements(Y,[c]);
        if len(collisions) == 0:
            Y.append(c);
            lex_c = compute_lex(Fq,Fq_star,n,vector(c[0,:]),vector(c[1,:]));
            lex_Y.append(lex_c);


    #Write found subcodes into csv file
    
    subcodes_2_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w)+"_" + str(test_id) + "_list_subcodes_2.csv";
    list_for_csv = [];
    for i in range(len(Y)):
        a = Y[i][0,:];
        b = Y[i][1,:];
        list_for_csv.append(vector(a));
        list_for_csv.append(vector(b));
        
    with open(subcodes_2_csv_name, 'w') as f1:
        writefile = csv.writer(f1)
        writefile.writerows(list_for_csv)
    
    #Write lex values into csv file
    lex_2_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w)+"_" + str(test_id) + "_list_lex_subcodes_2.csv";
    list_for_csv = [];
    for i in range(len(lex_Y)):
        a = lex_Y[i][0,:];
        b = lex_Y[i][1,:];
        list_for_csv.append(vector(a));
        list_for_csv.append(vector(b));
    
    with open(lex_2_csv_name, 'w') as f1:
        writefile = csv.writer(f1)
        writefile.writerows(list_for_csv)
    
    print("We have done finding codewords...");
    

