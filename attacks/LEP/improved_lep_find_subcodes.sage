reset();

import numpy as np



##Lee and Brickell ISD algorithm
##Input:
# - Fq: finite field;
# - G: generator matrix
# - n: code length;
# - k: code dimension;
# - w: searched weight
# - Perm_set: set of length-n permutations
# - op_mode: 0 if you want only codewords with weight w, 1 if you want weight <= w

#The algorithm outputs three values:
# - first value: 0 if RREF is not possible, 1 if it has been computed
# - second value: 0 if no codeword has been found, 1 if ISD has found a codeword with the desired weight
# - third value: found codeword, or 0 (in case ISD has failed)

def lee_brickell_isd_Fq(Fq,G,n,k,w,Perm_set,op_mode):


    Fq_list = Fq.list();
    q = len(Fq_list);
    P = Perm_set.random_element(); #Random permutation
    Gp = G[:,P]; #Apply permutation
    Gp_0 = Gp[:,0:k]; 
    
    #If Gp_0 is singular, report failure
    if rank(Gp_0) < k:
        del Gp;
        del Gp_0;
        return 0,0,0;
    else:
        
        #Gp_0 is non singular, compute RREF
        Gp_0_inv = Gp_0^-1;
        new_G_1 = Gp_0_inv*Gp[:,k:n];
        
        flag_weight = 0; #flag_weight = 1 when we find codewor with desired weight
        
        c = vector(GF(2),n);
        
        #Go through all vectors with weight 2 and first set entry equal to 1
        #We directly compute the associated linear combination
        for i0 in range(0,k-1):
            
            g0 = new_G_1[i0,:];
            for i1 in range(i0+1,k):
                g1 = new_G_1[i1,:];
                
                #Compute all codewords in the form g0 + u*g1 (with u \neq 0)
                for i in range(1,q):
                    u = Fq_list[i];
                    cw = g0 + u*g1;
                    
                    this_w = len(cw.support())+2; #Hamming weight of the found codeword
                    if this_w <= w:
                        
                        
                        if op_mode == 0:
                            if this_w == w:
                                
                                #Compute full codeword and reverse permutation
                                full_cw = vector(Fq,n);
                                full_cw[i0] = 1;
                                full_cw[i1] = u;
                                
                                for j in range(n-k):                                    
                                    full_cw[k+j] = cw[0,j];
                                output_cw = vector(Fq,n);
                                for j in range(n):
                                    output_cw[P[j]] = full_cw[j];
                                    
                                    
                                ok = 1; #ok = 1, ISD has been successful
                                del Gp;
                                del Gp_0;
                                del Gp_0_inv;
                                del new_G_1;
                                return 1,1,output_cw;
                        else:
                            
                            #Compute full codeword and reverse permutation
                            full_cw = vector(Fq,n);
                            full_cw[i0] = 1;
                            full_cw[i1] = u;
                            
                            for j in range(n-k):
                                full_cw[k+j] = cw[0,j];
                            output_cw = vector(Fq,n);
                            for j in range(n):
                                output_cw[P[j]] = full_cw[j];

                            del Gp;
                            del Gp_0;
                            del Gp_0_inv;
                            del new_G_1;
                            return 1,1,output_cw;    
                        
        #ISD has failed: no valid codeword has been found                
        del Gp;
        del Gp_0;
        return 1,0,0;


##################################################################

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




##################################################################

#Compare lex_1 and lex_2, corresponding to c_1 and c_2 respectively, and returns the smaller one
def codewords_compare_lex(lex_1,lex_2,c_1, c_2, n):
   
   
    for i in range(n):
        #print(lex_1[i], lex_2[i]);
        #print(lex_1[i] < lex_2[i]);
        if lex_1[i] < lex_2[i]:
         #   print("c_1");
            return lex_1, c_1;
        else:
            if lex_1[i] > lex_2[i]:
          #      print("c_2");
                return lex_2, c_2;

    return lex_1, c_1;        
   
   
##################################################################

#Compute lex for vector c
def codewords_compute_lex( Fq, Fq_star, q, c, n ):
   
   
    lex_c = vector( Fq,n);
    for i in range(n):
        lex_c[i] = Fq_star[q-2];
       
   
    c_for_lex = vector(Fq,n);
   
    #Go through all scalar multiples, and return the one with smallest lexicographic value
    for u in Fq_star:
        this_c = u*c;
       
        this_lex_c = this_c.list();
        this_lex_c.sort();
       
        #Choosing which codeword has smaller lex value
        lex_c, c_for_lex = codewords_compare_lex(lex_c, this_lex_c, c_for_lex, this_c, n);
   
    return vector(c_for_lex), lex_c;    

##################################################################

import csv
load('list_sorting.sage');
load('subcodes_lex.sage');

#Code and algorithm parameters
n = 40; #code length
k = 20; #code dimension
q = 7; #finite field size
L_prime = 10; #number of codewords from each code
w_prime = 12; #codewords weight
w = 19 #subcodes support size

test_id_max = 100; # max number of generated codes; each is employed as a seed to generate the code
max_num_calls = 100; #max number of successful ISD calls

op_mode = 0; #op_mode = 0 means that we only want codewords with weight w
    
##Set finite field and group of permutations 
Fq = GF(q);
Fq_star = Set(Set(Fq)[1:q]);
Sn = Permutations(range(n));
    
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
    perm_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w) + "_" +str(w_prime) + "_" + str(test_id) + "_monomial_perm.csv";

    with open(perm_csv_name, 'w') as f1:
        writefile = csv.writer(f1)
        writefile.writerows(P_list)
    
    coeffs_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w) + "_" +str(w_prime) +  "_" + str(test_id) + "_monomial_coeffs.csv";

    with open(coeffs_csv_name, 'w') as f1:
        writefile = csv.writer(f1)
        writefile.writerows(coeffs)

    #Find subcodes from first code. Proceed until i) the number of distinct found subcodes reaches L, or ii) the total number of successful calls 
    #reaches max_num_calls (without considering different representations of the same subcode)
    num_found = 0; #number of successful ISD calls
    X_prime = []; #set with found subcodes
#    lex_X = []; #Lex values of the found subcodes
    while (len(X_prime)<L_prime)&(num_found<max_num_calls):
    
        ok_found = 0;
        while ok_found == 0:
            ok1,ok_found, c = lee_brickell_isd_Fq(Fq,G1,n,k,w_prime,Sn,op_mode);
        num_found += 1;
    
        print("Test id = "+str(test_id)+", Code 1: Num ISD successful calls = "+str(num_found)+", Num codewords = "+str(len(X_prime)));
        print("[n, k, q] = "+str([n,k,q])+", [w, w'] = "+str([w, w_prime])+", L' "+str(L_prime));
        #Avoid codewords duplicates 
        c_supp = vector(c.support());        
        c = c[c_supp[0]]^(-1)*c;
        
      
        #See if the codeword is a new one, and eventually put it into X
        collisions = colliding_elements(X_prime,[c]);
        if len(collisions) == 0:
            X_prime.append(c);
#            lex_c = codewords_compute_lex(Fq,Fq_star,n,vector(c[0,:]),vector(c[1,:]),w);
#            lex_X.append(lex_c);
    
    
    #Write number of found codewords into csv file
#    L_prime_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w) + "_" +str(w_prime) + "_" + str(test_id) + "_L_prime_1.csv";
    L_prime_values = [];
    L_prime_values.append([len(X_prime)]);
#    a.append(len(X_prime));
#    with open(L_prime_csv_name, 'w') as f1:
#        writefile = csv.writer(f1)
#        writefile.writerows(a)
    
    #Combine pairs of the found codewords to produce subcodes with support size w
    
    X = []; #set with found subcodes
    lex_X = []; #Lex values of the found subcodes
    
    for i in range(len(X_prime)-1):
        a = X_prime[i];
        a_supp = vector(a).support();
        for j in range(i+1,len(X_prime)):
            b = X_prime[j];
            b_supp = vector(b).support();
            
            #Check the support size of subocode generated by {a ; b}
            c_supp = Set(a_supp).union(Set(b_supp));
            if len(c_supp) == w:
                c = matrix(Fq,[a,b]);
                c = c.echelon_form();
        
      
                #See if the subcode is a new one, and eventually put it into X
                
                collisions = colliding_elements(X,[c]);
                if len(collisions) == 0:
                    X.append(c);
                    lex_c = compute_lex(Fq,Fq_star,n,vector(c[0,:]),vector(c[1,:]));
                    lex_X.append(lex_c);
            
    #Write found subcodes into csv file
    
    subcodes_1_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w)+"_" +str(w_prime) + "_" + str(test_id) + "_list_subcodes_1.csv";
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
    lex_1_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w)+"_"  +str(w_prime) + "_" +  str(test_id) + "_list_lex_subcodes_1.csv";
    list_for_csv = [];
    for i in range(len(lex_X)):
        a = lex_X[i][0,:];
        b = lex_X[i][1,:];
        list_for_csv.append(vector(a));
        list_for_csv.append(vector(b));
    
    with open(lex_1_csv_name, 'w') as f1:
        writefile = csv.writer(f1)
        writefile.writerows(list_for_csv)  
        
    ############################################################    
    ######Repeat the whole thing for the second code    
    
    num_found = 0; #number of successful ISD calls
    Y_prime = []; #set with found subcodes
#    lex_X = []; #Lex values of the found subcodes
    while (len(Y_prime)<L_prime)&(num_found<max_num_calls):
    
        ok_found = 0;
        while ok_found == 0:
            ok1,ok_found, c = lee_brickell_isd_Fq(Fq,G2,n,k,w_prime,Sn,op_mode);
        num_found += 1;
    
        print("Test id = "+str(test_id)+", Code 2: Num ISD successful calls = "+str(num_found)+", Num codewords = "+str(len(Y_prime)));
        print("[n, k, q] = "+str([n,k,q])+", [w, w'] = "+str([w, w_prime])+", L' "+str(L_prime));
    
        #Avoid codewords duplicates 
        c_supp = vector(c.support());        
        c = c[c_supp[0]]^(-1)*c;
        
      
        #See if the codeword is a new one, and eventually put it into X
        collisions = colliding_elements(Y_prime,[c]);
        if len(collisions) == 0:
            Y_prime.append(c);
#            lex_c = codewords_compute_lex(Fq,Fq_star,n,vector(c[0,:]),vector(c[1,:]),w);
#            lex_X.append(lex_c);
    
    #Write number of found codewords into csv file
    L_prime_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w) + "_" +str(w_prime) + "_" + str(test_id) + "_L_prime.csv";
    L_prime_values.append([len(Y_prime)]);
    with open(L_prime_csv_name, 'w') as f1:
        writefile = csv.writer(f1)
        writefile.writerows(L_prime_values)
    
    #Combine pairs of the found codewords to produce subcodes with support size w
    
    Y = []; #set with found subcodes
    lex_Y = []; #Lex values of the found subcodes
    
    for i in range(len(Y_prime)-1):
        a = Y_prime[i];
        a_supp = vector(a).support();
        for j in range(i+1,len(Y_prime)):
            b = Y_prime[j];
            b_supp = vector(b).support();
            
            #Check the support size of subocode generated by {a ; b}
            c_supp = Set(a_supp).union(Set(b_supp));
            if len(c_supp) == w:
                c = matrix(Fq,[a,b]);
                c = c.echelon_form();
        
      
                #See if the subcode is a new one, and eventually put it into X
                
                collisions = colliding_elements(Y,[c]);
                if len(collisions) == 0:
                    Y.append(c);
                    lex_c = compute_lex(Fq,Fq_star,n,vector(c[0,:]),vector(c[1,:]));
                    lex_Y.append(lex_c);
            
            
    #Write found subcodes into csv file
    
    subcodes_2_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w)+"_" +str(w_prime) + "_" + str(test_id) + "_list_subcodes_2.csv";
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
    lex_2_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w)+"_"  +str(w_prime) + "_" +  str(test_id) + "_list_lex_subcodes_2.csv";
    list_for_csv = [];
    for i in range(len(lex_Y)):
        a = lex_Y[i][0,:];
        b = lex_Y[i][1,:];
        list_for_csv.append(vector(a));
        list_for_csv.append(vector(b));
    
    with open(lex_2_csv_name, 'w') as f1:
        writefile = csv.writer(f1)
        writefile.writerows(list_for_csv)  
