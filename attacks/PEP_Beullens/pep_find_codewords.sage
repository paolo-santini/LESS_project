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
    for i in range(n):
        Q[i,P[i]] = Fq_star.random_element();
   
    return Q;

##################################################################

#Compare lex_1 and lex_2, corresponding to c_1 and c_2 respectively, and returns the smaller one
def compare_lex(lex_1,lex_2,c_1, c_2, n):
   
   
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
def compute_lex( Fq, Fq_star, q, c, n ):
   
   
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
        lex_c, c_for_lex = compare_lex(lex_c, this_lex_c, c_for_lex, this_c, n);
   
    return vector(c_for_lex), lex_c;    

##################################################################

import csv
load('list_sorting.sage');


#Code and algorithm parameters
n = 50; #code length
k = 25; #code dimension
q = 5; #finite field size
w = 13; #codewords weight
L = 100; #number of codewords from each code

max_num_calls = 25; #max number of successful ISD calls
test_id_max = 10; # max number of generated codes; each is employed as a seed to generate the code

op_mode = 0; #op_mode = 0 means that we only want codewords with weight w
    
##Set finite field and group of permutations 
Fq = GF(q);
Fq_star = Set(Set(Fq)[1:q]);
Sn = Permutations(range(n));

Nw = binomial(n,w)*(q-1)^(w-1)*(q^k-1.)/(q^n-1.);
print("Expected number of codewords is " +str(Nw));

#Start generating codes and launching ISD to find codewords    
for test_id in range(test_id_max):
    
    #Use test_id as the seed for the simulation
    set_random_seed(test_id);

    #Generate random code with dimension k
    r = 0;
    while r<k:
        G1 = random_matrix(Fq,k,n);
        r = rank(G1);
    

    #Extract random permutation and permute the code
    P = Sn.random_element();
    P_list = [];
    for i in range(n):
        P_list.append([P[i]]);
    
    Q = matrix(Fq,n,n);
    for i in range(n):
        Q[i,P[i]] = Fq(1);
    
    G2 = G1*Q;

    
    #Write permutation into csv file
    perm_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w) + "_" + str(test_id) + "_permutation.csv";

    with open(perm_csv_name, 'w') as f1:
        writefile = csv.writer(f1)
        writefile.writerows(P_list)


    #Find codewords from first code. Proceed until i) the number of distinct found codewords reaches L, or ii) the total number of successful calls 
    #reaches max_num_calls (without considering scalar multiples of the same codeword)
    num_found = 0; #number of successful ISD calls
    X = []; #set with found codewords
    lex_X = []; #Lex values of the found codewords
    while (len(X)<L)&(num_found<max_num_calls):
    
        ok_found = 0;
        while ok_found == 0:
            ok1,ok_found, c = lee_brickell_isd_Fq(Fq,G1,n,k,w,Sn,op_mode);
        num_found += 1;
    
        print("Test id = "+str(test_id)+", Code 1: Num ISD successful calls = "+str(num_found)+", Num codewords = "+str(len(X)));
    

        c_for_lex, lex_c = compute_lex( Fq, Fq_star, q, c, n ); #compute lex
    
        #See if the codeword is a new one, and eventually put it into X
        collisions = colliding_elements(X,[c_for_lex]);
        if len(collisions) == 0:
            X.append(c_for_lex);
            lex_X.append(lex_c);


    #Write found codewords into csv file
    codewords_1_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w)+"_" + str(test_id) + "_list_codewords_1.csv";
    with open(codewords_1_csv_name, 'w') as f1:
        writefile = csv.writer(f1)
        writefile.writerows(X)
    
    #Write lex values into csv file
    lex_1_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w)+"_" + str(test_id) + "_list_lex_1.csv";
    with open(lex_1_csv_name, 'w') as f1:
        writefile = csv.writer(f1)
        writefile.writerows(lex_X)    
    
    
    
    #######Repeat the whole thing on the second code
    #Finding codewords from first code
    num_found = 0;
    Y = [];
    lex_Y = [];
    while (len(Y)<L)&(num_found<max_num_calls):
    
        ok_found = 0;
    
        while ok_found == 0:
            ok1,ok_found, c = lee_brickell_isd_Fq(Fq,G2,n,k,w,Sn,op_mode);
        num_found += 1;
    
        print("Test id = "+str(test_id)+", Code 2: Num ISD successful calls = "+str(num_found)+", Num codewords = "+str(len(Y)));
    
        c_for_lex, lex_c = compute_lex( Fq, Fq_star, q, c, n ); #compute lex
     
        #See if the codeword is a new one, and eventually put it into Y
        collisions = colliding_elements(Y,[c_for_lex]);
        if len(collisions) == 0:
            Y.append(c_for_lex);
            lex_Y.append(lex_c);


    #Write codewords into csv file
    codewords_1_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w)+"_" + str(test_id) + "_list_codewords_2.csv";
    with open(codewords_1_csv_name, 'w') as f1:
        writefile = csv.writer(f1)
        writefile.writerows(Y)
    
    #Write permutation into csv file
    lex_1_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w)+"_" + str(test_id) + "_list_lex_2.csv";
    with open(lex_1_csv_name, 'w') as f1:
        writefile = csv.writer(f1)
        writefile.writerows(lex_Y)   
    

    print("We have done finding codewords...");

