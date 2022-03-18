reset();

load('list_sorting.sage');

import csv

#Compute the gaussian binomial coefficient

def gauss_binomial(m,r,q):
    val = 1;
    for i in range(r):
        val = val*(1-q^(m-i))/(1-q^(i+1));
    return val;

##########################################

############Test of the monomial recovery attack

##This function tests the monomial recovery strategy
#For the sake of implicity, we only test the recovery of the permutation

#Input values:
# - Fq: finite field;
# - L: number of codewords for each code;
# - list_lex_1, list_lex_2: lists containing the lex values of the found subcodes
# - list_1, list_2: lists containing the found subcodes
# - Q: monomial

#The function returns:
# - the number of good collisions;
# - the number of bad collisions;
# - the estimate on P (the permutation part of Q)

def monomial_recovery(Fq, L, list_lex_1, list_lex_2, list_1, list_2 , Q):
   
    ###Find list of matching subcodes
    indexes = colliding_indexes(list_lex_1[0:L],list_lex_2[0:L]); #indexes of colliding subcodes
    
    list_matching_subcodes = []; #list with colliding codewords
    num_wrong = 0; #number of bad collisions
  #  print(len(indexes));
    for i in range(len(indexes)):
        
        c1 = list_1[indexes[i][0]];
        c2 = list_2[indexes[i][1]];
        
        list_matching_subcodes.append([c1,c2]); #append new pair of colliding codewords
        
        
        #Verify if c2 == C1*Q
        c1_mono = matrix(c1*Q).echelon_form();        
        if c1_mono != matrix(c2).echelon_form():
            num_wrong += 1;

    ###First, reconstruct the permutation part of Q
    my_P = ones_matrix(Fq,n);
    
    for i in range(len(list_matching_subcodes)):
        c1 = list_matching_subcodes[i][0];
        c2 = list_matching_subcodes[i][1];
        for i0 in range(n):
            val1 = (c1[0,i0]==0)&(c1[1,i0]==0);
            for i1 in range(n):
                val2 = (c2[0,i1]==0)&(c2[1,i1]==0);
                if (val1!=val2):
                    my_P[i0,i1] = 0;
    
       
    return len(list_matching_subcodes), num_wrong, my_P

#######################################################################################

#Arrange data from csv into list of two-dimensional matrices
def list_from_cvs(Fq,n,file_name):
    
    ##Read raw data in csv file
    with open(file_name,'rU') as f1:
        raw_list=list( csv.reader(f1) )
    

    #Build list of two-dimensional subcodes from raw list

    parsed_list = [];

    for i in range(0, len(raw_list),  2):
        
        c = matrix(Fq,2,n);
        this_a = raw_list[i];
        this_b = raw_list[i+1];
        for j in range(n):
            c[0,j] = Fq( this_a[j] );
            c[1,j] = Fq( this_b[j] );
        
        parsed_list.append( c );
        
    return parsed_list;
#######################################################################################



#Code and algorithm parameters
n = 30; #code length
k = 20; #code dimension
q = 13; #finite field size
w = 10; #codewords weight
L = 100; #number of codewords from each code

test_id_max = 100; # max number of generated codes; each is employed as a seed to generate the code

#Prepare to analyze data
test_id_vec = list(range(test_id_max));
tensor_L_values = [];
tensor_correct_ratios = [];
tensor_num_matching_list = [];
tensor_num_wrong_list = [];

max_L = 10000000000000000000000000000000000000;

for test_id in test_id_vec:
    
    print("Doing seed = "+str(test_id));
    
    #Names of csv files we need for the simulation
    perm_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w) + "_" + str(test_id) + "_monomial_perm.csv";
    coeffs_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w) + "_" + str(test_id) + "_monomial_coeffs.csv";          
    subcodes_1_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w)+"_" + str(test_id) + "_list_subcodes_1.csv";
    lex_1_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w)+"_" + str(test_id) + "_list_lex_subcodes_1.csv";
    subcodes_2_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w)+"_" + str(test_id) + "_list_subcodes_2.csv";
    lex_2_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w)+"_" + str(test_id) + "_list_lex_subcodes_2.csv";

    
    #Read and build monomial
    Fq = GF(q);

    #Read and build permutation
    with open(perm_csv_name,'rU') as f1:
        P=list( csv.reader(f1) )
        
    #Read scaling coefficients
    with open(coeffs_csv_name,'rU') as f1:
        coeffs=list( csv.reader(f1) )
        
    #Build monomial
    secret_permutation = [];
    secret_coeffs = [];
    for i in range(n):
        secret_permutation.append(ZZ(P[i][0]));
        secret_coeffs.append(coeffs[i][0]);

    Q = matrix(Fq,n,n);
    secret_P_matrix = matrix(Fq,n);
    for i in range(n):
        secret_P_matrix[i,secret_permutation[i]] = Fq(1);
        Q[i,secret_permutation[i]] = secret_coeffs[i];    

    
    #Arrange csv files into proper lists of two-dimensional subcodes
    found_subcodes_1 = list_from_cvs(Fq,n,subcodes_1_csv_name);
    found_subcodes_2 = list_from_cvs(Fq,n,subcodes_2_csv_name);
    found_lex_1 = list_from_cvs(Fq,n,lex_1_csv_name);
    found_lex_2 = list_from_cvs(Fq,n,lex_2_csv_name);
    
    

    ######################################################################################
    L_values = [];
    correct_ratios = [];
    num_matching_list = [];
    num_wrong_list = [];
    
    this_L = min(len(found_lex_1), len(found_lex_2)); #maximum value of L for the current pair of codes
    if this_L < max_L:
        max_L = this_L;
    
    #Test the attack for different values of L
    for L in range(1,min(len(found_lex_1), len(found_lex_2))+1 ):    
        
        
        num_matching, num_wrong, my_P = monomial_recovery( Fq, L, found_lex_1, found_lex_2, found_subcodes_1, found_subcodes_2, Q );    
        
        #pr_correct = 1 if found permutation is correct, otherwise pr_correct = 0
        if my_P == secret_P_matrix:
            pr_correct = 1;
        else:
            pr_correct = 0;
        
        #Update data
        L_values.append(L);
        correct_ratios.append(pr_correct);
        num_matching_list.append(num_matching);
        num_wrong_list.append(num_wrong);
        
        
    #Saving results
    tensor_L_values.append(L_values);
    tensor_correct_ratios.append(correct_ratios); 
    tensor_num_matching_list.append(num_matching_list);
    tensor_num_wrong_list.append(num_wrong_list);
    

#Analyze found results: average over the results of all tests

emp_M_prime_list = [];
emp_M_second_list = [];

for i in range(max_L):
    L = i+1;
    emp_M_prime = 0;
    emp_M_second = 0;
    for j in range(test_id_max):
        num_matching_list = tensor_num_matching_list[j];
        num_wrong_list = tensor_num_wrong_list[j];
        emp_M_prime += num_matching_list[i]-num_wrong_list[i];
        emp_M_second += num_wrong_list[i];
    
    emp_M_prime_list.append(emp_M_prime/test_id_max*1.);
    emp_M_second_list.append(emp_M_second/test_id_max*1.);
    
####################################      

Nw = binomial(n,w)*((q^2-1)^w)/((q^2-1)*(q^2-q))*gauss_binomial(k,2,q)/gauss_binomial(n,2,q);    
print("Values of M' as a function of L");    
print("Format: (L, emp M')");
for i in range(0,len(emp_M_prime_list),1):
#    L = i+1;
#    th_M_prime = L^2/Nw;
    print("("+str(i+1)+" , "+str(emp_M_prime_list[i]*1.)+")");
#    print("("+str(i+1)+", "+str(emp_M_prime_list[i]*1.)+")");

print("Values of M' as a function of L");    
print("Format: (L, th M')");
for i in range(0,len(emp_M_prime_list),1):
    L = i+1;
    th_M_prime = L^2/Nw;
    print("("+str(i+1)+" , "+str(th_M_prime*1.)+")");
#    print("("+str(i+1)+", "+str(emp_M_prime_list[i]*1.)+")");


print("Values of M'' as a function of L");    
print("Format: [L, emp M'', th M'']");
tw = binomial(n,w)*factorial(w)*(q-1)^(w-1)*gaussian_binomial(k,2,q)/gaussian_binomial(n,2,q);    

#pr_same = (1/(q^2-1))^(w-1);




for i in range(0,len(emp_M_second_list),1):
    L = i+1;
    th_M_prime = L^2/Nw;
    th_M_second = tw*(L^2-th_M_prime)/Nw;
    print("("+str(i+1)+" , "+str(emp_M_second_list[i]*1.)+" , "+str(th_M_second*1.)+")");
    
