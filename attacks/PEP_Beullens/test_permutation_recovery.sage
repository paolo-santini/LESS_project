reset();

load('list_sorting.sage');
import csv


############Test of the permutation recovery attack

##This function tests the permutation recovery strategy 
#Input values:
# - Fq: finite field;
# - L: number of codewords for each code;
# - list_lex_1, list_lex_2: lists containing the lex values of the found codewords
# - list_1, list_2: lists containing the found codewords
# - Q: permutation

#The function returns:
# - the number of good collisions;
# - the number of bad collisions;
# - the estimate on Q

def permutation_recovery(Fq, L, list_lex_1, list_lex_2, list_1, list_2 , Q):
   
    ###Find list of matching codewords
    indexes = colliding_indexes(list_lex_1[0:L],list_lex_2[0:L]); #indexes of colliding codewords
    
    list_matching_codewords = []; #list with colliding codewords
    num_wrong = 0; #number of bad collisions

    for i in range(len(indexes)):
        
        c1 = list_1[indexes[i][0]];
        c2 = list_2[indexes[i][1]];
        
        list_matching_codewords.append([c1,c2]); #append new pair of colliding codewords
        
        
        #Verify if c2 == C1*Q
        c1_perm = matrix(c1*Q);        
        if c1_perm != matrix(c2):
            num_wrong += 1;

    ###Reconstruct permutation matrix
    my_Q = ones_matrix(Fq,n);

    for i in range(len(list_matching_codewords)):
        c1 = list_matching_codewords[i][0];
        c2 = list_matching_codewords[i][1];
        for x in range(n):
            val_1 = c1[x];
            for y in range(n):
                if c2[y]!= val_1:
                    my_Q[x,y] = 0;
       
    return len(list_matching_codewords), num_wrong, my_Q

#######################################################################################
#######################################################################################


#Code and algorithm parameters
n = 50; #code length
k = 25; #code dimension
q = 5; #finite field size
w = 13; #codewords weight
L = 100; #number of codewords from each code

max_num_calls = 150; #max number of successful ISD calls
test_id_max = 20; # max number of generated codes; each is employed as a seed to generate the code

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
    perm_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w)+ "_" + str(test_id) + "_permutation.csv";
    codewords_1_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w)+ "_" + str(test_id) + "_list_codewords_1.csv";
    lex_1_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w)+ "_" + str(test_id) + "_list_lex_1.csv";
    codewords_2_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w)+ "_" + str(test_id) + "_list_codewords_2.csv";
    lex_2_csv_name = "CodeEquivalence_results/csv_"+str(q)+"_"+str(n)+"_"+str(k)+"_"+str(w)+ "_" + str(test_id) + "_list_lex_2.csv";



    Fq = GF(q);

    #Read and build permutation
    with open(perm_csv_name,'rU') as f1:
        P=list( csv.reader(f1) )
    
    secret_permutation = [];
    for i in range(n):
        secret_permutation.append(ZZ(P[i][0]));

    Q = matrix(ZZ,n,n);
    for i in range(n):
        Q[i,secret_permutation[i]] = 1;    


    ##Read codewords and lex values
    with open(codewords_1_csv_name,'rU') as f1:
        list_codewords_1=list( csv.reader(f1) )

    with open(lex_1_csv_name,'rU') as f1:
        list_lex_1=list( csv.reader(f1) )    
    
    with open(codewords_2_csv_name,'rU') as f1:
        list_codewords_2 = list( csv.reader(f1) )
    
    with open(lex_2_csv_name,'rU') as f1:
        list_lex_2 = list( csv.reader(f1) )  
    

    #Build lists with found codewords and lex values

    ##Code 1
    found_codewords_1 = [];
    found_lex_1 = [];

    for i in range( len(list_lex_1) ):
        
        c = vector(ZZ,n);
        this_c = list_codewords_1[i];
        for j in range(n):
            c[j] = ZZ( this_c[j] );
        
        found_codewords_1.append( c );
    
        c = vector(ZZ,n);
        this_c = list_lex_1[i];
        for j in range(n):
            c[j] = ZZ( this_c[j] );
        
        found_lex_1.append( c );
    
    ##Code 2
    found_codewords_2 = [];
    found_lex_2 = [];

    for i in range( len(list_lex_2) ):
    
        c = vector(ZZ,n);
        this_c = list_codewords_2[i];
        for j in range(n):
            c[j] = ZZ( this_c[j] );
        
        found_codewords_2.append( c );
    
    
        c = vector(ZZ,n);
        this_c = list_lex_2[i];
        for j in range(n):
            c[j] = ZZ( this_c[j] );
        
        found_lex_2.append( c );    
    

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
        
        
        num_matching, num_wrong, my_Q = permutation_recovery( Fq, L, found_lex_1, found_lex_2, found_codewords_1, found_codewords_2, Q );    
        
        #pr_correct = 1 if found permutation is correct, otherwise pr_correct = 0
        if my_Q == Q:
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
    
    #Displaying results
#    print("\addplot[dotdash, black, mark size=2pt,line width=1pt]coordinates{");
#    for j in range(len(L_values)):
#        print("("+str(L_values[j])+", "+str(correct_ratios[j])+")");
#    print("};");    


#Analyze found results

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
print("Values of M' as a function of L");    
for i in range(0,len(emp_M_prime_list),1):
    print("("+str(i+1)+", "+str(emp_M_prime_list[i]*1.)+")");
    
print("             ");        
print("Values of M'' as a function of L");    
for i in range(0,len(emp_M_prime_list),1):
    print("("+str(i+1)+", "+str(emp_M_second_list[i]*1.)+")");    
    
print("             ");        
print("Success probability as a function of L");        
pr_correct = vector(RR,max_L);
for L in range(max_L):
    pr = 0;
    for id in range(test_id_max):
         pr += tensor_correct_ratios[id][L];
    pr_correct[L] = pr/test_id_max*1.;
    
for i in range(max_L):
    print("("+str(i+1)+", "+str(pr_correct[i]*1.)+")");
