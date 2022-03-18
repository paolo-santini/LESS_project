##################################################################
import numpy;

def compare_matrix_lex(G1,G2,n):
   
    w1 = n-G1[0,:].list().count(0);
    w2 = n-G2[0,:].list().count(0);
    if w1 < w2:
        return G1;
    if w1>w2:
        return G2;
   
    w1 = n-G1[1,:].list().count(0);
    w2 = n-G2[1,:].list().count(0);
    if w1 < w2:
        return G1;
    if w1>w2:
        return G2;
   
    #If arrived up to this point, it means that the matrices have the same row weights
    #We consider the second row of each matrix
    for j in range(n):
        if G1[1,j]<G2[1,j]:
            return G1;
        else:
            if G1[1,j]>G2[1,j]:
                return G2;
    return G1;

#######################################


def minimal_matrix(Fq,n,a,b):

    ###reducing all elements in the first row    
    lex_G = matrix([a,b]);
    first_row_support = vector(lex_G[0,:]).support();
    for i in first_row_support:
        lex_G[1,i] = lex_G[0,i]^-1*lex_G[1,i];
        lex_G[0,i] = Fq(1);

    #Find null columns
    null_columns_indexes = [];
    seminull_columns_indexes = [];
    other_columns_indexes = [];
    for i in range(n):
        if lex_G[0,i] == 0:
            if lex_G[1,i] == 0:
                null_columns_indexes.append(i);
            else:
                seminull_columns_indexes.append(i);
        else:
            other_columns_indexes.append(i);


    new_lex_G = matrix(Fq,2,n);
    for i in range(len(null_columns_indexes)):
        new_lex_G[0,i] = lex_G[0,null_columns_indexes[i]];
        new_lex_G[1,i] = lex_G[1,null_columns_indexes[i]];


    for i in range(len(seminull_columns_indexes)):
        new_lex_G[0,len(null_columns_indexes)+i] = lex_G[0,seminull_columns_indexes[i]];
        new_lex_G[1,len(null_columns_indexes)+i] = Fq(1);

    for i in range(len(other_columns_indexes)):
        new_lex_G[0,len(null_columns_indexes)+len(seminull_columns_indexes)+i] = lex_G[0,other_columns_indexes[i]];
        new_lex_G[1,len(null_columns_indexes)+len(seminull_columns_indexes)+i] = lex_G[1,other_columns_indexes[i]];    
#    print(new_lex_G);
#    print(" = ");
    other_columns_indexes = range(len(null_columns_indexes)+len(seminull_columns_indexes),n);

    ##Apply permutation to find
    second_row_elements = new_lex_G[1,other_columns_indexes];
  #  print(second_row_elements);
    P = numpy.argsort(second_row_elements);
    second_new_lex_G = copy(new_lex_G);
    for j in range( len(other_columns_indexes ) ):
        second_new_lex_G[1,other_columns_indexes[j]] = second_row_elements[0,P[0,j]];
  #  print(second_new_lex_G);
    return second_new_lex_G;
   

######################################

#def fast_compute_lex(Fq,Fq_star,n,a,b,w):
   
   
#    G = ones_matrix(Fq,2,n);
    
    #Find minimum codewords in the code
#    min_codewords = [];
#    min_w = min(len(a.support()), len(b.support()));
#    almost_min_codewords = [a,b]; 
#    for u in Fq_star:
#        c = a+u*b;
#        wc = n - c.list().count(0);
#        if wc<=min_w:
#            if wc<min_w:
#                min_w = wc;
#                min_codewords = [c];
#                almost_min_codewords = min_codewords;
#            else:
#                min_codewords.append(c);
               
#    num_min_codewords = len(min_codewords);        
#    if num_min_codewords < 2:
#        for c in almost_min_codewords:
#            min_codewords.append(c);
#    num_min_codewords = len(min_codewords);      
#    print(min_codewords)
       
#    for i in range(num_min_codewords-1):
#        a = min_codewords[i];
#        for j in range(i+1,num_min_codewords):
#            b = min_codewords[j];
#            for u in Fq_star:
               
#                new_G = minimal_matrix(Fq,n,a,u*b);
#                G = compare_matrix_lex(G,new_G,n);
               
#                new_G = minimal_matrix(Fq,n,u*a,b);
#                G = compare_matrix_lex(G,new_G,n);
               

#    return G;

#########################################################

###slow version
#def compute_lex(Fq,Fq_star,n,a,b,w):
   
   
#    G = ones_matrix(Fq,2,n);
#    S = matrix(Fq,2,2);
    
#    for u00 in Fq:
#        for u01 in Fq:
#            new_a = u00*a+u01*b;
#            for u10 in Fq:
#                for u11 in Fq:

#                    new_b = u10*a+u11*b;
                    
#                    S[0,0] = u00;
#                    S[0,1] = u01;
#                    S[1,0] = u10;
#                    S[1,1] = u11;
#                    if rank(S) == 2:
#                        new_G = minimal_matrix(Fq,n,new_a,new_b);
                            
                        #   print("---------------->");
                        #   print(new_G)
                        #   print("---------------->")
                        #   print(G)
#                        G = compare_matrix_lex(G,new_G,n);
                        #   print(G)
                        #  print(G);
                        #   print("-------------")
               

#    return G;

######################################à

###slow version
def compute_lex(Fq,Fq_star,n,a,b):
   
   
    G = ones_matrix(Fq,2,n);
    S = matrix(Fq,2,2);
    
    u00 = Fq(1);
    for u01 in Fq:
        new_a = a+u01*b;
        for u10 in Fq:
            for u11 in Fq:

                new_b = u10*a+u11*b;
                
                S[0,0] = u00;
                S[0,1] = u01;
                S[1,0] = u10;
                S[1,1] = u11;
                if rank(S) == 2:
                    new_G = minimal_matrix(Fq,n,new_a,new_b);
                        
                    G = compare_matrix_lex(G,new_G,n);
    u00 = Fq(0);
    for u01 in Fq:
        new_a = u01*b;
        for u10 in Fq:
            for u11 in Fq:

                new_b = u10*a+u11*b;
                
                S[0,0] = u00;
                S[0,1] = u01;
                S[1,0] = u10;
                S[1,1] = u11;
                if rank(S) == 2:
                    new_G = minimal_matrix(Fq,n,new_a,new_b);
                        
                    G = compare_matrix_lex(G,new_G,n);
     

    return G;

#q = 11;
#k = 2;
#n = 10;

#num_test = 10^3;
#Fq = GF(q);
#Fq_star = Set(Set(Fq)[1:q]);
#num_err = 0;
#while num_err == 0:
#for i in range(num_test):
#    r = 0;
#    while r<2:
#        G = random_matrix(Fq,2,n);
#        r = rank(G);
    
#    if compute_lex(Fq,Fq_star,n,vector(G[0,:]), vector(G[1,:]),10)!=fast_compute_lex(Fq,Fq_star,n,vector(G[0,:]), vector(G[1,:]),10):
#        num_err += 1;
        
#    print("Num errors = "+str(num_err)+", num tests = "+str(i+1));
