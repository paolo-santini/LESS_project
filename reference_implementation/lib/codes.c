/**
 *
 * Reference ISO-C11 Implementation of LESS.
 *
 * @version 1.0 (February 2022)
 *
 * @author Alessandro Barenghi <alessandro.barenghi@polimi.it>
 *
 * This code is hereby placed in the public domain.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 **/

#include "monomial_mat.h"
#include "fq_arith.h"
#include "codes.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>
#define ASSERT_LIMITS(x) assert((x>=0)&&(x<Q))
// 
/* right-multiplies a generator by a monomial */
void generator_monomial_mul(generator_mat_t *res,
                            const generator_mat_t *G,
                            const monomial_t * const monom){
   for(int src_col_idx = 0; src_col_idx < N; src_col_idx++){
       for(int row_idx = 0; row_idx < K; row_idx++){        
            res->values[row_idx][monom->permutation[src_col_idx]] =
                  fq_red( (FQ_DOUBLEPREC) G->values[row_idx][src_col_idx] * 
                          (FQ_DOUBLEPREC) monom->coefficients[src_col_idx] );
       }
   }
}


static inline
void swap_cols_and_rescale_dest(generator_mat_t* M,
                                POSITION_T dst_pos_in_M, 
                                FQ_ELEM src[K],
                                FQ_ELEM scaling_factor){
   FQ_DOUBLEPREC tmp;
   for(int i=0; i<K; i++){
       tmp  = src[i];
       src[i] = M->values[i][dst_pos_in_M];  
       M->values[i][dst_pos_in_M] = fq_red( (FQ_DOUBLEPREC)tmp * 
                                            (FQ_DOUBLEPREC)scaling_factor );
   }
}


/* TODO FIX: right-multiplies a generator by a monomial in-place*/
void generator_monomial_mul_in_place(generator_mat_t *G,
                                     const monomial_t * const monom){
   monomial_t local_monom;
   memcpy(&local_monom, monom, sizeof(monomial_t));
   for(POSITION_T src_col_idx = 0; src_col_idx < N; src_col_idx++){
       FQ_ELEM temp_col[K] = {0};
       POSITION_T tgt_col_idx;
       FQ_ELEM tgt_coeff;
       
       for(int row_idx=0; row_idx < K; row_idx++){
           temp_col[row_idx] = (FQ_DOUBLEPREC) G->values[row_idx][src_col_idx];
       }
       
       tgt_col_idx =  local_monom.permutation[src_col_idx]; 
       tgt_coeff   = local_monom.coefficients[src_col_idx];
       if (tgt_col_idx == src_col_idx) continue;
       while(tgt_col_idx != src_col_idx){
           swap_cols_and_rescale_dest(G, 
                                      tgt_col_idx, 
                                      temp_col, 
                                      tgt_coeff);
           FQ_ELEM tmp_coeff;
           tmp_coeff = tgt_coeff;
           tgt_coeff = local_monom.coefficients[tgt_col_idx];
           local_monom.coefficients[tgt_col_idx] = tmp_coeff;
           
           POSITION_T tmp_idx;
           tmp_idx = local_monom.permutation[tgt_col_idx] ; 
           local_monom.permutation[tgt_col_idx] = tgt_col_idx;
           tgt_col_idx = tmp_idx;
       }
       swap_cols_and_rescale_dest(G, 
                                  src_col_idx, 
                                  temp_col, 
                                  tgt_coeff);       
       local_monom.coefficients[tgt_col_idx] = tgt_coeff; 
       local_monom.permutation[tgt_col_idx] = tgt_col_idx;    
   }
}


/* right-multiplies a generator by a monomial, input generator in compact form */
void rref_generator_monomial_mul(generator_mat_t *res,
                            const generator_mat_t *G,
                            const monomial_t * const monom){
   for(int src_col_idx = 0; src_col_idx < N; src_col_idx++){
       for(int row_idx = 0; row_idx < K; row_idx++){
            res->values[row_idx][monom->permutation[src_col_idx]] =
                  fq_red( (FQ_DOUBLEPREC) G->values[row_idx][src_col_idx] * 
                          (FQ_DOUBLEPREC) monom->coefficients[src_col_idx] );
       }
   }
}


static inline
void swap_rows(FQ_ELEM r[N],FQ_ELEM s[N]){
   FQ_ELEM tmp;
   for(int i=0; i<N; i++){
       tmp = r[i];
       r[i] = s[i];
       s[i] = tmp;
   }
}

int generator_gausselim(generator_mat_t * G,
                        uint8_t is_pivot_column[N]){
    for(int row_to_reduce = 0; row_to_reduce < K; row_to_reduce++){
        int pivot_row = row_to_reduce;
        /*start by searching the pivot in the col = row*/
        int pivot_column = row_to_reduce;
        while( (pivot_column < N) && 
               (G->values[pivot_row][pivot_column] == 0) ){

           while ( (pivot_row < K) && 
                   (G->values[pivot_row][pivot_column] == 0) ){
              pivot_row++;
           }
           if(pivot_row >= K){ /*entire column tail swept*/ 
              pivot_column++; /* move to next col */
              pivot_row = row_to_reduce; /*starting from row to red */
           }
        } 
        if ( pivot_column >=N ){
                return 0; /* no pivot candidates left, report failure */           
        }
//         if (pivot_column != row_to_reduce) fprintf(stderr,"aaa")
        is_pivot_column[pivot_column] = 1; /* pivot found, mark the column*/
        assert(G->values[pivot_row][pivot_column] != 0);
        
        /* if we found the pivot on a row which has an index > pivot_column 
         * we need to swap the rows */
        if (row_to_reduce != pivot_row) {
           swap_rows(G->values[row_to_reduce],G->values[pivot_row]);
        }
        pivot_row = row_to_reduce; /* row with pivot now in place */
        
        assert(G->values[pivot_row][pivot_column] != 0);
        
        /* Compute rescaling factor */
        FQ_DOUBLEPREC scaling_factor = fq_inv(G->values[pivot_row][pivot_column]);
        assert(scaling_factor != 0);
        ASSERT_LIMITS(scaling_factor);
        /* rescale pivot row to have pivot = 1. Values at the left of the pivot
         * are already set to zero by previous iterations */
        for(int i = pivot_column; i < N; i++){
            ASSERT_LIMITS(G->values[pivot_row][i]);
            G->values[pivot_row][i] = fq_red( (FQ_DOUBLEPREC) scaling_factor *
                                (FQ_DOUBLEPREC) (G->values[pivot_row][i]) );
        }     
        assert(G->values[pivot_row][pivot_column] == 1); /* pivot row now rescaled*/

        /* Subtract the now placed and reduced pivot rows, from the others, 
         * after rescaling it */
        for(int row_idx = 0; row_idx < K; row_idx++){
            if (row_idx != pivot_row){
                FQ_DOUBLEPREC multiplier = G->values[row_idx][pivot_column];
                /* all elements before the pivot in the pivot row are null, no need to
                 * subtract them from other rows. */
                for(int col_idx = 0; col_idx < N; col_idx++){
                    FQ_DOUBLEPREC tmp;
                    tmp = fq_red( (FQ_DOUBLEPREC) multiplier * 
                                  (FQ_DOUBLEPREC) G->values[pivot_row][col_idx] );
                    ASSERT_LIMITS(tmp);
                    tmp = (FQ_DOUBLEPREC) Q + (FQ_DOUBLEPREC) G->values[row_idx][col_idx] - tmp;
                    tmp = fq_red(tmp);
                    ASSERT_LIMITS(tmp);

                    G->values[row_idx][col_idx] = tmp;
                }
                assert(G->values[row_idx][pivot_column] == 0); /* test removal from pivot column */
            }
        }
    }
    return 1;
}

void generator_rref_compact(rref_generator_mat_t* compact,
                          const generator_mat_t* const full,
                          const uint8_t is_pivot_column[N] ){
    int dst_col_idx = 0;
    for (int src_col_idx = 0; src_col_idx < N; src_col_idx++){
        if(!is_pivot_column[src_col_idx]){
            for (int row_idx = 0; row_idx < K; row_idx++){
                    compact->values[row_idx][dst_col_idx] = full->values[row_idx][src_col_idx];
            }
            compact->column_pos[dst_col_idx] = src_col_idx;
            dst_col_idx++;
        }
    }
}

void generator_rref_expand(generator_mat_t * full,
                         const rref_generator_mat_t * const compact){
   int placed_dense_cols = 0;
   for (int col_idx = 0; col_idx < N; col_idx++){
       if ( (placed_dense_cols< N-K) && 
            (col_idx == compact->column_pos[placed_dense_cols])) {
           /* non-pivot column, restore one full column */
           for (int row_idx = 0; row_idx < K; row_idx++){
                full->values[row_idx][col_idx] = compact->values[row_idx][placed_dense_cols];
           }
           placed_dense_cols++;
       } else {
           /* regenerate the appropriate pivot column */
           for (int row_idx = 0; row_idx < K; row_idx++){
                    full->values[row_idx][col_idx] = (row_idx == col_idx-placed_dense_cols);
           }
       }
   }  
}

/* samples a random generator matrix */
void generator_rnd(generator_mat_t *res){
    for(int i = 0; i < K; i++){
        for(int j = 0; j < N ; j++){
            res->values[i][j]=rand_range_q();
        }
    }
}

/* pretty_print for full generator matrices */
void generator_pretty_print_name(char* name, const generator_mat_t * const G){
    fprintf(stderr,"%s = M([",name);
    for(int i = 0; i < K-1 ; i++ ){
       fprintf(stderr,"[");
       for(int j = 0; j < N-1; j++){
           fprintf(stderr,"%u, ",G->values[i][j]);
       }
       fprintf(stderr,"%u ],\n",G->values[i][N-1]);
    }
    fprintf(stderr,"[");
    for(int j = 0; j < N-1; j++){
        fprintf(stderr,"%u, ",G->values[K-1][j]);
    }
    fprintf(stderr,"%u ] ])\n",G->values[K-1][N-1]);
}

/* pretty_print for generator matrices in row-reduced echelon form*/
void generator_rref_pretty_print_name(char* name, const rref_generator_mat_t * const G){
    fprintf(stderr,"%s =\n[",name);
    for(int i = 0; i < K-1 ; i++ ){
       fprintf(stderr,"[");
       for(int j = 0; j < (N-K)-1; j++){
           fprintf(stderr,"%u, ",G->values[i][j]);
       }
       fprintf(stderr,"%u ],\n",G->values[i][(N-K)-1]);
    }
    fprintf(stderr,"[");
    for(int j = 0; j < (N-K)-1; j++){
        fprintf(stderr,"%u, ",G->values[K-1][j]);
    }
    fprintf(stderr,"%u ] ]\n",G->values[K-1][(N-K)-1]);
    fprintf(stderr,"column_pos = \n [ ");
    for(int x=0;x < K ;x++){fprintf(stderr," %d ",G->column_pos[x]); }
    fprintf(stderr,"]\n");
       
}
