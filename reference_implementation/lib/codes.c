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

/* right-multiplies a generator by a monomial */
void generator_monomial_mul(generator_mat_t *res,
                            const generator_mat_t *G,
                            const monomial_t * const monom){
   for(int src_col_idx = 0; src_col_idx < N; src_col_idx++){
       for(int row_idx = 0; row_idx < K; row_idx++){
            res->values[row_idx][monom->permutation[src_col_idx]] =
                  fq_red( (FQ_DOUBLEPREC) G->values[row_idx][src_col_idx] * (FQ_DOUBLEPREC) monom->coefficients[src_col_idx]);
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

/* computes the RREF of the generator matrix in place with the first k columns
 * as minor returns 1 on success, 0 on failure. Conventionally RREF(G) = [I | V]
 * In case a singular minor is found, the RREF procedure stops. *
 * No information leakage occurs in upper layers as the matrix is discarded */

int generator_gausselim(generator_mat_t * G){
    for(int row_to_reduce = 0; row_to_reduce < K; row_to_reduce++){
        int pivot_row = row_to_reduce;
        int pivot_column = row_to_reduce;
        while ( (pivot_row < K) &&
                (G->values[pivot_row][pivot_column] == 0) ){
            pivot_row++;
        }
        if (pivot_row == K) {
            return 0; /* no pivots left */
        }

        assert(G->values[pivot_row][pivot_column] != 0);
        FQ_DOUBLEPREC scaling_factor = fq_inv(G->values[pivot_row][pivot_column]);
        assert(scaling_factor!=0);
        /* rescale pivot row first */
        for(int i = pivot_column; i < N; i++){
            G->values[pivot_row][i] = fq_red( scaling_factor *
                                (FQ_DOUBLEPREC) (G->values[pivot_row][i]) );
        }
        assert(G->values[pivot_row][pivot_column] == 1);
        if (row_to_reduce != pivot_row) {
           swap_rows(G->values[row_to_reduce],G->values[pivot_row]);
        }
        for(int row_idx = 0; row_idx < K; row_idx++){
            if (row_idx != row_to_reduce){
                FQ_DOUBLEPREC multiplier = G->values[row_idx][pivot_column];
                for(int col_idx = pivot_column; col_idx < N; col_idx++){
                    G->values[row_idx][col_idx] =
                       fq_red( ( (FQ_DOUBLEPREC)Q + (FQ_DOUBLEPREC) G->values[row_idx][col_idx]) -
                               fq_red( (FQ_DOUBLEPREC) multiplier *
                                       (FQ_DOUBLEPREC) G->values[row_to_reduce][col_idx] )
                             );
                }
                assert(G->values[row_idx][pivot_column] == 0);
            }
        }
    }
    return 1;
}

/* samples a random monomial matrix */
void generator_rnd(generator_mat_t *res){
    for(int i = 0; i < K; i++){
        for(int j = 0; j < N ; j++){
            res->values[i][j]=rand_range_n();
        }
    }
}

/* pretty_print for monomial matrices */
void generator_pretty_print(const generator_mat_t * const G){
    fprintf(stderr,"G = [");
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
    fprintf(stderr,"%u ] ]\n",G->values[K-1][N-1]);
}
