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

#include "parameters.h"
#include "fq_arith.h"
#include "monomial_mat.h"
#include "utils.h"
#include <stdio.h>
#include <assert.h>



void monomial_mat_seed_rnd(monomial_seed_t *res){
   randombytes(res->value, SEED_LENGTH_BYTES);
}


void monomial_mat_seed_expand(monomial_t *res,
                              const monomial_seed_t * const seed){
    SHAKE_STATE_STRUCT shake_monomial_state = {0};
    initialize_prng(&shake_monomial_state,seed->value);
    for(int i = 0; i < N; i++){
        res->coefficients[i] = fq_star_rnd_state(&shake_monomial_state);
        res->permutation[i] = i;
    }
    /* FY shuffle on the permutation */
    POSITION_T tmp;
    for(int i = 0; i < N; i++){
        POSITION_T rnd = rand_range_n_state(&shake_monomial_state);

        tmp = res->permutation[i];
        res->permutation[i] = res->permutation[rnd];
        res->permutation[rnd] = tmp;
    }    
}


/* samples a random perm matrix */
void monomial_mat_rnd(monomial_t *res){
    for(int i = 0; i < N; i++){
        res->coefficients[i] = fq_star_rnd();
        assert(res->coefficients[i] !=0);
        res->permutation[i] = i;
    }
    /* FY shuffle on the permutation */
    POSITION_T tmp;
    for(int i = 0; i < N; i++){
        POSITION_T rnd = rand_range_n();

        tmp = res->permutation[i];
        res->permutation[i] = res->permutation[rnd];
        res->permutation[rnd] = tmp;
    }
}

void monomial_mat_mul(monomial_t *res,
                  const monomial_t * const A,
                  const monomial_t * const B){
    for(int i = 0; i < N; i++){
         res->permutation[i] = B->permutation[A->permutation[i]];
         res->coefficients[i] = fq_red(
                                (FQ_DOUBLEPREC) A->coefficients[i] *
                                (FQ_DOUBLEPREC) B->coefficients[A->permutation[i]] );
    }
}

void monomial_mat_inv(monomial_t *res,
                      const monomial_t * const to_invert){
    for(int i = 0; i < N; i++){
        res->permutation[to_invert->permutation[i]] = i;
        res->coefficients[to_invert->permutation[i]] = fq_inv(to_invert->coefficients[i]);        
    }
}

/* yields the identity matrix */
void monomial_mat_id(monomial_t *res){
    for(int i = 0; i < N; i++){
        res->permutation[i] = i;
        res->coefficients[i] = 1;
    }
}

/* pretty_print for monomial matrices */
void monomial_mat_pretty_print(const monomial_t * const to_print){
    fprintf(stderr,"perm = [");
    for(int i = 0; i < N-1; i++){
        fprintf(stderr,"%03u, ",to_print->permutation[i]);
    }
    fprintf(stderr,"%03u ]\n",to_print->permutation[N-1]);
    fprintf(stderr,"coeffs = [");
    for(int i = 0; i < N-1; i++){
        fprintf(stderr,"%03u, ",to_print->coefficients[i]);
    }
    fprintf(stderr,"%03u ]\n",to_print->coefficients[N-1]);    
}

void monomial_mat_pretty_print_name(char* name, const monomial_t *to_print){
    fprintf(stderr,"%s = [",name);
    for(int i = 0; i < N-1; i++){
        fprintf(stderr,"%03u, ",to_print->permutation[i]);
    }
    fprintf(stderr,"%03u ]\n",to_print->permutation[N-1]);
    fprintf(stderr,"coeffs = [");
    for(int i = 0; i < N-1; i++){
        fprintf(stderr,"%03u, ",to_print->coefficients[i]);
    }
    fprintf(stderr,"%03u ]\n",to_print->coefficients[N-1]);    
}

void monomial_mat_print_exp_name(char* name,const monomial_t *to_print){
    FQ_ELEM mu[N][N]= {{0}};
  
    for(int i = 0; i < N; i++){
        mu[to_print->permutation[i]][i] = to_print->coefficients[i];
    }
    
    fprintf(stderr,"%s = Mon([",name);
    for(int i = 0; i < N-1 ; i++ ){
       fprintf(stderr,"[");
       for(int j = 0; j < N-1; j++){
           fprintf(stderr,"%u, ",mu[i][j]);
       }
       fprintf(stderr,"%u ],\n",mu[i][N-1]);
    }
    fprintf(stderr,"[");
    for(int j = 0; j < N-1; j++){
        fprintf(stderr,"%u, ",mu[N-1][j]);
    }
    fprintf(stderr,"%u ] ])\n",mu[N-1][N-1]);
}

