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

#pragma once
#include "parameters.h"
#include "assert.h"
#include "utils.h"

#define NUM_BITS_Q (BITS_TO_REPRESENT(Q))


/* GCC actually inlines and vectorizes Barrett's reduction already.
 * Backup implementation for less aggressive compilers follows

#define BARRETT_MU  (((uint32_t)1<<(2*NUM_BITS_Q))/Q)
#define BARRETT_MASK ( ((FQ_DOUBLEPREC)1 << (NUM_BITS_Q+3))-1 )

static inline
FQ_ELEM fq_red(FQ_DOUBLEPREC a) {
    FQ_DOUBLEPREC q_1, q_2, q_3;
    q_1 = a >> (NUM_BITS_Q);
    q_2 = q_1 * BARRETT_MU;
    q_3 = q_2 >> (NUM_BITS_Q);
    FQ_DOUBLEPREC r_1;
    r_1 = (a & BARRETT_MASK) - ( (q_3*Q) & BARRETT_MASK);
    r_1 = r_1 & BARRETT_MASK;
    FQ_ELEM r_2;
    FQ_DOUBLEPREC need_to_red;  
    need_to_red = r_1 >= Q;
    r_1 = r_1-Q*need_to_red; // not needed for 127
    need_to_red = r_1 >= Q;
    r_2 = r_1-Q*need_to_red;
    return r_2;
} */

static inline
FQ_ELEM fq_red(FQ_DOUBLEPREC x){
    return ((FQ_DOUBLEPREC) Q+x) % (FQ_DOUBLEPREC) Q;
}

/* Fermat's method for inversion employing r-t-l square and multiply, 
 * unrolled for actual parameters */
static
FQ_ELEM fq_inv(FQ_ELEM x){
    assert(x != 0);
    FQ_DOUBLEPREC xlift;
    xlift = x;
    FQ_DOUBLEPREC accum = 1;
    /* No need for square and mult always, Q-2 is public*/
    uint32_t exp = Q-2;
    while(exp){
        if(exp & 1) {
            accum = fq_red(accum*xlift);
        }
        xlift = fq_red(xlift*xlift);
        exp >>= 1;
    }
    return fq_red(accum);
}

static inline
FQ_ELEM fq_star_rnd(){
   const uint32_t mask = ( (uint32_t) 1 << BITS_TO_REPRESENT(Q-2) ) - 1;
   uint32_t rnd_value;
   do {
      rnd_value = csprng();
      rnd_value = mask & rnd_value;
   } while (rnd_value > Q-2);

    return rnd_value+1;
}

void compute_systematic_form(FQ_ELEM to_gausselim[N][K]);
