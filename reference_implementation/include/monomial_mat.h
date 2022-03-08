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

 /* Structure representing a monomial matrix */
typedef struct {
    /* coefficients listed in order of appearance columnwise */
    FQ_ELEM coefficients[N];
    /* stores the coefficient of the destination position, i.e., assuming
     * a GQ product, the i-th element is the destination column index of the i-th
     * column of G */
    POSITION_T permutation[N];
} monomial_t;

/* multiplies two monomial matrices*/
void monomial_mat_mul(monomial_t *res,
                      const monomial_t * const A,
                      const monomial_t * const B);

/* computes the inverse of the monomial matrix */
void monomial_mat_inv(monomial_t *res,
                      monomial_t *to_invert);

/* samples a random monomial matrix */
void monomial_mat_rnd(monomial_t *res);

/* yields the identity matrix */
void monomial_mat_id(monomial_t *res);

/* pretty_print for monomial matrices */
void monomial_mat_pretty_print(const monomial_t *res);
