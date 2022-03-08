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

typedef struct {
    FQ_ELEM values[K][N];
} generator_mat_t;

typedef struct {
    FQ_ELEM values[K][N-K];
} sf_generator_mat_t;

/* multiplies a monomial matrix by a generator matrix */
void generator_monomial_mul(generator_mat_t *res,
                            const generator_mat_t *G,
                            const monomial_t * const monom);

/* computes the normal form of the generator matrix with the first k columns as minor
 * returns 1 on success, 0 on failure, computation is done in-place */
int generator_gausselim(generator_mat_t *res);

/* samples a random monomial matrix */
void generator_rnd(generator_mat_t *res);

void generator_pretty_print(const generator_mat_t * const G);
