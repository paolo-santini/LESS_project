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
#include "monomial_mat.h"
#include "codes.h"
#include <assert.h>

uint32_t csprng();

static inline
uint32_t rand_range_n(){
   const uint32_t mask = ( (uint32_t) 1 << BITS_TO_REPRESENT(N-1)) - 1;
   uint32_t rnd_value;
   do {
      rnd_value = csprng();
      rnd_value = mask & rnd_value;
   } while (rnd_value >= N);
   assert(rnd_value < N);
   return rnd_value;
}

void hash( uint8_t digest[DENSE_HASH_LENGTH],
             const char * const m,
             const uint64_t mlen,
             generator_mat_t G_tilde[T]);

void parse_digest( uint8_t parsed_digest[PARSED_DIGEST_LEN],
                   const uint8_t digest[DENSE_HASH_LENGTH]);
