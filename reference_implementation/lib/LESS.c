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

#include "LESS.h"
#include "utils.h"
#include <assert.h>
#include <string.h>
#include <stdio.h>

void LESS_keygen(prikey_t * SK,
                 pubkey_t * PK){
    /* the first monomial matrix is an identity */
    monomial_mat_id(&SK->private_Q[0]);
    monomial_mat_id(&SK->private_Q_inv[0]);
    do {
        generator_rnd(&PK->Full_G);
        memcpy(&PK->SF_G[0], &PK->Full_G, sizeof(generator_mat_t) );
    } while ( generator_gausselim(&PK->SF_G[0]) == 0);

    for(int i = 1; i < NUM_KEYPAIRS; i++){
        do {
           monomial_mat_rnd(&SK->private_Q[i]);
           generator_monomial_mul(&PK->SF_G[i],
                                  &PK->Full_G,
                                  &SK->private_Q[i]);
        } while ( generator_gausselim(&PK->SF_G[i]) == 0);
        monomial_mat_inv(&SK->private_Q_inv[i], &SK->private_Q[i]);
    }
}

void LESS_sign(const prikey_t * SK,
              const pubkey_t * PK,
              const char * const m,
              const uint64_t mlen,
              sig_t * sig) {
    monomial_t Q_tilde[T];
    generator_mat_t G_tilde[T];
    for(int i = 0; i < T; i++){
        do {
           monomial_mat_rnd(&Q_tilde[i]);
           generator_monomial_mul(&G_tilde[i],
                                  &PK->Full_G,
                                  &Q_tilde[i]);
        } while (generator_gausselim(&G_tilde[i]) == 0);
    }

   hash(sig->digest, m, mlen, G_tilde);

   uint8_t parsed_digest[PARSED_DIGEST_LEN] = {0};
   parse_digest(parsed_digest,sig->digest);
   for(int i = 0; i < T; i++){
      monomial_mat_mul(&sig->mu[i],
                       &SK->private_Q_inv[parsed_digest[i]],
                       &Q_tilde[i]);
   }
}

int LESS_verify(const pubkey_t * const PK,
                const char * const m,
                const uint64_t mlen,
                const sig_t * const sig){
    uint8_t parsed_digest[PARSED_DIGEST_LEN] = {0};
    generator_mat_t G_hat[T];

    parse_digest(parsed_digest,sig->digest);
    for(int i = 0; i < T; i++){
        generator_monomial_mul(&G_hat[i],
                               &PK->SF_G[parsed_digest[i]],
                               &sig->mu[i]);
        int is_gausselim_ok;
        is_gausselim_ok = generator_gausselim(&G_hat[i]);
        assert(is_gausselim_ok);
    }

    uint8_t check_digest[DENSE_HASH_LENGTH] = {0};
    hash(check_digest,m,mlen,G_hat);
    int does_hash_match = 0;
    does_hash_match = (memcmp(check_digest,sig->digest,DENSE_HASH_LENGTH) == 0);
    return does_hash_match;
}
