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

#include <stdlib.h>
#include <stdint.h>
#include "rng.h"
#include "sha3.h"
#include "utils.h"
#include <string.h>

uint32_t csprng(){
    uint32_t x;
    randombytes((unsigned char*) &x, sizeof(uint32_t));
    return x;
}

void hash( uint8_t digest[DENSE_HASH_LENGTH],
             const char * const m,
             const uint64_t mlen,
             generator_mat_t G_tilde[T]){
    unsigned char message_buffer[ mlen + T*sizeof(generator_mat_t) ];
    memcpy(message_buffer, G_tilde, T*sizeof(generator_mat_t) );
    memcpy(message_buffer + T*sizeof(generator_mat_t), m, mlen);
    sha3_256(message_buffer, mlen+T*sizeof(generator_mat_t), digest);
}

/* parses a digest expanding it according to the LESS variant requirements
 * it either expands bitwise the digest or generates a fixed-weight string
 * ell bit sized elements */
void parse_digest( uint8_t parsed_digest[PARSED_DIGEST_LEN],
                   const uint8_t digest[DENSE_HASH_LENGTH]){
#if defined(LESSF)
    AES_XOF_struct parser_ctx;
    unsigned char empty_diversifier[8] = {0};
    seedexpander_init(&parser_ctx,
                      (unsigned char*) digest,
                      empty_diversifier,
                      (unsigned int) (1<<20));
    int placed_elements = 0;
    uint16_t rnd_buf;
    while (placed_elements < W){
        uint8_t value;
        POSITION_T pos;
        do {
            seedexpander(&parser_ctx, (unsigned char *) &rnd_buf, sizeof(uint16_t));
            value = rnd_buf & ( ( (uint16_t)1 << BITS_TO_REPRESENT(NUM_KEYPAIRS-1) )-1);
            pos   = rnd_buf >> BITS_TO_REPRESENT(NUM_KEYPAIRS-1) ;
            pos = pos & (((POSITION_T)1 << BITS_TO_REPRESENT(PARSED_DIGEST_LEN-1))-1);
        } while ( (value >= NUM_KEYPAIRS) ||
                  (pos  >= PARSED_DIGEST_LEN) ||
                  (parsed_digest[pos] != 0) );
        parsed_digest[pos] = value;
        placed_elements++;
    }
#else
#error LESS with variable length digest is not supported
#endif
}
