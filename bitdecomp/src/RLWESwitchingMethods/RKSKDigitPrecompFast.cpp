
#include <utility>
#include <cassert>

#include "RLWESwitchingMethods/RKSKDigitPrecompFast.h"


#define BLOCK_SIZE 8
#define MAX_DIGITS 8
#define ALIGNMENT 64

alignas(ALIGNMENT) UINTNAT bbuffer[BLOCK_SIZE];

RKSKDigitPrecompFast::RKSKDigitPrecompFast(std::shared_ptr<RingGSWCryptoParams> params, uint64_t L_bits, uint32_t max_ct) :
        m_params(std::move(params)), m_max_ct(max_ct), m_L(1ull << L_bits), m_L_bits(L_bits)
    {
}

RKSKDigitPrecompFast::RKSKDigitPrecompFast(std::shared_ptr<RingGSWCryptoParams> params, uint64_t L_bits,
                                           uint64_t L_LWE_bits, uint64_t n_LWE, uint32_t max_ct) :
        m_params(std::move(params)), m_max_ct(max_ct), m_L(1ull << L_bits), m_L_bits(L_bits), m_LWE_n(n_LWE),
        m_LWE_L(1ull << L_LWE_bits), m_LWE_L_bits(L_LWE_bits), m_use_intermediate(true)

                                           {
}

void RKSKDigitPrecompFast::Keygen(NativeVector &sk_in, NativePoly &sk_out) {

    auto Qnat = m_params->GetQ();
    auto Q = Qnat.ConvertToInt<UINTNAT>();
    auto N = m_params->GetN();
    auto dgg = m_params->GetDgg();
    auto PolyParams = m_params->GetPolyParams();

    sk_out.SetFormat(EVALUATION);
    NativeVector intermediate_key;

    uint64_t Qbits = 128 - clz_u128(Q) - 1;
    m_digits = Qbits / m_L_bits + (Qbits % m_L_bits != 0);

    // Can we even use the fast one ?
    auto Nbits = 128 - clz_u128(N) - 1;
    auto digit_bits = 128 - clz_u128(std::max(m_digits, m_LWE_digits)) - 1;

    auto acc_bits = Qbits + Nbits + digit_bits + 1; // safety factor
    assert(acc_bits < 64);
    // needed, otherwise clang thinks we"re not using acc_bits
    ((void)(acc_bits));
    // check succeeded

    DiscreteUniformGeneratorImpl<NativeVector> dug;
    dug.SetModulus(Q);

    auto n_pad = m_params->GetN() / m_max_ct;
    for(uint32_t i = 0; i < m_max_ct; i++) {
        auto Xfi = NativePoly(PolyParams, COEFFICIENT, true);
        Xfi[i * n_pad] = 1;
        Xfi.SetFormat(EVALUATION);
        m_offset_polys.push_back(Xfi);
    }

    if (m_use_intermediate) {

        // generate intermediate keyswitching key
        intermediate_key = dug.GenerateVector(m_LWE_n);


        m_LWE_digits = Qbits / m_LWE_L_bits + (Qbits % m_LWE_L_bits != 0);

        // allocate memory
        m_LWE_ksk_a = static_cast<UINTNAT*>(aligned_alloc(ALIGNMENT,  N * (m_LWE_n) * m_LWE_digits * (m_LWE_L - 1) * sizeof(UINTNAT)));
        m_LWE_ksk_b = static_cast<UINTNAT*>(aligned_alloc(ALIGNMENT,  N * m_LWE_digits * (m_LWE_L - 1) * sizeof(UINTNAT)));

        // Generate first ksk
        for(uint32_t i = 0; i < N; i++) {
            auto ski = sk_in[i];
            for(uint32_t j = 0; j < m_LWE_digits; j++) {
                for(uint32_t k = 1; k < m_LWE_L; k++) {
                    auto skijk = ski.ModMul(k, Qnat);
                    auto vec_a = dug.GenerateVector(m_LWE_n);
                    auto vec_b = skijk.ModAdd(0, Qnat); //dgg.GenerateInteger(Qnat),Qnat);

                    auto offset_ijk = i * m_LWE_n * m_LWE_digits * (m_LWE_L - 1) + j * m_LWE_n * (m_LWE_L - 1) + (k - 1) * m_LWE_n;

                    // encrypt and store the A component at the same time
                    for(uint32_t kk = 0; kk < m_LWE_n; kk++) {
                        m_LWE_ksk_a[offset_ijk + kk] = vec_a[kk].ConvertToInt<UINTNAT>();
                        auto askk = vec_a[kk].ModMul(intermediate_key[kk], Q);
                        vec_b.ModAddEq(askk, Q);
                    }

                    // store the b component
                    m_LWE_ksk_b[(k - 1) + j * (m_LWE_L - 1) + i * m_LWE_digits * (m_LWE_L - 1)] = vec_b.ConvertToInt<UINTNAT>();
                }
                ski.ModMulEq(m_LWE_L, Qnat);
            }
        }


    } else {
        intermediate_key = sk_in;
    }

    auto int_N = intermediate_key.GetLength();
    auto chunksize = N / m_max_ct;

    auto Xfi = NativePoly(PolyParams, COEFFICIENT, true);
    for(uint32_t j = 0; j < chunksize; j++) {
        Xfi[j] = 1;
    }
    Xfi.SetFormat(EVALUATION);

    m_LWE_out_buffer = static_cast<uint64_t*>(aligned_alloc(ALIGNMENT,  N * sizeof(uint64_t)));
    m_ksk = static_cast<UINTNAT*>(aligned_alloc(ALIGNMENT, 2 * N * N * m_digits * (m_L - 1) * sizeof(UINTNAT)));
    m_out_buffer = static_cast<uint64_t*>(aligned_alloc(ALIGNMENT, 2 * N * m_max_ct * sizeof (uint64_t)));

    for(uint32_t i = 0; i < int_N; i++) {
        NativeInteger ski = intermediate_key[i];
        NativeInteger factor = 1;
        auto offset_i = i * m_digits * (m_L - 1);
        /* go other every power of the basis */
        for(uint32_t j = 0; j < m_digits; j++) {
            auto skij = ski.ModMul(factor, Q);
            auto offset_j = j * (m_L - 1);
            for(uint32_t k = 1; k < m_L; k++) {
                auto skijk = skij.ModMul(k, Q);

                //NativePoly a_poly(PolyParams, EVALUATION, true);
                NativePoly a_poly(dug, PolyParams, EVALUATION);

                NativePoly b_poly(PolyParams, COEFFICIENT, true);
                //NativePoly b_poly(dgg, PolyParams, COEFFICIENT);
                b_poly.SetFormat(EVALUATION);
                b_poly += skijk * Xfi;

                b_poly += a_poly * sk_out;

                auto col_idx = offset_i + offset_j + k - 1;

                // write
                for(uint32_t kk = 0; kk < N; kk++) {
                    m_ksk[col_idx * 2 * N + kk] = a_poly[kk].ConvertToInt<UINTNAT>();
                    m_ksk[col_idx * 2 * N + kk + N] = b_poly[kk].ConvertToInt<UINTNAT>();
                }
            }
            factor.MulEq(m_L);
        }
    }

}

template<uint32_t n_digits> void MergeDigits(uint64_t* __restrict outputs, UINTNAT** i_addr, uint64_t N2) {
    static_assert(n_digits > 1, "digits must be greater 1");

    // small tiling buffer

    for(uint32_t i = 0; i < N2 / BLOCK_SIZE; i++) {

        for(uint32_t k = 0; k < BLOCK_SIZE; k++) {
            bbuffer[k] = i_addr[0][i * BLOCK_SIZE + k];
        }
        for(uint32_t dd = 1; dd < n_digits; dd++) {
            for(uint32_t k = 0; k < BLOCK_SIZE; k++) {
                bbuffer[k] += i_addr[dd][i * BLOCK_SIZE + k];
            }
        }
        for(uint32_t k = 0; k < BLOCK_SIZE; k++) {
            outputs[i * BLOCK_SIZE + k] += bbuffer[k];
        }
    }
}

#define SWITCH_CASE(i) case i: MergeDigits<i>(output, dgbuf, n_elem); break;
void GenericSwitch(uint64_t* __restrict output, UINTNAT** dgbuf, uint32_t n_digits, uint64_t n_elem) {

    assert(n_elem % BLOCK_SIZE == 0);
    n_elem = BLOCK_SIZE * (n_elem / BLOCK_SIZE);
    // Magic

        switch (n_digits) {
            SWITCH_CASE(2)
            SWITCH_CASE(3)
            SWITCH_CASE(4)
            SWITCH_CASE(5)
            SWITCH_CASE(6)
            SWITCH_CASE(7)
            SWITCH_CASE(8)
        }

}

RLWECiphertext RKSKDigitPrecompFast::DoKeyswitch(const std::vector<LWECiphertext> &cts) {

    auto repeat = m_max_ct / cts.size();

    auto PolyParams = m_params->GetPolyParams();
    auto Q = m_params->GetQ().ConvertToInt<UINTNAT>();
    auto N = m_params->GetN();
    auto factor = N / m_max_ct;


    NativePoly DistrPoly = NativePoly(PolyParams, COEFFICIENT, true);

    for(uint32_t i = 0; i < repeat; i++) {
        DistrPoly[i * factor] = 1;
    }
    DistrPoly.SetFormat(EVALUATION);


    auto a_poly = NativePoly(PolyParams, EVALUATION, true);
    auto b_poly = NativePoly(PolyParams, EVALUATION, true);

    NativePoly PadPoly = NativePoly(PolyParams, COEFFICIENT, true);
    for(uint32_t j = 0; j < factor; j++) {
        PadPoly[j] = 1;
    }
    PadPoly.SetFormat(EVALUATION);

    UINTNAT mask = m_L - 1;
    uint64_t block_offset = 2 * N * m_digits * (m_L - 1);

    UINTNAT* dgbuf[MAX_DIGITS];
    uint64_t ndig = 0;

    std::vector<uint64_t> intermediate_B;
    std::memset(m_out_buffer, 0, m_max_ct * 2 * N * sizeof(uint64_t));

    for(uint32_t cti = 0; cti < cts.size(); cti++) {
        auto* current_out = m_out_buffer + cti * 2 * N;
        auto& ctA = cts[cti]->GetA();
        uint64_t B = 0;

        if (m_use_intermediate) {
            std::memset(m_LWE_out_buffer, 0, N * sizeof(uint64_t));
            // First we switch to intermediate LWE
            auto lwe_block_offset = m_LWE_n * m_LWE_digits * (m_LWE_L - 1);
            auto LWE_mask = m_LWE_L - 1;

            //B = cts[cti]->GetB().ConvertToInt<uint64_t>();

            for(uint32_t i = 0; i < N; i++) {
                auto di = ctA[i].ConvertToInt<UINTNAT>();
                auto* chunk = m_LWE_ksk_a + i * lwe_block_offset;

                ndig = 0;
                for(uint32_t j = 0; j < MAX_DIGITS; j++) {
                    auto dig = di & LWE_mask;
                    if (dig) {
                        dgbuf[ndig++] = chunk + (j * (m_LWE_L - 1) + (dig - 1)) * m_LWE_n;
                        B += m_LWE_ksk_b[(dig - 1) + j * (m_LWE_L - 1) + i * m_LWE_digits * (m_LWE_L - 1)];
                    }
                    di >>= m_LWE_L_bits;
                }

                if (ndig == 1) {
                    for(uint32_t k = 0; k < m_LWE_n / BLOCK_SIZE; k++) {
                        for(uint32_t kk = 0; kk < BLOCK_SIZE; kk++) {
                            m_LWE_out_buffer[kk + k * BLOCK_SIZE] += dgbuf[0][kk + k * BLOCK_SIZE];
                        }
                    }
                }

                if (ndig > 1) {
                    GenericSwitch(m_LWE_out_buffer, dgbuf, ndig, m_LWE_n);
                }
            }

            auto Bneg = (Q - (B % Q)) % Q;
            B = (cts[cti]->GetB().ConvertToInt<UINTNAT>() + Bneg) % Q;
            // negate and reduce
            for(uint32_t i = 0; i < m_LWE_n; i++) {
                m_LWE_out_buffer[i] = (Q - (m_LWE_out_buffer[i] % Q)) % Q;
            }

        } else {
            for(uint32_t kk = 0; kk < N; kk++)
                m_LWE_out_buffer[kk] = ctA[kk].ConvertToInt<UINTNAT>();
            B = cts[cti]->GetB().ConvertToInt<UINTNAT>();
        }

        intermediate_B.push_back(B);

        auto next_dim = m_use_intermediate ? m_LWE_n : N;
        for(uint32_t i = 0; i < next_dim; i++) {
            auto di = m_LWE_out_buffer[i];
            auto* chunk = m_ksk + i * block_offset;
            ndig = 0;

            // decompose first
            for(uint32_t j = 0; j < MAX_DIGITS; j++) {
                auto dig = di & mask;
                if (dig) {
                    dgbuf[ndig++] = chunk + (j * (m_L - 1) + (dig - 1)) * 2 * N;
                }
                di >>= m_L_bits;
            }

            // we save a call / jump
            if (ndig == 1) {
                for(uint32_t k = 0; k < 2 * N / BLOCK_SIZE; k++) {
                    for(uint32_t kk = 0; kk < BLOCK_SIZE; kk++) {
                        current_out[kk + k * BLOCK_SIZE] += dgbuf[0][kk + k * BLOCK_SIZE];
                    }
                }
            }

            if (ndig > 1) {
                GenericSwitch(current_out, dgbuf, ndig, 2 * N);
            }
        }

    }

    // rebuild
    NativePoly acc_a(PolyParams, EVALUATION, true);
    NativePoly acc_b(PolyParams, EVALUATION, true);
    for(uint32_t i = 0; i < cts.size(); i++) {
        auto col = m_out_buffer + i * 2 * N;

        NativeInteger bi = intermediate_B[i] % Q;
        auto tmp = m_offset_polys[i * repeat] * DistrPoly;
        tmp *= PadPoly;
        acc_b += (tmp *= bi);

        for(uint32_t j = 0; j < N; j++) {
            a_poly[j] = col[j] % Q;
            b_poly[j] = col[j + N] % Q;
        }

        a_poly *= m_offset_polys[i * repeat] * DistrPoly;
        b_poly *= m_offset_polys[i * repeat] * DistrPoly;


        acc_a -= a_poly;
        acc_b -= b_poly;

    }

    std::vector<NativePoly> entries = {acc_a, acc_b};

    return std::make_shared<RLWECiphertextImpl>(entries);

}

RKSKDigitPrecompFast::~RKSKDigitPrecompFast() {

    std::free(m_ksk);
    std::free(m_out_buffer);
    std::free(m_LWE_out_buffer);

    if (m_use_intermediate) {
        std::free(m_LWE_ksk_a);
        std::free(m_LWE_ksk_b);
    }

}
