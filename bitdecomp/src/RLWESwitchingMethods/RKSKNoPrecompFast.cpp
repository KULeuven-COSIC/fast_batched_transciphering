

#include "RLWESwitchingMethods/RKSKNoPrecompFast.h"

#define BLOCK_SIZE 64
#define MAX_DIGITS 8
#define ALIGNMENT 64

alignas(ALIGNMENT) uint64_t bbuffer[BLOCK_SIZE];

// Absolutely cursed, do not use for anything serious
template <class T>
class vectorView : public std::vector<T>
{
public:
    vectorView() {
        this->_M_impl._M_start = this->_M_impl._M_finish = this->_M_impl._M_end_of_storage = NULL;
    }

    vectorView(T* sourceArray, int arraySize)
    {
        this->_M_impl._M_start = sourceArray;
        this->_M_impl._M_finish = this->_M_impl._M_end_of_storage = sourceArray + arraySize;
    }

    ~vectorView() {
        this->_M_impl._M_start = this->_M_impl._M_finish = this->_M_impl._M_end_of_storage = NULL;
    }

    void wrapArray(T* sourceArray, int arraySize)
    {
        this->_M_impl._M_start = sourceArray;
        this->_M_impl._M_finish = this->_M_impl._M_end_of_storage = sourceArray + arraySize;
    }
};



RKSKNoPrecompFast::RKSKNoPrecompFast(std::shared_ptr<RingGSWCryptoParams> params, uint64_t L_bits, uint32_t max_ct) :
        m_params(std::move(params)), m_max_ct(max_ct), m_L(1ull << L_bits), m_L_bits(L_bits)
{
}

RKSKNoPrecompFast::RKSKNoPrecompFast(std::shared_ptr<RingGSWCryptoParams> params, uint64_t L_bits,
                                   uint64_t L_LWE_bits, uint64_t n_LWE, uint32_t max_ct) :
        m_params(std::move(params)), m_max_ct(max_ct), m_L(1ull << L_bits), m_L_bits(L_bits), m_LWE_n(n_LWE),
        m_LWE_L(1ull << L_LWE_bits), m_LWE_L_bits(L_LWE_bits), m_use_intermediate(true)

{
}

void RKSKNoPrecompFast::Keygen(NativeVector &sk_in, NativePoly &sk_out) {

    auto Qnat = m_params->GetQ();
    auto Q = Qnat.ConvertToInt<uint64_t>();
    auto N = m_params->GetN();
    auto dgg = m_params->GetDgg();
    auto PolyParams = m_params->GetPolyParams();

    sk_out.SetFormat(EVALUATION);
    NativeVector intermediate_key;

    uint64_t Qbits = 128 - clz_u128(Q) - 1;
    m_digits = Qbits / m_L_bits + (Qbits % m_L_bits != 0);

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

        // Generate first ksk
        for(uint32_t i = 0; i < N; i++) {
            std::vector<std::vector<LWECiphertext>> power_digits;
            auto ski = sk_in[i];
            for(uint32_t j = 0; j < m_LWE_digits; j++) {
                std::vector<LWECiphertext> digits;
                auto vec_a = dug.GenerateVector(m_LWE_n);
                    // noiseless
                auto vec_b = ski.ModAdd(0, Qnat);
                    // with noise
                    // auto vec_b = ski.ModAdd(dgg.GenerateInteger(Q), Qnat);

                    // encrypt and store the A component at the same time
                for(uint32_t kk = 0; kk < m_LWE_n; kk++) {
                    auto askk = vec_a[kk].ModMul(intermediate_key[kk], Q);
                    vec_b.ModAddEq(askk, Q);
                }
                digits.push_back(std::make_shared<LWECiphertextImpl>(vec_a, vec_b));
                power_digits.push_back(digits);
                ski.ModMulEq(m_LWE_L, Qnat);
            }
            m_LWE_ksk.push_back(power_digits);
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

    m_ksk_raw = static_cast<uint64_t*>(aligned_alloc(ALIGNMENT, 2 * int_N * m_digits * int_N * sizeof(uint64_t))); // [2 * int_N * m_digits * int_N];

    for(uint32_t i = 0; i < int_N; i++) {
        NativeInteger ski = intermediate_key[i];
        NativeInteger factor = 1;
        /* go other every power of the basis */
        std::vector<std::vector<RLWECiphertext>> rlwe_power_digits;
        for(uint32_t j = 0; j < m_digits; j++) {
            auto skij = ski.ModMul(factor, Q);

            std::vector<RLWECiphertext> rlwe_digits;

            NativePoly a_poly(dug, PolyParams, EVALUATION);

            // Noiseless
            NativePoly b_poly(PolyParams, COEFFICIENT, true);

            // With noise
            //NativePoly b_poly(dgg, PolyParams, COEFFICIENT);

            b_poly.SetFormat(EVALUATION);
            b_poly += skij * Xfi;

            b_poly += a_poly * sk_out;
            for(int kk = 0; kk < int_N; kk+=BLOCK_SIZE) {

                auto start_a = m_ksk_raw + kk * int_N * m_digits;
                auto start_b = m_ksk_raw + (kk + int_N) * int_N * m_digits;

                auto inner_a = start_a + i * m_digits * BLOCK_SIZE;
                auto inner_b = start_b + i * m_digits * BLOCK_SIZE;

                auto write_a = inner_a + j * BLOCK_SIZE;
                auto write_b = inner_b + j * BLOCK_SIZE;

                for(int kkk = 0; kkk < BLOCK_SIZE; kkk++) {
                    write_a[kkk] = a_poly[kk + kkk].ConvertToInt<uint64_t>();
                    write_b[kkk] = b_poly[kk + kkk].ConvertToInt<uint64_t>();
                }
            }

            //std::vector<NativePoly> entries = {a_poly, b_poly};
            //rlwe_digits.push_back(std::make_shared<RLWECiphertextImpl>(entries));

            //rlwe_power_digits.push_back(rlwe_digits);
            factor.MulEq(m_L);
        }
        m_ksk.push_back(rlwe_power_digits);
    }

}

RLWECiphertext RKSKNoPrecompFast::DoKeyswitch(const std::vector<LWECiphertext> &cts) {



    auto Q = m_params->GetQ().ConvertToInt<uint64_t>();
    auto N = m_params->GetN();
    auto factor = N / m_max_ct;

    auto PolyParams = m_params->GetPolyParams();
    auto a_poly = NativePoly(PolyParams, EVALUATION, true);
    auto b_poly = NativePoly(PolyParams, EVALUATION, true);

    NativePoly DistrPoly = NativePoly(PolyParams, COEFFICIENT, true);

    auto repeat = m_max_ct / cts.size();
    for(uint32_t i = 0; i < repeat; i++) {
        DistrPoly[i * factor] = 1;
    }
    DistrPoly.SetFormat(EVALUATION);

    NativePoly PadPoly = NativePoly(PolyParams, COEFFICIENT, true);
    for (uint32_t j = 0; j < factor; j++) {
        PadPoly[j] = 1;
    }
    PadPoly.SetFormat(EVALUATION);

    uint64_t mask = m_L - 1;

    auto next_dim = m_use_intermediate ? m_LWE_n : N;
    NativeVector LWE_A = NativeVector(next_dim, Q);
    NativeInteger LWE_B;

    NativePoly acc_a(PolyParams, EVALUATION, true);
    NativePoly acc_b(PolyParams, EVALUATION, true);

    auto mu = NativeInteger(Q).ComputeMu();

    for (uint32_t cti = 0; cti < cts.size(); cti++) {
        auto &ctA = cts[cti]->GetA();

        if (m_use_intermediate) {
            LWE_A *= 0;
            LWE_B = cts[cti]->GetB();

            auto LWE_mask = m_LWE_L - 1;

            for (uint32_t i = 0; i < N; i++) {
                auto &i_entries = m_LWE_ksk[i];

                auto di = ctA[i].ConvertToInt<uint64_t>();

                for (uint32_t j = 0; j < m_digits; j++) {
                    auto dig = di & LWE_mask;
                    if (dig) {
                        auto &lwe_ijk = i_entries[j][0];
                        auto &lwe_ijk_A = lwe_ijk->GetA();
                        auto &lwe_ijk_B = lwe_ijk->GetB();

                        LWE_A.ModSubEq(lwe_ijk_A.ModMul(dig));
                        LWE_B.ModSubEq(lwe_ijk_B.ModMul(dig, Q, mu), Q, mu);
                    }
                    di >>= m_LWE_L_bits;
                }

            }

        } else {
            LWE_A = cts[cti]->GetA();
            LWE_B = cts[cti]->GetB();
        }

        std::vector<uint64_t> a_buffer(2 * N);
        std::vector<uint16_t> in_digits(next_dim * m_digits);
        for(uint64_t i = 0; i < next_dim; i++) {
            auto di = LWE_A[i].ConvertToInt<uint64_t>();
            for (uint64_t j = 0; j < m_digits; j++) {
                in_digits.push_back(di & mask);
                di >>= m_L_bits;
            }
        }

        a_poly.AutomorphismTransform(3);

        std::fill(a_buffer.begin(), a_buffer.end(), 0);

        for(uint64_t kk = 0; kk < 2 * N; kk+=BLOCK_SIZE) {
            std::memset(bbuffer, 0, sizeof(bbuffer));

            for(uint64_t i = 0; i < next_dim; i++) {
                for(uint64_t j = 0; j < m_digits; j++) {
                    auto dij = in_digits[i * m_digits + j];
                    auto target = m_ksk_raw + kk * N * m_digits + i * m_digits * BLOCK_SIZE + j * BLOCK_SIZE;
                        for (uint64_t kkk = 0; kkk < BLOCK_SIZE; kkk++) {
                            bbuffer[kkk] += dij * target[kkk];
                        }
                }
            }

            for(uint64_t kkk = 0; kkk < BLOCK_SIZE; kkk++) {
                a_buffer[kk + kkk] = bbuffer[kkk];
            }

        }

        for(int k = 0; k < N; k++) {
            a_poly[k] = a_buffer[k] % m_params->GetQ();
            b_poly[k] = a_buffer[k + N] % m_params->GetQ();
        }

        a_poly *= m_offset_polys[cti * repeat] * DistrPoly;
        b_poly *= m_offset_polys[cti * repeat] * DistrPoly;

        acc_b += LWE_B * PadPoly * DistrPoly * m_offset_polys[cti * repeat];

        acc_a -= a_poly;
        acc_b -= b_poly;

    }

    std::vector<NativePoly> entries = {acc_a, acc_b};

    return std::make_shared<RLWECiphertextImpl>(entries);
}