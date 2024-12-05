

#include "RLWESwitchingMethods/RKSKNoPrecomp.h"


RKSKNoPrecomp::RKSKNoPrecomp(std::shared_ptr<RingGSWCryptoParams> params, uint64_t L_bits, uint32_t max_ct) :
        m_params(std::move(params)), m_max_ct(max_ct), m_L(1ull << L_bits), m_L_bits(L_bits)
{
}

RKSKNoPrecomp::RKSKNoPrecomp(std::shared_ptr<RingGSWCryptoParams> params, uint64_t L_bits,
                                   uint64_t L_LWE_bits, uint64_t n_LWE, uint32_t max_ct) :
        m_params(std::move(params)), m_max_ct(max_ct), m_L(1ull << L_bits), m_L_bits(L_bits), m_LWE_n(n_LWE),
        m_LWE_L(1ull << L_LWE_bits), m_LWE_L_bits(L_LWE_bits), m_use_intermediate(true)

{
}

void RKSKNoPrecomp::Keygen(NativeVector &sk_in, NativePoly &sk_out) {

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

    std::vector<int64_t> temp_buffer(2 * int_N);
    for(uint32_t i = 0; i < int_N; i++) {
        NativeInteger ski = intermediate_key[i];
        NativeInteger factor = 1;
        /* go other every power of the basis */
        std::vector<std::vector<RLWECiphertext>> rlwe_power_digits;
        for(uint32_t j = 0; j < m_digits; j++) {
            auto skij = ski.ModMul(factor, Q);

            std::vector<RLWECiphertext> rlwe_digits;

            //NativePoly a_poly(dug, PolyParams, EVALUATION);
            NativePoly a_poly(PolyParams, EVALUATION, true);

            // Noiseless
            NativePoly b_poly(PolyParams, COEFFICIENT, true);

            // With noise
            //NativePoly b_poly(dgg, PolyParams, COEFFICIENT);

            b_poly.SetFormat(EVALUATION);
            b_poly += skij;

            b_poly += a_poly * sk_out;

            std::vector<NativePoly> entries = {a_poly, b_poly};
            rlwe_digits.push_back(std::make_shared<RLWECiphertextImpl>(entries));

            rlwe_power_digits.push_back(rlwe_digits);
            factor.MulEq(m_L);
        }
        m_ksk.push_back(rlwe_power_digits);
    }

}

RLWECiphertext RKSKNoPrecomp::DoKeyswitch(const std::vector<LWECiphertext> &cts_in) {

    auto repeat = m_max_ct / cts_in.size();
    std::vector<LWECiphertext> cts;
    for (auto &ct: cts_in) {
        for (uint32_t i = 0; i < repeat; i++)
            cts.push_back(ct);
    }

    auto Q = m_params->GetQ().ConvertToInt<uint64_t>();
    auto N = m_params->GetN();
    auto factor = N / cts.size();

    auto PolyParams = m_params->GetPolyParams();
    auto a_poly = NativePoly(PolyParams, EVALUATION, true);
    auto b_poly = NativePoly(PolyParams, EVALUATION, true);

    NativePoly PadPoly = NativePoly(PolyParams, COEFFICIENT, true);
    for (uint32_t j = 0; j < factor; j++) {
        PadPoly[j] = 1;
    }
    PadPoly.SetFormat(EVALUATION);

    uint64_t mask = m_L - 1;

    auto next_dim = m_use_intermediate ? m_LWE_n : N;
    NativeVector LWE_A = NativeVector(next_dim, Q);
    NativeInteger LWE_B;

    (void) LWE_B;

    NativePoly acc_a(PolyParams, EVALUATION, true);
    NativePoly acc_b(PolyParams, COEFFICIENT, true);

    auto mu = NativeInteger(Q).ComputeMu();
    (void) mu;

    for(uint64_t kk = 0; kk < cts.size(); kk++) {
        acc_b[kk * factor] = cts[kk]->GetB();
    }
    acc_b.SetFormat(EVALUATION);

    for (uint64_t i = 0; i < next_dim; i++) {
        auto &i_rlwe_entries = m_ksk[i];

        for (uint32_t j = 0; j < m_digits; j++) {
            NativePoly digit_poly = NativePoly(PolyParams, COEFFICIENT, true);
            for(uint64_t kk = 0; kk < cts.size(); kk++) {
                auto di = cts[kk]->GetA()[i].ConvertToInt<uint64_t>();
                auto digK = ((di >> (m_L_bits * j)) & mask);
                digit_poly[kk * factor] = digK;
            }
            digit_poly.SetFormat(EVALUATION);

            //if (dig) {
                auto &rlwe_ijk = i_rlwe_entries[j][0];
                auto A = rlwe_ijk->GetElements()[0];
                auto B = rlwe_ijk->GetElements()[1];

                acc_a -= digit_poly * (A);
                acc_b -= digit_poly * (B);
            //}
        }

    }

    if (factor > 0) {
        acc_a *= PadPoly;
        acc_b *= PadPoly;
    }

    std::vector<NativePoly> entries = {acc_a, acc_b};

    return std::make_shared<RLWECiphertextImpl>(entries);
}