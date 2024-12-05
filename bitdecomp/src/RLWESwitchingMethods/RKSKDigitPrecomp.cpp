
#include "RLWESwitchingMethods/RKSKDigitPrecomp.h"


RKSKDigitPrecomp::RKSKDigitPrecomp(std::shared_ptr<RingGSWCryptoParams> params, uint64_t L_bits, uint32_t max_ct) :
        m_params(std::move(params)), m_max_ct(max_ct), m_L(1ull << L_bits), m_L_bits(L_bits)
{
}

RKSKDigitPrecomp::RKSKDigitPrecomp(std::shared_ptr<RingGSWCryptoParams> params, uint64_t L_bits,
                                           uint64_t L_LWE_bits, uint64_t n_LWE, uint32_t max_ct) :
        m_params(std::move(params)), m_max_ct(max_ct), m_L(1ull << L_bits), m_L_bits(L_bits), m_LWE_n(n_LWE),
        m_LWE_L(1ull << L_LWE_bits), m_LWE_L_bits(L_LWE_bits), m_use_intermediate(true)

{
}

void RKSKDigitPrecomp::Keygen(NativeVector &sk_in, NativePoly &sk_out) {

    auto Qnat = m_params->GetQ();
    auto Q = Qnat.ConvertToInt<UINTNAT>();
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
                for(uint32_t k = 1; k < m_LWE_L; k++) {
                    auto skijk = ski.ModMul(k, Qnat);
                    auto vec_a = dug.GenerateVector(m_LWE_n);
                    // noiseless
                    auto vec_b = skijk.ModAdd(0, Qnat);
                    // with noise
                    // auto vec_b = skijk.ModAdd(dgg.GenerateInteger(Q), Qnat);

                    // encrypt and store the A component at the same time
                    for(uint32_t kk = 0; kk < m_LWE_n; kk++) {
                        auto askk = vec_a[kk].ModMul(intermediate_key[kk], Q);
                        vec_b.ModAddEq(askk, Q);
                    }
                    digits.push_back(std::make_shared<LWECiphertextImpl>(vec_a, vec_b));
                }
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

    for(uint32_t i = 0; i < int_N; i++) {
        NativeInteger ski = intermediate_key[i];
        NativeInteger factor = 1;
        /* go other every power of the basis */
        std::vector<std::vector<RLWECiphertext>> rlwe_power_digits;
        for(uint32_t j = 0; j < m_digits; j++) {
            auto skij = ski.ModMul(factor, Q);

            std::vector<RLWECiphertext> rlwe_digits;
            for(uint32_t k = 1; k < m_L; k++) {
                auto skijk = skij.ModMul(k, Q);

                NativePoly a_poly(dug, PolyParams, EVALUATION);

                // Noiseless
                NativePoly b_poly(PolyParams, COEFFICIENT, true);

                // With noise
                //NativePoly b_poly(dgg, PolyParams, COEFFICIENT);

                b_poly.SetFormat(EVALUATION);
                b_poly += skijk * Xfi;

                b_poly += a_poly * sk_out;
                std::vector<NativePoly> entries = {a_poly, b_poly};
                rlwe_digits.push_back(std::make_shared<RLWECiphertextImpl>(entries));
            }
            rlwe_power_digits.push_back(rlwe_digits);
            factor.MulEq(m_L);
        }
        m_ksk.push_back(rlwe_power_digits);
    }

}

RLWECiphertext RKSKDigitPrecomp::DoKeyswitch(const std::vector<LWECiphertext> &cts_in) {

    auto repeat = m_max_ct / cts_in.size();
    std::vector<LWECiphertext> cts;
    for(auto& ct : cts_in) {
        for(uint32_t i = 0; i < repeat; i++)
            cts.push_back(ct);
    }

    auto Q = m_params->GetQ().ConvertToInt<UINTNAT>();
    auto N = m_params->GetN();
    auto factor = N / m_max_ct;

    auto PolyParams = m_params->GetPolyParams();
    auto a_poly = NativePoly(PolyParams, EVALUATION, true);
    auto b_poly = NativePoly(PolyParams, EVALUATION, true);

    NativePoly PadPoly = NativePoly(PolyParams, COEFFICIENT, true);
    for(uint32_t j = 0; j < factor; j++) {
        PadPoly[j] = 1;
    }
    PadPoly.SetFormat(EVALUATION);

    UINTNAT mask = m_L - 1;

    auto next_dim = m_use_intermediate ? m_LWE_n : N;
    NativeVector LWE_A = NativeVector(next_dim, Q);
    NativeInteger LWE_B;

    NativePoly acc_a(PolyParams, EVALUATION, true);
    NativePoly acc_b(PolyParams, EVALUATION, true);

    for(uint32_t cti = 0; cti < cts.size(); cti++) {
        auto& ctA = cts[cti]->GetA();

        if (m_use_intermediate) {
            LWE_A *= 0;
            LWE_B = cts[cti]->GetB();

            auto LWE_mask = m_LWE_L - 1;

            for(uint32_t i = 0; i < N; i++) {
                auto& i_entries = m_LWE_ksk[i];

                auto di = ctA[i].ConvertToInt<UINTNAT>();

                for(uint32_t j = 0; j < m_digits; j++) {
                    auto dig = di & LWE_mask;
                    if (dig) {
                        auto& lwe_ijk = i_entries[j][dig - 1];
                        auto& lwe_ijk_A = lwe_ijk->GetA();
                        auto& lwe_ijk_B = lwe_ijk->GetB();

                        LWE_A.ModSubEq(lwe_ijk_A);
                        LWE_B.ModSubEq(lwe_ijk_B, Q);
                    }
                    di >>= m_LWE_L_bits;
                }

            }

        } else {
            LWE_A = cts[cti]->GetA();
            LWE_B = cts[cti]->GetB();
        }

        for(uint64_t i = 0; i < next_dim; i++) {
            auto di = LWE_A[i].ConvertToInt<UINTNAT>();
            auto& i_rlwe_entries = m_ksk[i];
            for(uint32_t j = 0; j < m_digits; j++) {
                auto dig = di & mask;
                if (dig) {
                    auto& rlwe_ijk = i_rlwe_entries[j][dig - 1];
                    auto& A = rlwe_ijk->GetElements()[0];
                    auto& B = rlwe_ijk->GetElements()[1];

                    a_poly += A;
                    b_poly += B;
                }
                di >>= m_L_bits;
            }

        }

        acc_a -= a_poly *= m_offset_polys[cti] ;
        acc_b += (LWE_B * PadPoly - b_poly) * m_offset_polys[cti];

        a_poly.SetValuesToZero();
        b_poly.SetValuesToZero();
    }

    std::vector<NativePoly> entries = {acc_a, acc_b};

    return std::make_shared<RLWECiphertextImpl>(entries);

}
