//
// Created by ******** on 3/21/24.
// We implement the iterative version of Algorithm 2 in https://eprint.iacr.org/2020/015.pdf
//
#include "RLWESwitchingMethods/RKSKDecimation.h"

#include <utility>

RKSKDecimation::RKSKDecimation(std::shared_ptr<RingGSWCryptoParams> params, uint64_t L_bits, uint32_t max_ct) : m_params(std::move(params)),
        m_max_ct(max_ct), m_L(1ull << L_bits), m_L_bits(L_bits)
{

}

void RKSKDecimation::Keygen(NativeVector &sk_in, lbcrypto::NativePoly &sk_out) {
    (void) sk_in;
    (void) m_max_ct;

    sk = sk_out;
    sk.SetFormat(EVALUATION);

    auto Qnat = m_params->GetQ();
    auto Q = Qnat.ConvertToInt<uint64_t>();
    auto N = m_params->GetN();
    auto dgg = m_params->GetDgg();
    auto PolyParams = m_params->GetPolyParams();
    sk_out.SetFormat(EVALUATION);

    auto n_pad = m_params->GetN() / m_max_ct;
    for(uint32_t i = 0; i < m_max_ct; i++) {
        auto Xfi = NativePoly(PolyParams, COEFFICIENT, true);
        Xfi[i * n_pad] = 1;
        Xfi.SetFormat(EVALUATION);
        m_offset_polys.push_back(Xfi);
    }

    auto Qbits = 0ull;
    while ((1ull << Qbits) < Q){
        Qbits++;
    }

    m_digits = Qbits / m_L_bits + (Qbits % m_L_bits != 0);

    DiscreteUniformGeneratorImpl<NativeVector> dug;
    dug.SetModulus(Q);

    NativeInteger I2 = 2;
    I2.ModExpEq(Q - 2, Q);

    NativeInteger scale = N;
    scale.ModExpEq(Q - 2, Q);

    // Generate automorphism keys
    for(int aut = 1; aut < N; aut*=2) {
        uint32_t auto_idx = N / aut + 1;
        NativePoly sk_t = sk_out.AutomorphismTransform(auto_idx);
        NativeInteger factor = 1;

        std::vector<RLWECiphertext> automorphism_key_aut;
        for(int i = 0; i < m_digits; i++) {
            auto A = NativePoly(dug,PolyParams, EVALUATION);
            auto E = NativePoly(dgg, PolyParams, COEFFICIENT);
            E.SetFormat(EVALUATION);
            //auto A = NativePoly(PolyParams, EVALUATION, true);
            auto B = sk_t * factor;
            B += A * sk_out;
            B += E;

            std::vector<NativePoly> AB = {A, B};
            automorphism_key_aut.emplace_back(std::make_shared<RLWECiphertextImpl>(AB));
            factor.ModMulEq(m_L, Q);
        }
        scale.ModMulEq(I2, Q);
        m_auto_keys.push_back(automorphism_key_aut);
    }

}
uint32_t reverse32(uint32_t x)
{
    // See: https://stackoverflow.com/questions/9144800/c-reverse-bits-in-unsigned-integer
    x = ((x >> 1) & 0x55555555u) | ((x & 0x55555555u) << 1);
    x = ((x >> 2) & 0x33333333u) | ((x & 0x33333333u) << 2);
    x = ((x >> 4) & 0x0f0f0f0fu) | ((x & 0x0f0f0f0fu) << 4);
    x = ((x >> 8) & 0x00ff00ffu) | ((x & 0x00ff00ffu) << 8);
    x = ((x >> 16) & 0xffffu) | ((x & 0xffffu) << 16);
    return x;
}

RLWECiphertext RKSKDecimation::DoKeyswitch(const std::vector<LWECiphertext> &cts) {

    auto PolyParams = m_params->GetPolyParams();
    auto Q = m_params->GetQ().ConvertToInt<UINTNAT>();
    (void) Q;
    auto N = m_params->GetN();
    auto factor = N / cts.size();

    auto acc_a = NativePoly(PolyParams, EVALUATION, true);
    auto acc_b = NativePoly(PolyParams, EVALUATION, true);

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
    //PadPoly *= DistrPoly;

    NativeInteger scale = N;
    scale.ModExpEq(Q - 2, Q);
    int32_t first_stage_depth = 0;
    NativePoly A_red;
    NativePoly B_red;

    if (cts.size() > 1) {

        std::vector<LWECiphertext> cts_brev;
        /* DANGER: We assume that the number of samples is always a power of 2 */
        uint32_t ct_size_u32 = cts.size();
        uint32_t log2_ct_size = 0;
        while (ct_size_u32 >>= 1) ++log2_ct_size;

        for(uint32_t i = 0; i < cts.size(); i++) {
            auto rev_i_32 = reverse32(i);
            auto rev_i = rev_i_32 >> (32 - log2_ct_size);
            cts_brev.push_back(cts[rev_i]);
        }

        first_stage_depth = std::round(std::log2(cts.size()));
        // setup
        std::vector<NativePoly> reduced_rlwe;
        for(auto& ct : cts_brev) {

            //std::cout << scale << std::endl;
            NativeVector lwe_a = ct->GetA();
            NativeInteger lwe_b = ct->GetB();

            // first we recover the "true" (a) polynomial
            auto true_a = NativePoly(PolyParams, COEFFICIENT, true);

            true_a.SetValues(lwe_a, COEFFICIENT);
            true_a.SetFormat(EVALUATION);
            true_a = true_a.Transpose();
            true_a *= scale;

            auto true_b = NativePoly(PolyParams, COEFFICIENT, true);
            true_b[0] = lwe_b;
            true_b.SetFormat(EVALUATION);
            true_b *= scale;

            reduced_rlwe.push_back(true_a);
            reduced_rlwe.push_back(true_b);
        }

        for(int32_t i = 0; i < first_stage_depth; i++) {

            NativePoly shift = NativePoly(PolyParams, COEFFICIENT, true);
            shift[N >> (i + 1)] = 1;
            shift.SetFormat(EVALUATION);

            std::vector<NativePoly> reduced_new;
            std::vector<NativePoly> to_switch;

            for(int32_t j = 0; j < reduced_rlwe.size(); j+=4) {

                auto& A0 = reduced_rlwe[j];
                auto& B0 = reduced_rlwe[j+1];

                auto& A1 = (reduced_rlwe[j+2] *= shift);
                auto& B1 = (reduced_rlwe[j+3] *= shift);

                auto tmp0 = reduced_rlwe[j] - reduced_rlwe[j+2];
                auto tmp1 = reduced_rlwe[j+1] - reduced_rlwe[j+2+1];

                auto auto_a = tmp0.AutomorphismTransform((1 << (i + 1)) + 1);
                auto auto_b = tmp1.AutomorphismTransform((1 << (i + 1)) + 1);

                auto autom2 = DoRLWEToRLWESwitch(m_auto_keys[m_auto_keys.size() - i - 1], auto_a, auto_b);

                auto A2p = (A0 += A1) += autom2->GetElements()[0];
                auto B2p = (B0 += B1) += autom2->GetElements()[1];

                reduced_new.push_back(A2p);
                reduced_new.push_back(B2p);
            }

            reduced_rlwe.clear();
            reduced_rlwe = reduced_new;

        }

        A_red = reduced_rlwe[0];
        B_red = reduced_rlwe[1];
    } else {
        NativeVector lwe_a = cts[0]->GetA();
        NativeInteger lwe_b = cts[0]->GetB();

        // first we recover the "true" (a) polynomial
        auto true_a = NativePoly(PolyParams, COEFFICIENT, true);

        true_a.SetValues(lwe_a, COEFFICIENT);
        true_a.SetFormat(EVALUATION);
        true_a = true_a.Transpose();
        true_a *= scale;

        auto true_b = NativePoly(PolyParams, COEFFICIENT, true);
        true_b[0] = lwe_b;
        true_b.SetFormat(EVALUATION);
        true_b *= scale;

        A_red = true_a;
        B_red = true_b;
    }

    for(int32_t i = 0; i < (m_auto_keys.size() - first_stage_depth); i++) {
        auto auto_a = A_red.AutomorphismTransform((N >> i) + 1);
        auto auto_b = B_red.AutomorphismTransform((N >> i) + 1);

        auto rlwe_switched = DoRLWEToRLWESwitch(m_auto_keys[i], auto_a, auto_b);

        A_red += rlwe_switched->GetElements()[0];
        B_red += rlwe_switched->GetElements()[1];
    }

    /*
    for(uint32_t cti = 0; cti < cts.size(); cti++) {

        // Note, do not use ModInverse, it does not return the correct inverse
        // So here, we simply use small fermat i.e. a^{p-2} = a^{-1} mod p when p is prime
        //std::cout << scale << std::endl;
        NativeVector lwe_a = cts[cti]->GetA();
        NativeInteger lwe_b = cts[cti]->GetB();

        // first we recover the "true" (a) polynomial
        auto true_a = NativePoly(PolyParams, COEFFICIENT, true);

        true_a.SetValues(lwe_a, COEFFICIENT);
        true_a.SetFormat(EVALUATION);
        true_a = true_a.Transpose();
        true_a *= scale;

        auto true_b = NativePoly(PolyParams, COEFFICIENT, true);
        true_b[0] = lwe_b;
        true_b.SetFormat(EVALUATION);
        true_b *= scale;

        // next step remove coefficients via decimation
        for(int i = 0; i < m_auto_keys.size(); i++) {
            auto auto_a = true_a.AutomorphismTransform((N >> i) + 1);
            auto auto_b = true_b.AutomorphismTransform((N >> i) + 1);

            auto rlwe_switched = DoRLWEToRLWESwitch(m_auto_keys[i], auto_a, auto_b);

            true_a += rlwe_switched->GetElements()[0];
            true_b += rlwe_switched->GetElements()[1];
        }

        true_a *= PadPoly;
        true_b *= PadPoly;

        acc_a += (true_a *= m_offset_polys[cti * repeat]);
        acc_b += (true_b *= m_offset_polys[cti * repeat]);
    }
    */

    if (factor != 1) {
        A_red *= PadPoly;
        B_red *= PadPoly;
    }

    std::vector<NativePoly> coefs = {A_red, B_red};
    return std::make_shared<RLWECiphertextImpl>(coefs);
}

RLWECiphertext
RKSKDecimation::DoRLWEToRLWESwitch(std::vector<RLWECiphertext> &ksk, NativePoly &A_in, NativePoly &B_in) {

    NativePoly b_out = B_in;
    b_out.SetFormat(EVALUATION);

    NativePoly a_out = NativePoly(m_params->GetPolyParams(), EVALUATION, true);

    NativePoly a_in =A_in;
    a_in.SetFormat(COEFFICIENT);

    NativePoly digit_poly;


    auto mask = m_L - 1;

    for(int i = 0; i< m_digits; i++) {

        digit_poly = NativePoly(m_params->GetPolyParams(), COEFFICIENT, true);
        for(int j = 0; j < m_params->GetN(); j++) {
            auto aj = a_in[j].ConvertToInt<uint64_t>();
            digit_poly[j] = aj & mask;
            a_in[j] = aj >> m_L_bits;
        }
        digit_poly.SetFormat(EVALUATION);
        auto d1 = digit_poly;


        auto ksk_a = ksk[i]->GetElements()[0];
        auto ksk_b = ksk[i]->GetElements()[1];

        b_out -= (digit_poly *= ksk_b);
        a_out -= (d1 *= ksk_a);
    }


    std::vector<NativePoly> coefs = {a_out, b_out};
    return std::make_shared<RLWECiphertextImpl>(coefs);
}
