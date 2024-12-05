
#include "RLWESwitchingMethods/RKSKDecimationFast.h"

#define BLOCK_SIZE 64
#define MAX_DIGITS 8
#define ALIGNMENT 64

RKSKDecimationFast::RKSKDecimationFast(std::shared_ptr<RingGSWCryptoParams> params, uint64_t L_bits, uint32_t max_ct) : m_params(std::move(params)),
                                                                                                                m_max_ct(max_ct), m_L(1ull << L_bits), m_L_bits(L_bits)
{

}

void RKSKDecimationFast::Keygen(NativeVector &sk_in, lbcrypto::NativePoly &sk_out) {

    (void) sk_in;

    // For debugging purposes, do not deploy !
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


    uint64_t Qbits = 128 - clz_u128(Q) - 1;
    uint64_t Nbits = 128 - clz_u128((uint64_t)N) - 1;
    m_digits = Qbits / m_L_bits + (Qbits % m_L_bits != 0);

    DiscreteUniformGeneratorImpl<NativeVector> dug;
    dug.SetModulus(Q);

    m_auto_keys_raw = static_cast<uint64_t*>(aligned_alloc(ALIGNMENT, Nbits * m_digits * 2 * N * sizeof(uint64_t)));
    m_output_buffer = static_cast<uint64_t*>(aligned_alloc(ALIGNMENT, m_max_ct * 2 * N * sizeof(uint64_t)));
    m_decomp_buffer = static_cast<uint64_t*>(aligned_alloc(ALIGNMENT, m_max_ct * m_digits * N * sizeof(uint64_t)));

    // Generate automorphism keys
    std::vector<uint64_t> tmp_buffer(2 * N);
    for(int aut = 1; aut < N; aut*=2) {
        uint32_t auto_idx = N / aut + 1;
        NativePoly sk_t = sk_out.AutomorphismTransform(auto_idx);
        NativeInteger factor = 1;

        std::vector<RLWECiphertext> automorphism_key_aut;
        for(int i = 0; i < m_digits; i++) {
            //auto A = NativePoly(dug,PolyParams, EVALUATION);
            auto A = NativePoly(PolyParams, EVALUATION, true);
            auto B = sk_t * factor;
            B += A * sk_out;

            auto err = NativePoly(dgg, PolyParams, EVALUATION);
            //B += err;

            for(int j = 0; j < N; j++) {
                tmp_buffer[j] = A[j].ConvertToInt<uint64_t>();
                tmp_buffer[j+N] = B[j].ConvertToInt<uint64_t>();
            }

            std::memcpy(m_auto_keys_raw + aut * m_digits * 2 * N + i * 2 * N, tmp_buffer.data(), 2 * N * sizeof(uint64_t));

            factor.ModMulEq(m_L, Q);
        }
        m_auto_keys.push_back(automorphism_key_aut);
    }

}

template<uint32_t n_digits> void DecomposeMultipleVectorsProto(uint64_t* __restrict output_buffer, 
        const uint64_t* __restrict input_buffer, uint32_t N, uint32_t n_cts, uint64_t shift, uint64_t mask) {
    for (uint32_t j = 0; j < n_digits; j++) {
        auto digit_block = output_buffer + j * n_cts * N;
        for(uint32_t ct_idx = 0; ct_idx < n_cts; ct_idx++) {
            auto ct_row = digit_block + ct_idx * N;
            
            for (uint32_t i = 0; i < N; i++) {
                ct_row[i] = (input_buffer[ct_idx + i * n_cts * N] >> shift) & mask;
            }
        }
    }
}

void DecomposeMultipleVectors(uint64_t* __restrict output_buffer,
                              const uint64_t* __restrict input_buffer, uint32_t N, uint32_t n_cts, uint32_t digits, uint64_t shift, uint64_t mask) {
    switch (digits) {
        case 2:
            DecomposeMultipleVectorsProto<2>(output_buffer, input_buffer, N, n_cts, shift, mask);
            break;
        case 3:
            DecomposeMultipleVectorsProto<3>(output_buffer, input_buffer, N, n_cts, shift, mask);
            break;
        case 4:
            DecomposeMultipleVectorsProto<4>(output_buffer, input_buffer, N, n_cts, shift, mask);
            break;
        case 5:
            DecomposeMultipleVectorsProto<5>(output_buffer, input_buffer, N, n_cts, shift, mask);
            break;
        case 6:
            DecomposeMultipleVectorsProto<6>(output_buffer, input_buffer, N, n_cts, shift, mask);
            break;
        case 7:
            DecomposeMultipleVectorsProto<7>(output_buffer, input_buffer, N, n_cts, shift, mask);
            break;
        case 8:
            DecomposeMultipleVectorsProto<8>(output_buffer, input_buffer, N, n_cts, shift, mask);
            break;
        case 9:
            DecomposeMultipleVectorsProto<9>(output_buffer, input_buffer, N, n_cts, shift, mask);
            break;
        default:
            std::cerr << "No matching decomp case..." << std::endl;
            return;
    }
}    


RLWECiphertext RKSKDecimationFast::DoKeyswitch(const std::vector<LWECiphertext> &cts) {

    auto PolyParams = m_params->GetPolyParams();
    auto Q = m_params->GetQ().ConvertToInt<UINTNAT>();
    auto N = m_params->GetN();
    auto factor = N / m_max_ct;
    auto mask = m_L - 1;

    (void)mask;

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
    PadPoly *= DistrPoly;

    NativeInteger scale = N;
    scale.ModExpEq(Q - 2, Q);

    // Note, do not use ModInverse, it does not return the correct inverse
    // So here, we simply use small fermat i.e. a^{p-2} = a^{-1} mod p when p is prime
    // start by recovering the "true" value of the "a" polynomial

    auto ct_stride = cts.size();


    return std::make_shared<RLWECiphertextImpl>(std::vector<NativePoly>(ct_stride));
}

RKSKDecimationFast::~RKSKDecimationFast() {
    free(m_auto_keys_raw);
    free(m_output_buffer);
    free(m_decomp_buffer);
}
