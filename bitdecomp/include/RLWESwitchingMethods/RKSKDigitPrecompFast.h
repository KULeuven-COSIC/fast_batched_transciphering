
#ifndef CPP_IMPL_RKSKDIGITPRECOMPFAST_H
#define CPP_IMPL_RKSKDIGITPRECOMPFAST_H

#include "RLWEKeyswitchInterface.h"

class RKSKDigitPrecompFast : public RKSKInterface {

public:

    RKSKDigitPrecompFast(std::shared_ptr<RingGSWCryptoParams> params, uint64_t L_bits, uint64_t L_LWE_bits, uint64_t n_LWE, uint32_t max_ct);

    RKSKDigitPrecompFast(std::shared_ptr<RingGSWCryptoParams> params, uint64_t L_bits, uint32_t max_ct);

    void Keygen(NativeVector &sk_in, lbcrypto::NativePoly &sk_out) override;

    RLWECiphertext DoKeyswitch(const std::vector<LWECiphertext> &ct_in) override;

    ~RKSKDigitPrecompFast() override;

private:

    std::shared_ptr<RingGSWCryptoParams> m_params;

    uint32_t m_max_ct;

    uint64_t m_L;
    uint64_t m_L_bits;
    uint64_t m_digits;

    uint64_t m_LWE_n = 0;
    uint64_t m_LWE_L = 0;
    uint64_t m_LWE_L_bits = 0;
    uint64_t m_LWE_digits = 0;

    UINTNAT* m_ksk;
    uint64_t* m_out_buffer;

    UINTNAT* m_LWE_ksk_a = nullptr;
    UINTNAT* m_LWE_ksk_b = nullptr;
    uint64_t* m_LWE_out_buffer = nullptr;

    bool m_use_intermediate = false;

    std::vector<NativePoly> m_offset_polys;

};

#endif //CPP_IMPL_RKSKDIGITPRECOMPFAST_H
