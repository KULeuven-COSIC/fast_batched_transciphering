
#ifndef CPP_IMPL_RKSKDIGITPRECOMP_H
#define CPP_IMPL_RKSKDIGITPRECOMP_H

#include "RLWEKeyswitchInterface.h"

class RKSKDigitPrecomp : public RKSKInterface {

public:

    RKSKDigitPrecomp(std::shared_ptr<RingGSWCryptoParams> params, uint64_t L_bits, uint64_t L_LWE_bits, uint64_t n_LWE, uint32_t max_ct);

    RKSKDigitPrecomp(std::shared_ptr<RingGSWCryptoParams> params, uint64_t L_bits, uint32_t max_ct);

    void Keygen(NativeVector &sk_in, lbcrypto::NativePoly &sk_out) override;

    RLWECiphertext DoKeyswitch(const std::vector<LWECiphertext> &ct_in) override;

    ~RKSKDigitPrecomp() = default;

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

    std::vector<std::vector<std::vector<RLWECiphertext>>> m_ksk;

    std::vector<std::vector<std::vector<LWECiphertext>>> m_LWE_ksk;

    bool m_use_intermediate = false;

    std::vector<NativePoly> m_offset_polys;

};


#endif //CPP_IMPL_RKSKDIGITPRECOMP_H
