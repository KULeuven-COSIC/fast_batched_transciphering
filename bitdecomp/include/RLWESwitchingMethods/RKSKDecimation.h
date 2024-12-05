
#ifndef CPP_IMPL_RKSKDECIMATION_H
#define CPP_IMPL_RKSKDECIMATION_H

#include "RLWEKeyswitchInterface.h"

class RKSKDecimation : public RKSKInterface {

public:

    RKSKDecimation(std::shared_ptr<RingGSWCryptoParams> params, uint64_t L_bits, uint32_t max_ct);

    void Keygen(NativeVector &sk_in, lbcrypto::NativePoly &sk_out) override;

    RLWECiphertext DoKeyswitch(const std::vector<LWECiphertext> &ct_in) override;

    ~RKSKDecimation() = default;

private:

    RLWECiphertext DoRLWEToRLWESwitch(std::vector<RLWECiphertext> &ksk, NativePoly &A_in, NativePoly &B_in);

    std::vector<std::vector<RLWECiphertext>> m_auto_keys;
    std::shared_ptr<RingGSWCryptoParams> m_params;
    std::vector<NativePoly> m_offset_polys;

    NativePoly sk;

    uint32_t m_max_ct;

    uint64_t m_L;
    uint64_t m_L_bits;
    uint64_t m_digits;

};

#endif //CPP_IMPL_RKSKDECIMATION_H
