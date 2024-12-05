
#ifndef CPP_IMPL_RLWEKEYSWITCHING_H
#define CPP_IMPL_RLWEKEYSWITCHING_H

#include "RLWESwitchingMethods/RLWEKeyswitchInterface.h"

#include "RLWESwitchingMethods/RKSKNoPrecomp.h"
#include "RLWESwitchingMethods/RKSKNoPrecompFast.h"
#include "RLWESwitchingMethods/RKSKDigitPrecompFast.h"
#include "RLWESwitchingMethods/RKSKDigitPrecomp.h"
#include "RLWESwitchingMethods/RKSKDecimation.h"

enum RLWEKeyswitchingMethod {
    NO_DIGIT_PRECOMP,
    NO_DIGIT_PRECOMP_FAST,
    SINGLE_DIGIT_PRECOMP,
    SINGLE_DIGIT_PRECOMP_FAST,
    DECIMATION_TRANSFORM
};

class RLWEKeyswitchingKey {

public:

    RLWEKeyswitchingKey() = default;

    RLWEKeyswitchingKey(const std::shared_ptr<RingGSWCryptoParams>& params, uint64_t L_bits, uint32_t max_ct, RLWEKeyswitchingMethod method) : m_method(method) {
        switch (m_method) {
            case NO_DIGIT_PRECOMP:
                m_ksk_impl = std::make_unique<RKSKNoPrecomp>(params, L_bits, max_ct);
                break;
            case NO_DIGIT_PRECOMP_FAST:
                m_ksk_impl = std::make_unique<RKSKNoPrecompFast>(params, L_bits, max_ct);
                break;
            case SINGLE_DIGIT_PRECOMP:
                m_ksk_impl = std::make_unique<RKSKDigitPrecomp>(params, L_bits, max_ct);
                break;
            case SINGLE_DIGIT_PRECOMP_FAST:
                m_ksk_impl = std::make_unique<RKSKDigitPrecompFast>(params, L_bits, max_ct);
                break;
            case DECIMATION_TRANSFORM:
                m_ksk_impl = std::make_unique<RKSKDecimation>(params, L_bits, max_ct);
                break;
        }
    }

    RLWEKeyswitchingKey(const std::shared_ptr<RingGSWCryptoParams> &params, uint64_t L_bits,
                                       uint64_t L_LWE_bits, uint64_t n_LWE, uint32_t max_ct, RLWEKeyswitchingMethod method) : m_method(method) {
        switch (m_method) {
            case NO_DIGIT_PRECOMP:
                m_ksk_impl = std::make_unique<RKSKNoPrecomp>(params, L_bits, L_LWE_bits, n_LWE, max_ct);
                break;
            case NO_DIGIT_PRECOMP_FAST:
                m_ksk_impl = std::make_unique<RKSKNoPrecompFast>(params, L_bits, L_LWE_bits, n_LWE, max_ct);
                break;
            case SINGLE_DIGIT_PRECOMP:
                m_ksk_impl = std::make_unique<RKSKDigitPrecomp>(params, L_bits, L_LWE_bits, n_LWE, max_ct);
                break;
            default:
                case SINGLE_DIGIT_PRECOMP_FAST:
                m_ksk_impl = std::make_unique<RKSKDigitPrecompFast>(params, L_bits, L_LWE_bits, n_LWE, max_ct);
                break;

        }
    }

    void Keygen(NativeVector sk_in, NativePoly sk_out) {
        m_ksk_impl->Keygen(sk_in, sk_out);
        m_keygen_called = true;
    }

    RLWECiphertext DoKeyswitch(const std::vector<LWECiphertext> & cts) {
        if (!m_keygen_called) {
            std::cerr << "Keygen needs to be called before performing the key switch !" << std::endl;
            std::exit(-1);
        }
        return m_ksk_impl->DoKeyswitch(cts);
    }

private:

    RLWEKeyswitchingMethod m_method;

    std::unique_ptr<RKSKInterface> m_ksk_impl;

    bool m_keygen_called = false;

};

#endif //CPP_IMPL_RLWEKEYSWITCHING_H
