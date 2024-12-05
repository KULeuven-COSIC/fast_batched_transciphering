
#ifndef CPP_IMPL_RLWEKEYSWITCHINTERFACE_H
#define CPP_IMPL_RLWEKEYSWITCHINTERFACE_H

#include "binfhecontext.h"

using namespace lbcrypto;

using UINTNAT = std::conditional<NATIVEINT == 64, uint64_t, uint32_t>::type;
#if NATIVEINT == 32
inline uint32_t clz_u128(uint64_t X) {
    uint32_t bits = 63;
    for(uint32_t i = 0; i < 64; i++) {
        if ((X >> (bits - i)) & 1) {
            return 64+i-1;
        }
    }
    return 128;
}
#endif

class RKSKInterface {

public:
    virtual void Keygen(NativeVector& sk_in, NativePoly& sk_out) = 0;

    virtual RLWECiphertext DoKeyswitch(const std::vector<LWECiphertext>& ct_in) = 0;

    virtual ~RKSKInterface() = default;

};

#endif //CPP_IMPL_RLWEKEYSWITCHINTERFACE_H
