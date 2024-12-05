
#ifndef CPP_IMPL_BITINCREASEPRECISION_H
#define CPP_IMPL_BITINCREASEPRECISION_H

#include "binfhecontext.h"

#include "LUTEvalParams.h"

using namespace lbcrypto;

namespace LUTEval {

    LWECiphertext BitIncreasePrecision(std::vector<LWECiphertext>& lwe_in, uint32_t target_space, LUTEvalParams& keys, LWECiphertext msb = nullptr, LWE_STATE insert_before = PRE_FIRST_MODSWITCH);

}

#endif //CPP_IMPL_BITINCREASEPRECISION_H
