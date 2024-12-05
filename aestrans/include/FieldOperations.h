
#ifndef CPP_IMPL_FIELDOPERATIONS_H
#define CPP_IMPL_FIELDOPERATIONS_H

#include "binfhecontext.h"

#include <bitset>
#include <cassert>

namespace AESTransciphering {

    using AESBlock128 = std::bitset<128>;
    using ByteVec = std::array<lbcrypto::LWECiphertext, 8>;
    using BitTensor = std::array<std::array<ByteVec, 4>, 4>;

    lbcrypto::LWECiphertext XOR(lbcrypto::LWECiphertext& a, lbcrypto::LWECiphertext& b);

    void XOREQ(lbcrypto::LWECiphertext& a, lbcrypto::LWECiphertext& b);

    ByteVec XOR_VEC(ByteVec& lhs, ByteVec& rhs);

    void XOR_VEC_EQ(ByteVec& lhs, ByteVec& rhs);

    template<uint8_t S> ByteVec GMUL_TMP(ByteVec& bits) {
        ByteVec acc;

        lbcrypto::LWECiphertext zero = std::make_shared<lbcrypto::LWECiphertextImpl>(bits[0]->GetA() * 0, bits[0]->GetB() * 0);

        for(uint32_t i = 0; i < 8; i++)
            acc[i] = zero;

        for(int i = 0; i < 8; i++) {
            if ((S >> i) & 1)  {
                ByteVec tmp_i;
                // shift
                for(int j = 0; j < 8 - i; j++)
                    tmp_i[i+j] = bits[j];
                for(int j = 0; j < i; j++)
                    tmp_i[j] = zero;
                // take care of carry
                for (int j = 7; j > 7 - i; j--) {
                    XOREQ(tmp_i[(j+i)%8],bits[j]);
                    XOREQ(tmp_i[(j+i+1)%8],bits[j]);
                    XOREQ(tmp_i[(j+i+3)%8],bits[j]);
                    XOREQ(tmp_i[(j+i+4)%8],bits[j]);
                    //tmp_i[(j+i) % 8] = XOR(tmp_i[(j+i)%8],bits[j]);
                    //tmp_i[(j+i+1) % 8] = XOR(tmp_i[(j+i+1)%8],bits[j]);
                    //tmp_i[(j+i+3) % 8] = XOR(tmp_i[(j+i+3)%8],bits[j]);
                    //tmp_i[(j+i+4) % 8] = XOR(tmp_i[(j+i+4)%8],bits[j]);
                }
                XOR_VEC_EQ(acc, tmp_i);
            }
        }

        return acc;
    }

    ByteVec GMUL(uint8_t scal, ByteVec& bits);

}

#endif //CPP_IMPL_FIELDOPERATIONS_H
