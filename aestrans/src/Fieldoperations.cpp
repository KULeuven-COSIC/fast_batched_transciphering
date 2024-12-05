
#include "FieldOperations.h"

namespace AESTransciphering {

    ByteVec XOR_VEC(ByteVec& lhs, ByteVec& rhs) {
        ByteVec out;
        for(uint32_t i = 0; i < lhs.size(); i++) {
            out[i] = XOR(lhs[i], rhs[i]);
        }
        return out;
    }

    void XOR_VEC_EQ(ByteVec& lhs, ByteVec& rhs) {
        for(uint32_t i = 0; i < lhs.size(); i++) {
            XOREQ(lhs[i], rhs[i]);
        }
    }


    lbcrypto::LWECiphertext XOR(lbcrypto::LWECiphertext& a, lbcrypto::LWECiphertext& b) {
        NativeVector new_a = a->GetA().ModAdd(b->GetA());
        NativeInteger new_b = a->GetB().ModAdd(b->GetB(), a->GetA().GetModulus());
        return std::make_shared<lbcrypto::LWECiphertextImpl>(new_a, new_b);
    }

    void XOREQ(lbcrypto::LWECiphertext& a, lbcrypto::LWECiphertext& b) {
        auto& A = a->GetA();
        auto& B = a->GetB();

        A.ModAddEq(b->GetA());
        B.ModAddEq(b->GetB(), b->GetModulus());

    }



    ByteVec GMUL(uint8_t scal, AESTransciphering::ByteVec &bits)
    {
        switch (scal) {
            case 1: return bits;
            case 2:
                return GMUL_TMP<2>(bits);
            case 3:
                return GMUL_TMP<3>(bits);
            case 9:
                return GMUL_TMP<9>(bits);
            case 11:
                return GMUL_TMP<11>(bits);
            case 13:
                return GMUL_TMP<13>(bits);
            case 14:
                return GMUL_TMP<14>(bits);
            default:
                assert(false);
                return ByteVec();
        }
    }
}

