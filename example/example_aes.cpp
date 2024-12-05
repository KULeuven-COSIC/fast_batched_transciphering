
#include <iostream>
#include "AESState.h"
#include "binfhecontext.h"

using namespace AESTransciphering;

int main() {

    // need to manually set n: (I: 512,II:600,III:448)
    uint32_t n = 448;
    uint32_t q = 1<<10;

    AESTransciphering::BitTensor T;

    srand(time(nullptr));
    DiscreteUniformGeneratorImpl<NativeVector> dug;
    dug.SetModulus(q);

    for(uint32_t i = 0; i < 4; i++) {
        for(uint32_t j = 0; j < 4; j++) {
            for(uint32_t k = 0; k < 8; k++) {
                NativeVector v = dug.GenerateVector(n);
                v *= 0;
                uint32_t bit = rand() % 2;
                NativeInteger b = (q >> (bit)) % q;
                T[i][j][k] = std::make_shared<lbcrypto::LWECiphertextImpl>(v,b);
            }
            //std::cout << acc << " ";
        }
    }

    std::vector<AESTransciphering::BitTensor> schedule;
    uint32_t n_rounds = 10;
    for(uint32_t r = 0; r < n_rounds; r++) {
        AESTransciphering::BitTensor key_r;
        for(uint32_t i = 0; i < 4; i++) {
            for(uint32_t j = 0; j < 4; j++) {
                for(uint32_t k = 0; k < 8; k++) {
                    NativeVector v = NativeVector(n, q);
                    v *= 0;
                    NativeInteger b = 0;
                    key_r[i][j][k] = std::make_shared<lbcrypto::LWECiphertextImpl>(v,b);
                }
            }
        }
        schedule.push_back(key_r);
    }

    AESTransciphering::HAES transcipherer(schedule);
    transcipherer.Transcipher(T);

    std::cout << T[0][0][0]->GetB() << std::endl;

}