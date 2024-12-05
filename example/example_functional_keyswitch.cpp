
#include "binfhecontext.h"
#include "RLWEKeyswitching.h"

#include <chrono>

using namespace lbcrypto;


int main() {

    auto ref_ctx = BinFHEContext();
    ref_ctx.GenerateBinFHEContext(SET_IV, GINX);
    auto sk_vec = ref_ctx.KeyGenN()->GetElement();

    auto rgswparams = ref_ctx.GetParams()->GetRingGSWParams();
    auto my_params = rgswparams;

    auto N = my_params->GetN();
    auto Q = my_params->GetQ();

    auto skN = NativePoly(my_params->GetPolyParams(), COEFFICIENT);
    skN.SetValues(sk_vec, COEFFICIENT);

    uint32_t capa_bits = 5;
    auto switchKey = RLWEKeyswitchingKey(rgswparams, 6, 1 << capa_bits, NO_DIGIT_PRECOMP);
    switchKey.Keygen(sk_vec, skN);

    DiscreteUniformGeneratorImpl<NativeVector> dug;
    dug.SetModulus(Q);

    srand(time(nullptr));
    std::vector<LWECiphertext> to_switch;
    uint32_t n_ct = 1 << capa_bits;
    for(uint32_t i = 0; i <n_ct; i++) {
        NativeVector a = dug.GenerateVector(N);

        auto val = (i);

        NativeInteger b = val;// * Q.DividedBy(16);
        for(uint32_t j = 0; j < N; j++) {
            b.ModAddEq(sk_vec.at(j).ModMul(a.at(j), Q), Q);
        }
        to_switch.push_back(std::make_shared<LWECiphertextImpl>(a,b));
    }
    std::cout << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    auto ct_out = switchKey.DoKeyswitch(to_switch);
    auto stop = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start).count();

    auto a = ct_out->GetElements()[0];
    a.SetFormat(EVALUATION);
    skN.SetFormat(EVALUATION);

    auto as = a * skN;
    auto dec = ct_out->GetElements()[1] - as;
    dec.SetFormat(COEFFICIENT);

    std::cout << dec << std::endl;
    std::cout << "Took " << elapsed << "ms." << std::endl;
}
