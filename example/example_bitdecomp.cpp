

#include <iostream>

#include "binfhecontext.h"

#include "LUTEvalParams.h"


using namespace lbcrypto;
using namespace LUTEval;


#include <chrono>
#include <cassert>
#include "RLWEKeyswitching.h"
#include "LUTEvaluator.h"

using namespace std::chrono;

NativeInteger quick_decrypt(LWECiphertext& ct, const NativeVector& s, NativeInteger modulus) {

    NativeInteger inner = ct->GetB();
    for(uint32_t i = 0; i < ct->GetA().GetLength(); i++) {
        auto prod = ct->GetA(i).ModMulFast(s[i], modulus);
        inner.ModSubEq(prod, modulus);
    }
    return inner;
}

int main() {

    // setup context
    LUTEvalParams params;

    uint32_t pt_bits = 8;
    uint32_t pt = 1 << pt_bits;

    auto ctx = params.comp_context;

    /* Choose random message in plaintext space */
    srand(time(nullptr));
    LWEPlaintext msg = rand() % pt;
    std::cout << "Input message was : " << msg << std::endl;

    /* Encrypt bits of the message */
    std::vector<LWECiphertext> bitvector;
    DiscreteUniformGeneratorImpl<NativeVector> dug;

    auto q = params.GetLWEParams()->Getq();

    dug.SetModulus(q);
    auto skv = params.sk->GetElement();
    skv.SwitchModulus(q);
    for(uint32_t i = 0; i < pt_bits; i++) {
        auto msgbit = (msg >> i) & 1;
        //auto ct = ctx.Encrypt(params.sk, msgbit, FRESH, 2);
        NativeVector a = dug.GenerateVector(params.GetLWEParams()->Getn());
        NativeInteger b = msgbit * q.DividedBy(2);
        for(uint32_t j = 0; j < params.sk->GetLength(); j++) {
            b.ModAddEq(a[j].ModMul(skv[j], q), q);
        }
        auto ct = std::make_shared<LWECiphertextImpl>(a,b);
        bitvector.push_back(ct);
    }

    std::function<uint32_t(uint32_t)> F = [](uint32_t x) {return (5 * x) ^ 41;};

    //// create LUT evaluator
    // setup rlwe ksk
    // If we use SET_IV set @method to DECIMATION_TRANSFORM
    auto RLWE_L_bits = 7;
    RLWEKeyswitchingKey RLWEKey = RLWEKeyswitchingKey(params.GetMonRGSWParams(), RLWE_L_bits, params.GetCapacity(), NO_DIGIT_PRECOMP);
    RLWEKey.Keygen(params.GetRLWEKeyVector(), params.GetRLWEKey());

    auto eval = LUTEval::LUTEvaluator(params, RLWEKey,F,pt_bits);

    auto start_time = high_resolution_clock::now();
    auto lut_out = eval.EvalLUT(bitvector, LUTEval::FINALIZED, LUTEval::FINALIZED);

    auto stop_time = high_resolution_clock::now();

    // Decrypt and print stats /
    std::vector<NativeInteger> decrypted;
    for(auto& ct : lut_out) {
        auto phase = quick_decrypt(ct, params.sk->GetElement(), params.GetLWEParams()->Getq());
        std::cout << "Phase is " << phase << std::endl;
        phase.MultiplyAndRoundEq(2, params.GetLWEParams()->Getq());
        phase.ModEq(2);
        decrypted.push_back(phase);
    }

    std::cout << "After reassembling the bits : " << std::endl;
    for(auto& v : decrypted)
        std::cout << v << " ";
    std::cout << std::endl;

    auto elapsed_time = duration_cast<milliseconds>(stop_time - start_time).count();
    std::cout << "Bit composition took " << elapsed_time << "ms." << std::endl;
}
