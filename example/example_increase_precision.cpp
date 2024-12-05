
#include <chrono>

#include "binfhecontext.h"
#include "BitIncreasePrecision.h"

using namespace std::chrono;
using namespace LUTEval;

NativeInteger quick_decrypt(LWECiphertext& ct, LWEPrivateKey s, NativeInteger modulus) {

    NativeInteger inner = ct->GetB();
    for(uint32_t i = 0; i < ct->GetA().GetLength(); i++) {
        auto prod = ct->GetA(i).ModMulFast(s->GetElement()[i], modulus);
        inner.ModSubEq(prod, modulus);
    }
    return inner;
}

int main() {

    LUTEvalParams params;

    uint32_t pt_bits = 4;
    uint32_t pt = 1 << pt_bits;

    auto ctx = params.comp_context;
    NativeInteger q = params.GetLWEParams()->Getq();
    NativeInteger delta = q.DivideAndRound(pt);

    /* Choose random message in plaintext space */
    srand(time(nullptr));
    LWEPlaintext msg = rand() % pt;

    /* Encrypt bits of the message */
    std::vector<LWECiphertext> bitvector;
    for(uint32_t i = 0; i < pt_bits; i++) {
        auto msgbit = (msg >> i) & 1;
        auto ct = ctx.Encrypt(params.sk, msgbit, FRESH, 2);
        bitvector.push_back(ct);
    }

    /* recompose bits into single integer */
    auto start_time = high_resolution_clock::now();
    auto recomposed = BitIncreasePrecision(bitvector, pt, params);
    auto stop_time = high_resolution_clock::now();

    /* Decrypt and print stats */
    auto dec = quick_decrypt(recomposed, params.sk, q);

    auto dec_pt = dec.MultiplyAndRound(pt, q);

    std::cout << "Input message was : " << msg << std::endl;
    std::cout << "After reassembling the bits : " << dec_pt << std::endl;

    auto err = (msg * delta).ModSub(dec, q);
    err = std::min<uint64_t>(q.ModSub(err, q).ConvertToInt(), err.ConvertToInt());
    std::cout << "Error from exact message : " <<  err << std::endl;

    auto elapsed_time = duration_cast<milliseconds>(stop_time - start_time).count();
    std::cout << "Bit composition took " << elapsed_time << "ms." << std::endl;

}