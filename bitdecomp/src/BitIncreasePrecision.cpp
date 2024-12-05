
#include "BitIncreasePrecision.h"

namespace LUTEval {

    LWECiphertext BitIncreasePrecision(std::vector<LWECiphertext>&lwe_in, uint32_t target_space, LUTEvalParams& keys, LWECiphertext MSB, LWE_STATE insert_before){

        auto LWEscheme = keys.GetLWEScheme();
        auto LWEparams = keys.GetLWEParams();

        auto RGSWParams = keys.GetCompRGSWParams();
        auto PolyParams = RGSWParams->GetPolyParams();

        auto Q = RGSWParams->GetQ();
        auto N = RGSWParams->GetN();

        NativeInteger scale_factor = Q.DivideAndRound(2 * target_space);

        NativePoly accumulator_poly = NativePoly(PolyParams,COEFFICIENT,true);
        for(uint32_t i = 0; i < (N >> 1); i++)
            accumulator_poly[i] = Q.ModSub(scale_factor, Q);
        for(uint32_t i = (N >> 1); i < N; i++)
            accumulator_poly[i] = scale_factor;

        accumulator_poly.SetFormat(EVALUATION);

        NativePoly zero = NativePoly(PolyParams, EVALUATION, true);

        std::vector<NativePoly> RLWEcontainer;

        NativePoly integer_acc_a = NativePoly(PolyParams, EVALUATION, true);
        NativePoly integer_acc_b = NativePoly(PolyParams, EVALUATION, true);
        NativePoly offset_poly = NativePoly(PolyParams, COEFFICIENT, true);
        offset_poly[0] = scale_factor;
        offset_poly.SetFormat(EVALUATION);

        auto accumulator_scheme = keys.GetAccScheme();

        for(auto& bit : lwe_in) {

            //auto ct_N = LWEscheme->ModSwitch(2 * N, bit);
            auto aN = bit->GetA();
            auto bN = bit->GetB() ;

            /* shift distribution */
            bN.MulEq((2 * LWEparams->GetN()) / LWEparams->Getq());
            //bN.ModSubEq(2 * keys.m_sk_sum, 2 * aN.GetModulus());

            auto b_monomial = NativePoly(PolyParams, COEFFICIENT, true);
            auto b_index = bN.ConvertToInt() % N;
            b_monomial[b_index] = 1;
            b_monomial.SetFormat(EVALUATION);
            if (bN >= N)
                b_monomial = -b_monomial;

            RLWEcontainer = {zero, accumulator_poly * b_monomial};
            RLWECiphertext acc_bit = std::make_shared<RLWECiphertextImpl>(RLWEcontainer);

            aN.SetModulus(2 * LWEparams->GetN());
            aN *= (2 * LWEparams->GetN()) / LWEparams->Getq();

            accumulator_scheme.EvalAcc(RGSWParams, keys.GetCompAccKey(), acc_bit, aN);

            integer_acc_a += acc_bit->GetElements()[0];
            integer_acc_b += acc_bit->GetElements()[1] + offset_poly;

            accumulator_poly *= 2;
            offset_poly *= 2;
        }

        integer_acc_b.SetFormat(COEFFICIENT);


        auto aT = integer_acc_a.Transpose();
        aT.SetFormat(COEFFICIENT);

        if (MSB == nullptr) {
            auto ct_out_N = std::make_shared<LWECiphertextImpl>(aT.GetValues(), integer_acc_b[0]);
            auto ct_out_MS = LWEscheme->ModSwitch(LWEparams->GetqKS(), ct_out_N);
            auto ct_out_KS = LWEscheme->KeySwitch(LWEparams, keys.comp_context.GetSwitchKey(), ct_out_MS);
            auto ct_out = LWEscheme->ModSwitch(LWEparams->Getq(), ct_out_KS);

            return ct_out;
        } else {

            switch (insert_before) {
                case PRE_FIRST_MODSWITCH: {
                    auto ct_out_N = std::make_shared<LWECiphertextImpl>(aT.GetValues() + MSB->GetA(),
                                                                        integer_acc_b[0].ModAdd(MSB->GetB(), Q));
                    auto ct_out_MS = LWEscheme->ModSwitch(LWEparams->GetqKS(), ct_out_N);
                    auto ct_out_KS = LWEscheme->KeySwitch(LWEparams, keys.comp_context.GetSwitchKey(), ct_out_MS);
                    auto ct_out = LWEscheme->ModSwitch(LWEparams->Getq(), ct_out_KS);
                    return ct_out;
                }
                case PRE_KEYSWITCH: {
                    auto Qks = LWEparams->GetqKS();
                    auto ct_out_N = std::make_shared<LWECiphertextImpl>(aT.GetValues(), integer_acc_b[0]);
                    auto ct_out_MS = LWEscheme->ModSwitch(LWEparams->GetqKS(), ct_out_N);
                    ct_out_MS->GetA().ModAddEq(MSB->GetA());
                    ct_out_MS->GetB().ModAddEq(MSB->GetB(), Qks);
                    auto ct_out_KS = LWEscheme->KeySwitch(LWEparams, keys.comp_context.GetSwitchKey(), ct_out_MS);
                    auto ct_out = LWEscheme->ModSwitch(LWEparams->Getq(), ct_out_KS);
                    return ct_out;
                }
                case PRE_FINAL_MODSWITCH: {
                    auto Qks = LWEparams->GetqKS();
                    auto ct_out_N = std::make_shared<LWECiphertextImpl>(aT.GetValues(), integer_acc_b[0]);
                    auto ct_out_MS = LWEscheme->ModSwitch(LWEparams->GetqKS(), ct_out_N);
                    auto ct_out_KS = LWEscheme->KeySwitch(LWEparams, keys.comp_context.GetSwitchKey(), ct_out_MS);
                    ct_out_KS->GetA().ModAddEq(MSB->GetA());
                    ct_out_KS->GetB().ModAddEq(MSB->GetB(), Qks);
                    auto ct_out = LWEscheme->ModSwitch(LWEparams->Getq(), ct_out_KS);
                    return ct_out;
                }
                default: {
                    auto ct_out_N = std::make_shared<LWECiphertextImpl>(aT.GetValues(), integer_acc_b[0]);
                    auto ct_out_MS = LWEscheme->ModSwitch(LWEparams->GetqKS(), ct_out_N);
                    auto ct_out_KS = LWEscheme->KeySwitch(LWEparams, keys.comp_context.GetSwitchKey(), ct_out_MS);
                    auto ct_out = LWEscheme->ModSwitch(LWEparams->Getq(), ct_out_KS);
                    ct_out->GetA().ModAddEq(MSB->GetA());
                    ct_out->GetB().ModAddEq(MSB->GetB(), ct_out->GetModulus());
                    return ct_out;
                }
            }
        }
    }
}
