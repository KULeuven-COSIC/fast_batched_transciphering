
#include <cassert>
#include <chrono>
#include "LUTEvaluator.h"
#include "BitIncreasePrecision.h"

namespace LUTEval {

    std::vector<LWECiphertext>
    LUTEvaluator::EvalLUT(std::vector<LWECiphertext> &cts, LUTEval::LWE_STATE input_format,
                          LUTEval::LWE_STATE output_format) {
        NativeInteger scale = (2u << (m_io_bits - m_LUT_params.int_param_pt_bits));

        auto cts_fixed = Finalize(cts, input_format, FINALIZED);

        auto mon = ComputeMonomial(cts, scale, input_format);

        auto outputs_QN = EvaluateDecisionTree(cts_fixed, mon);

        auto v = Finalize(outputs_QN, PRE_FIRST_MODSWITCH, output_format);


        //auto cts_fixed = Finalize(cts, input_format, FINALIZED);
        //auto mon = ComputeMonomial(cts_fixed, scale);
        //auto outputs_QN = EvaluateDecisionTree(cts_fixed, mon);
        //auto v = Finalize(outputs_QN, PRE_FIRST_MODSWITCH, output_format);


        return v;
    }

    RLWECiphertext
    LUTEvaluator::ComputeMonomial(std::vector<LWECiphertext> cts, NativeInteger scale, LWE_STATE input_state) {

        auto MonRGSWParams = m_LUT_params.GetMonRGSWParams();
        auto MonPolyParams = MonRGSWParams->GetPolyParams();
        auto MonLWEScheme = m_LUT_params.GetLWEScheme();
        auto MonLWEParams = m_LUT_params.GetLWEParams();

        /* compose first bits into integer */

        std::vector<LWECiphertext> to_compose(cts.begin(), cts.begin() + m_LUT_params.prec_param_pt_bits - 1);
        auto to_compose_final = Finalize(to_compose, input_state, FINALIZED);

        uint32_t target_space = 1u << m_LUT_params.prec_param_pt_bits;
        auto composed_bits = BitIncreasePrecision(to_compose_final, target_space, m_LUT_params, cts[m_LUT_params.prec_param_pt_bits - 1], input_state);

        auto N = MonRGSWParams->GetN();

        /* compute the (padded) monomial we want */
        NativePoly zero = NativePoly(MonPolyParams, EVALUATION, true);
        NativePoly monBN = NativePoly(MonPolyParams, COEFFICIENT, true);


        /* compute monomial for B to pre-rotate, but modswitch the value to N not 2 N */
        monBN[composed_bits->GetB().ConvertToInt() * (N / MonRGSWParams->Getq().ConvertToInt())] = 1;
        monBN.SetFormat(EVALUATION);
        monBN *= monP;

        std::vector<NativePoly> accP = {zero, monBN};
        auto acc = std::make_shared<RLWECiphertextImpl>(accP);

        NativeVector A = composed_bits->GetA();
        A.SetModulus(m_LUT_params.GetLWEParams()->Getq() * 2);

        m_LUT_params.GetAccScheme().EvalAcc(m_LUT_params.GetMonRGSWParams(), m_LUT_params.GetMonAccKey(), acc, A);


        return acc;
    }

    std::vector<LWECiphertext>
    LUTEvaluator::EvaluateDecisionTree(std::vector<LWECiphertext> &cts, RLWECiphertext &mono) {

        auto LWEParams = m_LUT_params.GetLWEParams();
        auto MuxRGSWParams = m_LUT_params.GetMuxRGSWParams();
        auto MuxPolyParams = MuxRGSWParams->GetPolyParams();
        auto MuxQ = MuxRGSWParams->GetQ();

        auto Muxq = MuxRGSWParams->Getq();

        auto MuxN = MuxRGSWParams->GetN();
        auto Bits = m_LUT_params.int_param_pt_bits;


        /* prepare decision tree evaluation */
        std::vector<std::vector<LWECiphertext>> possible_outputs = ComputePossibleOutputs(mono);

        /* Find how much needs to be extracted each time */
        uint32_t n_ct_left = (1 << (m_io_bits - Bits)) * m_io_bits;
        auto capability = m_LUT_params.GetCapacity();
        std::vector<uint32_t> extractionCount;

        n_ct_left >>= 1;
        while (n_ct_left > capability) {
            extractionCount.push_back(capability);
            n_ct_left >>= 1;
        }
        while (n_ct_left >= m_io_bits) {
            extractionCount.push_back(n_ct_left);
            n_ct_left >>= 1;
        }

        /******* Evaluate the tree ********/
        // We will shift our message space by q/4, so this needs to be taken into account and hence we will
        // multiply by N2Poly

        uint32_t bit_idx = Bits;
        auto factor = (2 * MuxN) / Muxq.ConvertToInt();


        while (bit_idx < m_io_bits) {
            LWECiphertext current_bit;
            // std::cout << "EXT " << extractionCount[bit_idx - IntBits] << std::endl;
            if (m_LUT_params.m_boot_before_mux) {
                /* In this branch we simply refresh the bit before muxing */
                NativeInteger scale_factor = MuxQ.DivideAndRound(4);

                NativePoly accumulator_poly = NativePoly(MuxPolyParams,COEFFICIENT,true);
                for(uint32_t i = 0; i < (MuxN >> 1); i++)
                    accumulator_poly[i] = MuxQ.ModSub(scale_factor, MuxQ);
                for(uint32_t i = (MuxN >> 1); i < MuxN; i++)
                    accumulator_poly[i] = scale_factor;

                accumulator_poly.SetFormat(EVALUATION);

                auto aN_ref = cts[bit_idx]->GetA();
                auto bN_ref = cts[bit_idx]->GetB() ;
                /* shift distribution */
                bN_ref.MulEq(2);

                auto b_monomial = NativePoly(MuxPolyParams, COEFFICIENT, true);
                auto b_index = bN_ref.ConvertToInt() % MuxN ;

                b_monomial[b_index] = 1;
                b_monomial.SetFormat(EVALUATION);
                if (bN_ref >= MuxN)
                    b_monomial = -b_monomial;

                std::vector<NativePoly> RLWEcontainer = {NativePoly(MuxPolyParams, EVALUATION, true), accumulator_poly * b_monomial};
                RLWECiphertext acc_bit = std::make_shared<RLWECiphertextImpl>(RLWEcontainer);

                aN_ref.SetModulus(2 * MuxN);
                aN_ref *= (2 * LWEParams->GetN()) / LWEParams->Getq();

                m_LUT_params.accumulator_scheme.EvalAcc(MuxRGSWParams, m_LUT_params.GetCompAccKey(), acc_bit, aN_ref);

                auto a_p = acc_bit->GetElements()[0];
                a_p.Transpose();
                a_p.SetFormat(COEFFICIENT);

                auto b_p = acc_bit->GetElements()[1];
                b_p.SetFormat(COEFFICIENT);

                NativeInteger bQ = b_p[0] + MuxQ.DividedBy(4);
                auto LWEScheme = m_LUT_params.comp_context.GetLWEScheme();
                auto ct_out_N = std::make_shared<LWECiphertextImpl>(a_p.GetValues(), bQ);
                auto ct_out_MS = LWEScheme->ModSwitch(LWEParams->GetqKS(), ct_out_N);
                auto ct_out_KS = LWEScheme->KeySwitch(LWEParams, m_LUT_params.comp_context.GetSwitchKey(), ct_out_MS);
                current_bit = LWEScheme->ModSwitch(LWEParams->Getq(), ct_out_KS);
            } else {
                current_bit = cts[bit_idx];
            }

            auto b2N = current_bit->GetB().ModAdd(Muxq.DividedBy(4), Muxq) * factor;

            NativePoly bP = NativePoly(MuxPolyParams, COEFFICIENT, true);
            bP[b2N.ConvertToInt() % MuxN] = b2N >= MuxN ? MuxQ.ModSub(1,MuxQ) : 1;
            bP.SetFormat(EVALUATION);

            bP *= N2Poly;

            auto packed_for_mux = PackForMux(possible_outputs);

            std::vector<std::vector<LWECiphertext>> new_outputs;

            for(auto& pair : packed_for_mux) {

                auto pair_out = PackedMux(pair.first, pair.second, bP, current_bit->GetA(), extractionCount[bit_idx - Bits]);
                for(uint32_t i = 0; i < pair_out.size(); i+=m_io_bits) {
                    std::vector<LWECiphertext> entry(pair_out.begin() + i, pair_out.begin() + (i + m_io_bits));
                    new_outputs.push_back(entry);
                }

            }

            possible_outputs.clear();
            possible_outputs = std::move(new_outputs);
            bit_idx++;
        }

        return possible_outputs[0];
    }

    std::vector<std::vector<LWECiphertext>> LUTEvaluator::ComputePossibleOutputs(lbcrypto::RLWECiphertext &mon) {

        auto monA = mon->GetElements()[0];
        auto monB = mon->GetElements()[1];


        std::vector<std::vector<LWECiphertext>> possible_outputs;
        for (auto& val_poly : m_LUTPoly) {
            std::vector<LWECiphertext> suffix_vals;
            for(auto& output_poly : val_poly) {
                //output_poly.SetFormat(EVALUATION);
                auto rotA = monA * output_poly;
                auto rotB = monB * output_poly;


                //auto dec = rotB - rotA * m_LUT_params.skNTT;
                //dec.SetFormat(COEFFICIENT);
                //std::cout << dec[0] << std::endl;


                rotA = rotA.Transpose();
                rotA.SetFormat(COEFFICIENT);
                rotB.SetFormat(COEFFICIENT);

                /*
                NativeInteger dec = rotB[0];
                for(uint32_t i = 0; i < N; i++) {
                    dec.ModSubEq(skN[i].ModMul(rotA[i], Q),Q);
                }
                std::cout << dec << " " << std::endl; */

                suffix_vals.push_back(std::make_shared<LWECiphertextImpl>(rotA.GetValues(), rotB[0]));
            }
            //std::cout << std::endl;
            possible_outputs.push_back(suffix_vals);
        }

        return possible_outputs;
    }


    std::vector<LWECiphertext> LUTEvaluator::PackedMux(std::vector<LWECiphertext> &ct0, lbcrypto::RLWECiphertext &ct1,
                                                       lbcrypto::NativePoly bNShifted, NativeVector aN,
                                                       uint32_t extraction_count) {
        auto params = m_LUT_params.GetMuxRGSWParams();


        //auto& pa = ct0->GetElements()[0];
        //auto& pb = ct0->GetElements()[1];

        auto& qa = ct1->GetElements()[0];
        auto& qb = ct1->GetElements()[1];


        /*
        auto IntQ = params->GetQ();
        auto skN = m_LUT_params.GetRLWEKey();
        auto dec0 = pb - pa * skN;
        auto dec1 = qb - qa * skN; */

        std::vector<NativePoly> MuxBuffer = {qa * bNShifted, qb * bNShifted};

        aN.SetModulus(2 * params->GetN());
        aN *= (2 * params->GetN()) / params->Getq();

        auto mux_acc = std::make_shared<RLWECiphertextImpl>(MuxBuffer);
        m_LUT_params.GetAccScheme().EvalAcc(params, m_LUT_params.GetMuxAccKey(), mux_acc, aN);

        auto extracted = MultiExtract(mux_acc, extraction_count);
        assert(extracted.size() == ct0.size());

        std::vector<LWECiphertext> output_data;
        for(uint32_t i = 0; i < extracted.size(); i++) {

            auto A0 = ct0[i]->GetA();
            auto B0 = ct0[i]->GetB();
            auto A1 = extracted[i]->GetA();
            auto B1 = extracted[i]->GetB();

            output_data.push_back(std::make_shared<LWECiphertextImpl>(A0+A1,B0.ModAdd(B1, params->GetQ())));
        }


        return output_data;
        //return MultiExtract(mux_acc, extraction_count);
    }


    void LUTEvaluator::GenerateLUTPolys(LUTFunc &F) {

        auto MonRGSWParams = m_LUT_params.GetMonRGSWParams();
        auto MonPolyParams = MonRGSWParams->GetPolyParams();

        auto prec_param_bits = m_LUT_params.prec_param_pt_bits;

        auto N = MonPolyParams->GetRingDimension();
        auto IntFactor = N >> prec_param_bits;

        for(uint32_t suffix = 0; suffix < (1u << (m_io_bits - prec_param_bits)); suffix++) {
            std::vector<NativePoly> suffix_polys(m_io_bits, NativePoly(MonPolyParams, COEFFICIENT, true));

            for(uint32_t prefix = 0; prefix < (1u << prec_param_bits); prefix++) {
                uint32_t x_in = prefix + (suffix << (prec_param_bits));
                uint32_t f_val = F(x_in);

                for(uint32_t i = 0; i < m_io_bits; i++) {
                    suffix_polys[i][(N - prefix * IntFactor) % N] = (f_val >> i) & 1;
                }

            }
            for(auto& p : suffix_polys)
                p.SetFormat(EVALUATION);

            m_LUTPoly.push_back(suffix_polys);
        }
    }

    std::vector<LWECiphertext> LUTEvaluator::MultiExtract(RLWECiphertext &ct_in, uint32_t amount) {


        auto params = m_LUT_params.GetMuxRGSWParams();
        auto PolyParams = params->GetPolyParams();


        auto aP = ct_in->GetElements()[0];
        auto bP = ct_in->GetElements()[1];

        std::vector<LWECiphertext> ct_out;
        auto max_capa = m_LUT_params.GetCapacity();
        auto rep = max_capa / amount;
        auto EPoly = ExtractorPoly;
        // We multiply a lot of polynomials and do a lot of INTTs here.
        // Technically, we only need 2. Right now we do 2 * amount INNTS (* 16 for all) that are not needed
        while (--rep)
            EPoly *= ExtractorPoly;

        for(uint32_t i = 0; i < amount; i++) {
            NativePoly aT = aP.Transpose();
            NativePoly bC(bP);
            aT.SetFormat(COEFFICIENT);
            bC.SetFormat(COEFFICIENT);

            if (i > 0) {
                aT = -aT;
                bC = -bC;
            }

            ct_out.push_back(std::make_shared<LWECiphertextImpl>(aT.GetValues(), bC[0]));
            aP *= EPoly;
            bP *= EPoly;
        }

        std::reverse(ct_out.begin() + 1, ct_out.end());
        return ct_out;
    }


    std::vector<std::pair<std::vector<LWECiphertext>, RLWECiphertext>>
    LUTEvaluator::PackForMux(std::vector<std::vector<LWECiphertext>> &cts) {

        auto MuxParams = m_LUT_params.GetMuxRGSWParams();
        // For different parameters, an initial modulus switch might be necessary

        auto capability = m_LUT_params.GetCapacity();
        auto n_groups = cts.size();
        auto n_groups2 = n_groups >> 1;
        auto n_out_tuples = std::max<uint64_t>((n_groups * m_io_bits) / (2 * capability), 1);
        auto n_pack = capability / (m_io_bits);
        if (n_groups <= n_pack)
            n_pack = n_groups2;

        std::vector<std::pair<std::vector<LWECiphertext>,RLWECiphertext>> output_data;
        //std::vector<std::pair<RLWECiphertext,RLWECiphertext>> output_data;

        for(uint32_t i = 0; i < n_out_tuples; i++) {

            // values for case bit = 0
            std::vector<LWECiphertext> pack_0;
            // values for case bit = 1
            std::vector<LWECiphertext> pack_1;

            auto offset_i = 2 * i * n_pack;
            for(uint32_t j = 0; j < n_pack; j++) {
                auto idx0 = offset_i + 2 * j;
                auto idx1 = offset_i + 2 * j + 1;
                pack_0.insert(pack_0.end(), cts[idx0].begin(), cts[idx0].end());
                pack_1.insert(pack_1.end(), cts[idx1].begin(), cts[idx1].end());
            }


            NativePoly halfshift = NativePoly(MuxParams->GetPolyParams(), COEFFICIENT, true);
            auto shiftK = MuxParams->GetN() - MuxParams->GetN() / (2 * pack_0.size());
            halfshift[shiftK] = NativeInteger(0).ModSub(1, MuxParams->GetQ());
            halfshift.SetFormat(EVALUATION);

            //RLWECiphertext p = m_RLWESwitchingKey.DoKeyswitch(pack_0);

            std::vector<LWECiphertext> plus;
            std::vector<LWECiphertext> minus;

            for(uint32_t j = 0; j < pack_0.size(); j++) {

                auto A0 = pack_0[j]->GetA();
                auto B0 = pack_0[j]->GetB();
                auto A1 = pack_1[j]->GetA();
                auto B1 = pack_1[j]->GetB();

                plus.push_back(std::make_shared<LWECiphertextImpl>(A0+A1,B0.ModAdd(B1, MuxParams->GetQ())));
                minus.push_back(std::make_shared<LWECiphertextImpl>(A0-A1,B0.ModSub(B1, MuxParams->GetQ())));

            }

            /*
            for(auto& ctp : pack_0) {
                NativeInteger dec = ctp->GetB();
                for(uint32_t k = 0; k < IntParams->GetN(); k++) {
                    dec.ModSubEq(skN[k].ModMul(ctp->GetA(k), IntParams->GetQ()),IntParams->GetQ());
                }
                std::cout << dec << " " << std::endl;
            }
            std::cout << std::endl;
            */

            //p->GetElements()[1] *= halfshift;
            //p->GetElements()[0] *= halfshift;

            //RLWECiphertext q = m_RLWESwitchingKey.DoKeyswitch(pack_1);
            RLWECiphertext  q = m_RLWESwitchingKey.DoKeyswitch(minus);

            q->GetElements()[0] *= halfshift;
            q->GetElements()[1] *= halfshift;

            output_data.emplace_back(plus, q);

        }

        return output_data;

    }

    std::vector<LWECiphertext>
    LUTEvaluator::Finalize(std::vector<LWECiphertext> &cts, LUTEval::LWE_STATE state_in,
                           LUTEval::LWE_STATE state_out) {

        if (state_in == state_out)
            return cts;

        assert(state_in < state_out);

        auto LWEScheme = m_LUT_params.GetLWEScheme();
        auto LWEParams = m_LUT_params.GetLWEParams();

        std::vector<LWECiphertext> outputs;
        std::vector<std::function<LWECiphertext(LWECiphertext&)>> steps = {
                [&](LWECiphertext& x) {return LWEScheme->ModSwitch(LWEParams->GetqKS(), x); },
                [&](LWECiphertext& x) {return LWEScheme->KeySwitch(LWEParams, m_LUT_params.comp_context.GetSwitchKey(), x); },
                [&](LWECiphertext& x) {return LWEScheme->ModSwitch(LWEParams->Getq(), x);}
        };


        state_out = std::min(state_out, PRE_FINAL_MODSWITCH);
        for(auto& ct : cts) {

            auto temp = std::make_shared<LWECiphertextImpl>(ct->GetA(),ct->GetB());
            for(uint32_t i = state_in; i <= state_out; i++) {
                temp = steps[i](temp);
            }
            outputs.push_back(temp);
        }

        return outputs;
    }


}
