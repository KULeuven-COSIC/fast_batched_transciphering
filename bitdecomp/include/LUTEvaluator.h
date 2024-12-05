
#ifndef CPP_IMPL_LUTEVALUATOR_H
#define CPP_IMPL_LUTEVALUATOR_H

#include "binfhecontext.h"
#include "RLWEKeyswitching.h"
#include "LUTEvalParams.h"


using namespace LUTEval;
using namespace lbcrypto;

namespace LUTEval {

    class LUTEvaluator {

    public:
        typedef std::function<uint32_t(uint32_t)> LUTFunc;

        // Describes the state of a LWE sample in the bootstrapping pipeline

        LUTEvaluator(LUTEvalParams& i_params, RLWEKeyswitchingKey& i_key, LUTFunc& i_F, uint32_t io_bits) : m_LUT_params(i_params), m_RLWESwitchingKey(i_key), m_io_bits(io_bits) {
            GenerateLUTPolys(i_F);

            /**** We precompute polynomials that we keep on using ****/
            auto IntParams = m_LUT_params.GetMonRGSWParams();
            auto IntPolyParams = IntParams->GetPolyParams();
            auto IntQ = IntParams->GetQ();
            auto IntN = IntParams->GetN();
            auto IntFactor2 = IntN >> (m_LUT_params.prec_param_pt_bits + 1);

            monP = NativePoly(IntPolyParams, COEFFICIENT, true);

            NativeInteger scale = (2u << (m_io_bits - m_LUT_params.int_param_pt_bits));

            NativeInteger Qt2 = IntQ.DivideAndRound(scale);
            for(uint32_t i = IntN - IntFactor2; i < IntN + IntFactor2 ; i++) {
                uint32_t idx = i % IntN;
                if (idx == i) {
                    monP[idx] = IntQ.ModSub(Qt2,IntQ);
                } else {
                    monP[idx] = Qt2;
                }
            }

            monP.SetFormat(EVALUATION);

            N2Poly = NativePoly(IntPolyParams, COEFFICIENT, true);
            N2Poly[IntN >> 1] = IntQ.ModSub(1, IntQ);
            N2Poly.SetFormat(EVALUATION);

            ExtractorPoly = NativePoly(IntPolyParams, COEFFICIENT, true);
            ExtractorPoly[IntN >> m_LUT_params.int_param_pt_bits] = 1;
            ExtractorPoly.SetFormat(EVALUATION);

        }

        /**
         * Main method to evaluate a LUT
         * @param cts the ciphertext bits whose composition equals the message the LUT depends upon
         * @param input_format Format for input see LWE_STATE
         * @param output_format Format for output see LWE state
         * @return a vector of bits describing the LUT output
         */
        std::vector<LWECiphertext> EvalLUT(std::vector<LWECiphertext> &cts, LUTEval::LWE_STATE input_format, LUTEval::LWE_STATE output_format);

    private:

        /**
         * Given a vector a ciphertexts containing bits b_0,...,b_k creates a RLWE sample containing a monomial
         * X^{\sum b_i 2^i} (with padding)
         * @param cts the ciphertext containing bits
         * @param scale the scaling factor to use in front of the monomial
         * @return a RLWE sample containing the monomial
         */
        RLWECiphertext ComputeMonomial(std::vector<LWECiphertext> cts, NativeInteger scale, LWE_STATE input_state);

        /**
         * Given subset of ciphertexts containing bits b_{k+1},...,b_{n - 1}, and a RLWE sample containing
         * a monomial generated by \ComputeMonomial where the exponent is induced by b_0,...,b_k
         * computes the correct LUT output by evaluating a decision tree
         * @param cts the ciphertexts containing the bits
         * @param mono the RLWE monomial
         * @return a vector of ciphertexts containing the output bits
         */
        std::vector<LWECiphertext> EvaluateDecisionTree(std::vector<LWECiphertext> &cts, RLWECiphertext &mono);

        /**
         * During the decision tree evaluation we need to chose between to packed sets of outputs. PackedMux returns
         * ct0 if the mux bit was 0, else ct1
         * @param ct0 the first set of packed outputs for the case where bit = 0
         * @param ct1 the second set of packed outputs for the case where bit = 1
         * @param bNShifted b component of the lwe ciphertext encoding the bit
         * @param aN a component of the lwe ciphertext encoding the bit
         * @param extraction_count number of slots to extract from the output
         * @return a vector of lwe ciphertext containing the values from ct0 if bit = 0 else the values from ct1
         */
        std::vector<LWECiphertext> PackedMux(std::vector<LWECiphertext> &ct0, RLWECiphertext &ct1, NativePoly bNShifted, NativeVector aN, uint32_t extraction_count);

        /**
         * Given a RLWE sample encoding a monomial, computes *all* possible outputs for that monomial
         * @param mon the RLWE sample
         * @return each possible set of ciphertext bits conforming to the monomial
         */
        std::vector<std::vector<LWECiphertext>> ComputePossibleOutputs(RLWECiphertext &mon);

        /**
         * Given a function F, generates polynomials F_{s, j} such that (X^m * F_{s, j})_0 = [F(s || m)]_j
         * Specific to this implementation: F_{s,j} will be sparse since we can use a padded monomial
         * @param F a function
         */
        void GenerateLUTPolys(LUTFunc& F);

        /**
         * On input a RLWE sample storing values at full capacity (i.e. using the full polynomial dimension), outputs
         * the \amount values as LWE ciphertexts
         * @param ct_in input RLWE sample
         * @param amount number of values stored in the RLWE sample
         * @return A vector of LWE samples containing the values
         */
        std::vector<LWECiphertext> MultiExtract(RLWECiphertext& ct_in, uint32_t amount);

        /**
         * As part of the LUT evaluation, we need to evaluate a decision tree. To this end,
         * we need to pack our set of LWE samples into 1 or more pairs of RLWE samples to model the cases
         * of a bit to later on select the right one.
         * @param cts vector of LWE ciphertext
         * @return vector of RLWE sample. The first element of a pair is correct for bit = 0, the other is correct for bit = 1
         */
        std::vector<std::pair<std::vector<LWECiphertext>, RLWECiphertext>> PackForMux(
                std::vector<std::vector<LWECiphertext>> &cts);

        /**
         * Helper function. In some cases, we wish to work with lwe samples before a key or modswitch.
         * @param cts the ciphertexts to finalize
         * @param state_in input state of the ciphertexts
         * @param state_out output state of the ciphertexts
         * @return cts in \state_out format
         */
        std::vector<LWECiphertext> Finalize(std::vector<LWECiphertext>& cts, LWE_STATE state_in, LWE_STATE state_out);

        /* polynomials we multiply with */
        std::vector<std::vector<NativePoly>> m_LUTPoly;
        /* LUT Evaluation / Bootstrapping parameters */
        LUTEvalParams& m_LUT_params;
        /* Keyswitching key */
        RLWEKeyswitchingKey& m_RLWESwitchingKey;
        /* Number of input and output bits */
        uint32_t m_io_bits;
        /*** Important Polynomials ***/
        /* Monomial Accumulator */
        NativePoly monP;
        /* Shift poly */
        NativePoly N2Poly;
        /* Extraction */
        NativePoly ExtractorPoly;

    };

};

#endif //CPP_IMPL_LUTEVALUATOR_H