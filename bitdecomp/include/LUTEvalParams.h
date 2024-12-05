
#ifndef CPP_IMPL_LUTEVALPARAMS_H
#define CPP_IMPL_LUTEVALPARAMS_H

#include "binfhecontext.h"
#include "RLWEKeyswitching.h"

using namespace lbcrypto;

#define USE_CGGI

#ifdef USE_CGGI
    #include "rgsw-acc-cggi.h"
    using AccumulatorScheme = RingGSWAccumulatorCGGI;
    const auto method = GINX;
#elif defined(USE_FHEW)
    #include "rgsw-acc-dm.h"
    const auto method = AP;
    using AccumulatorScheme = RingGSWAccumulatorDM;
#else
    // can't use that one since OPENFHE only applies the automorphism transform to the rhs of the accumulator...
    #include "rgsw-acc-lmkcdey.h"
    const auto method = LMKCDEY;
    using AccumulatorScheme = RingGSWAccumulatorLMKCDEY;
#endif


namespace LUTEval {

    enum LWE_STATE {
        PRE_FIRST_MODSWITCH = 0,
        PRE_KEYSWITCH = 1,
        PRE_FINAL_MODSWITCH = 2,
        FINALIZED = 3
    };

    struct LUTEvalParams {

        AccumulatorScheme accumulator_scheme;

        BinFHEContext comp_context;

        // use same LWE and RLWE key everywhere
        LWEPrivateKey sk;
        NativePoly skN, skNTT;
        bool m_boot_before_mux = false;

        std::shared_ptr<RingGSWCryptoParams> m_comp_rgsw_params;
        RingGSWACCKey m_comp_boot_key;

        std::shared_ptr<RingGSWCryptoParams> m_mon_rgsw_params;
        RingGSWACCKey m_mon_boot_key;

        std::shared_ptr<RingGSWCryptoParams> m_mux_rgsw_params;
        RingGSWACCKey m_mux_boot_key;

        /* for LUT eval */
        uint32_t int_param_pt_bits = 5;
        uint32_t prec_param_pt_bits = 5;

        NativeInteger m_sk_sum = 0;

        LUTEvalParams() : comp_context() {
            SetupKeys(SET_IV);
        }

        void SetupKeys(BINFHE_PARAMSET set) {

            uint32_t ham = 100;

            comp_context.GenerateBinFHEContext(set, method);
            auto sk_tmp = comp_context.KeyGen();
            NativeVector vec = sk_tmp->GetElement();
            auto ctr = 0;

            // make key sparse and compute sum...
            m_sk_sum = 0;
            for(uint32_t i = 0; i < vec.GetLength(); i++) {
                if (vec[i] != 0 and ctr < ham) {
                    ctr++;
                } else {
                    if (vec[i] != 0 and ctr >= ham) {
                        vec[i] = 0;
                    }
                }
                m_sk_sum.ModAddEq(vec[i], comp_context.GetParams()->GetLWEParams()->Getq());
            }
            //std::cout << "Sum of secret is " << m_sk_sum << std::endl;
            sk = std::make_shared<LWEPrivateKeyImpl>(vec);
            comp_context.BTKeyGen(sk);
            m_comp_rgsw_params = comp_context.GetParams()->GetRingGSWParams();
            m_comp_boot_key = comp_context.GetRefreshKey();

            auto baseG = m_comp_rgsw_params->GetBaseG();
            auto BTKey = comp_context.GetBTKeyMap()->at(baseG);
            skN = BTKey.skN;
            skN.SetFormat(COEFFICIENT);

            skNTT = skN;
            skNTT.SetFormat(EVALUATION);

            auto q = m_comp_rgsw_params->Getq();
            auto N = m_comp_rgsw_params->GetN();
            auto Q = m_comp_rgsw_params->GetQ();
            auto stddev = comp_context.GetParams()->GetRingGSWParams()->GetDgg().GetStd();

            uint32_t L_mon, L_mux;

            switch (set) {
                case lbcrypto::SET_I: {
                    L_mon = 1 << 2;
                    L_mux = 1 << 3;
                    break;
                }
                case lbcrypto::SET_II: {
                    L_mon = 1 << 5;
                    L_mux = 1 << 7;
                    m_boot_before_mux = true;
                    break;
                }
                case lbcrypto::SET_III: {
                    L_mon = 1 << 9;
                    L_mux = 1 << 12;
                    break;
                }
                case lbcrypto::SET_IV: {
                    L_mon = 1 << 17;
                    L_mux = 1 << 17;
                }
                default: {
                    L_mon = m_comp_rgsw_params->GetBaseG();
                    L_mux = m_comp_rgsw_params->GetBaseG();
                }
            }

            m_mon_rgsw_params = std::make_shared<RingGSWCryptoParams>(N,Q,q,L_mon,23, method, stddev);
            m_mux_rgsw_params = std::make_shared<RingGSWCryptoParams>(N,Q,q,L_mux,23, method, stddev);

            m_mon_boot_key = accumulator_scheme.KeyGenAcc(m_mon_rgsw_params, skNTT, sk);
            m_mux_boot_key = accumulator_scheme.KeyGenAcc(m_mux_rgsw_params, skNTT, sk);
        }

        const std::shared_ptr<LWEEncryptionScheme>& GetLWEScheme() {
            return comp_context.GetLWEScheme();
        }


        const std::shared_ptr<LWECryptoParams>& GetLWEParams() {
            return comp_context.GetParams()->GetLWEParams();
        }


        const std::shared_ptr<RingGSWCryptoParams> GetMonRGSWParams() {
            return m_mon_rgsw_params;
        }

        const std::shared_ptr<RingGSWCryptoParams> GetCompRGSWParams() {
            return m_comp_rgsw_params;
        }

        const std::shared_ptr<RingGSWCryptoParams> GetMuxRGSWParams() {
            return m_mux_rgsw_params;
        }

        const AccumulatorScheme GetAccScheme() {
            return accumulator_scheme;
        }

        ConstRingGSWACCKey GetMonAccKey() {
            return m_mon_boot_key;
        }

        ConstRingGSWACCKey GetCompAccKey() {
            return m_comp_boot_key;
        }

        ConstRingGSWACCKey GetMuxAccKey() {
            return m_mux_boot_key;
        }

        NativePoly GetRLWEKey() {
            return skNTT;
        }

        NativeVector GetRLWEKeyVector() {
            return skN.GetValues();
        }

        uint64_t GetCapacity() {
            return 1ull << int_param_pt_bits;
        }




    };

};

#endif //CPP_IMPL_LUTEVALPARAMS_H
