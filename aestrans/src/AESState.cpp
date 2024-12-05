
#include "AESState.h"
#include <chrono>


using namespace lbcrypto;
using namespace std::chrono;

namespace AESTransciphering {

    void HAES::Transcipher(AESTransciphering::BitTensor &Output) {

        if (decrypt) {
            // First step
            HAES::AddRoundKey(Output, n_rounds - 1);
            HAES::ShiftRowsInv(Output);
            HAES::SubBytesInv(Output);

            // First n - 1 rounds

            for(int i = 1; i < n_rounds - 1; i++) {
                auto round_start = high_resolution_clock::now();
                HAES::AddRoundKey(Output, n_rounds - i - 1);
                HAES::MixColumnsInv(Output);
                HAES::ShiftRowsInv(Output);
                HAES::SubBytesInv(Output);
                auto round_stop = high_resolution_clock::now();
                auto elapsed = duration_cast<milliseconds>(round_stop-round_start).count();
                std::cout << "[!] Round" << i << " took " << elapsed << "ms." << std::endl;
            }

            // Finalize
            HAES::AddRoundKey(Output, 0);
        } else {
            // First step
            HAES::AddRoundKey(Output, 0);

            // First n - 1 rounds
            for(int i = 1; i < n_rounds - 1; i++) {
                auto round_start = high_resolution_clock::now();
                HAES::SubBytes(Output);
                HAES::ShiftRows(Output);
                HAES::MixColumns(Output);
                HAES::AddRoundKey(Output, i);
                auto round_stop = high_resolution_clock::now();
                auto elapsed = duration_cast<milliseconds>(round_stop-round_start).count();
                std::cout << "[!] Round" << i << " took " << elapsed << "ms." << std::endl;
            }

            // Finalize
            HAES::SubBytes(Output);
            HAES::ShiftRows(Output);
            HAES::AddRoundKey(Output, n_rounds - 1);
        }


    }

    void HAES::Transcipher(AESTransciphering::BitTensor &Output, AESTransciphering::AESBlock128 Block) {

        // Initialize bittensor
        for(int i = 0; i < 4; i++) {
            for(int j = 0; j < 4; j++) {
                for(int k = 0; k < 8; k++) {
                    auto bit_idx = k + 8 * (j + 4 * i);
                    auto zero_vec = NativeVector(Output[0][0][0]->GetLength(), Output[0][0][0]->GetModulus(), 0);
                    auto b = NativeInteger(int(Block[bit_idx]));
                    Output[i][j][k] = std::make_shared<LWECiphertextImpl>(zero_vec, b);
                }
            }
        }

        if (decrypt) {
            // First step
            HAES::AddRoundKey(Output, n_rounds - 1);
            HAES::ShiftRowsInv(Output);
            HAES::SubBytesInv(Output);

            // First n - 1 rounds

            for(int i = 1; i < n_rounds - 1; i++) {
                auto round_start = high_resolution_clock::now();
                HAES::AddRoundKey(Output, n_rounds - i - 1);
                HAES::MixColumnsInv(Output);
                HAES::ShiftRowsInv(Output);
                HAES::SubBytesInv(Output);
                auto round_stop = high_resolution_clock::now();
                auto elapsed = duration_cast<milliseconds>(round_stop-round_start).count();
                std::cout << "[!] Round" << i << " took " << elapsed << "ms." << std::endl;
            }

            // Finalize
            HAES::AddRoundKey(Output, 0);
        } else {
            // First step
            HAES::AddRoundKey(Output, 0);

            // First n - 1 rounds
            for(int i = 1; i < n_rounds - 1; i++) {
                HAES::SubBytes(Output);
                HAES::ShiftRows(Output);
                HAES::MixColumns(Output);
                HAES::AddRoundKey(Output, i);
            }

            // Finalize
            HAES::SubBytes(Output);
            HAES::ShiftRows(Output);
            HAES::AddRoundKey(Output, n_rounds - 1);
        }


    }


    void HAES::ShiftRows(BitTensor& state) {

        for (int i = 0; i < 4; i++) {
            for(int j = 0; j < i; j++) {
                for(int k = 0; k < 8; k++) {
                    // the inner part represents a shift to the left by 1
                    // so, for row i we shift 1+1+...+1 = i times
                    std::swap(state[i][0][k],state[i][3][k]);
                    std::swap(state[i][0][k],state[i][1][k]);
                    std::swap(state[i][1][k],state[i][2][k]);
                }
            }
        }
    }

    void HAES::ShiftRowsInv(AESTransciphering::BitTensor &state) {
        for (int i = 0; i < 4; i++) {
            for(int j = 0; j < i; j++) {
                for(int k = 0; k < 8; k++) {
                    // the inner part represents a shift to the right by 1
                    // so, for row i we shift 1+1+...+1 = i times
                    std::swap(state[i][1][k],state[i][2][k]);
                    std::swap(state[i][0][k],state[i][1][k]);
                    std::swap(state[i][0][k],state[i][3][k]);
                }
            }
        }
    }

    void HAES::AddRoundKey(AESTransciphering::BitTensor &State, int idx) {
        auto round_key = round_keys.at(idx);
        for(int i = 0; i < 4; i++) {
            auto& state_row = State.at(i);
            auto& round_key_row = round_key.at(i);
            for(int j = 0; j < 4; j++) {
                auto& state_col = state_row.at(j);
                auto& round_key_col = round_key_row.at(j);
                for(int k = 0; k < 8; k++) {
                    XOREQ(state_col[k], round_key_col[k]);
                }
            }
        }
    }

    void HAES::MixColumnsInv(AESTransciphering::BitTensor &State) {

        uint8_t scalars[] = {0xe,0xb,0xd,0x9};
        ByteVec acc;

        //lbcrypto::LWECiphertext zero = std::make_shared<lbcrypto::LWECiphertextImpl>(bits[0]->GetA() * 0, bits[0]->GetB() * 0);

        std::array<ByteVec, 4> new_entries;
        for(uint32_t i = 0; i < 4; i++) {
            ByteVec& c0 = State[0][i];
            ByteVec& c1 = State[1][i];
            ByteVec& c2 = State[2][i];
            ByteVec& c3 = State[3][i];

            // create outputs
            for(uint32_t j = 0; j < 4; j++) {
                auto c01 = GMUL(scalars[((4 - j))%4], c0);
                auto c02 = GMUL(scalars[((4 - j) +1)%4], c1);
                auto c03 = GMUL(scalars[((4 - j)+2)%4], c2);
                auto c04 = GMUL(scalars[((4 - j)+3)%4], c3);

                new_entries[j] = XOR_VEC(c01, c02);
                XOR_VEC_EQ(new_entries[j], c03);
                XOR_VEC_EQ(new_entries[j], c04);
            }

            for(uint32_t j = 0; j < 4; j++)
                State[j][i] = new_entries[j];
        }

    }

    void HAES::MixColumns(AESTransciphering::BitTensor &State) {

        uint8_t scalars[] = {0x2,0x3,0x1,0x1};
        ByteVec acc;

        //lbcrypto::LWECiphertext zero = std::make_shared<lbcrypto::LWECiphertextImpl>(bits[0]->GetA() * 0, bits[0]->GetB() * 0);

        std::array<ByteVec, 4> new_entries;
        for(uint32_t i = 0; i < 4; i++) {
            ByteVec& c0 = State[0][i];
            ByteVec& c1 = State[1][i];
            ByteVec& c2 = State[2][i];
            ByteVec& c3 = State[3][i];

            // create outputs

            for(uint32_t j = 0; j < 4; j++) {

                auto c01 = GMUL(scalars[((4 - j))%4], c0);
                auto c02 = GMUL(scalars[((4 - j) +1)%4], c1);
                auto c03 = GMUL(scalars[((4 - j)+2)%4], c2);
                auto c04 = GMUL(scalars[((4 - j)+3)%4], c3);

                new_entries[j] = XOR_VEC(c01, c02);
                XOR_VEC_EQ(new_entries[j], c03);
                XOR_VEC_EQ(new_entries[j], c04);
            }

            for(uint32_t j = 0; j < 4; j++)
                State[j][i] = new_entries[j];
        }
    }

    void HAES::SubBytes(AESTransciphering::BitTensor &state) {
        for(uint32_t i = 0; i < 4; i++)
            for(uint32_t j = 0; j < 4; j++) {
                std::vector<LWECiphertext> bits(state[i][j].begin(), state[i][j].end());
                auto sbox_res = evaluator->EvalLUT(bits,LUTEval::FINALIZED,LUTEval::FINALIZED);
                std::copy_n(sbox_res.begin(),8,state[i][j].begin());
            }
    }

 void HAES::SubBytesInv(AESTransciphering::BitTensor &state) {
        auto eval = *evaluator;
        std::vector<LWECiphertext> bits(8);

        // Depending on parameters, change FINALIZED and modulus switch/ key switch later
        for(uint32_t i = 0; i < 4; i++) {
            auto& stateI = state[i];
            for(uint32_t j = 0; j < 4; j++) {
                auto& stateJ = stateI[j];
                std::copy_n(stateJ.begin(), 8, bits.begin());
                auto sbox_res = eval.EvalLUT(bits, LUTEval::FINALIZED, LUTEval::FINALIZED);
                std::copy_n(sbox_res.begin(), 8, stateJ.begin());
            }
        }
    }
}
