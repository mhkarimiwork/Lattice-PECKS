#ifndef LIBE_SCHEME_H
#define LIBE_SCHEME_H

#include "params.h"
#include "Sampling.h"
#include "FFT.h"
#include "Random.h"
#include "Algebra.h"

void GPV(RR_t * v, const RR_t * const c, const RR_t s, const SK_Data * const SKD);

void Keygen(ZZ_pX& PublicKey, ZZX* PrivateKey);
void CompletePrivateKey(mat_ZZ& B, const ZZX * const PrivateKey);
void CompleteSK(SK_Data * SKD, ZZX * SK);
void CompletePK(PK_Data * PKD, ZZ_pX PK);

void PECKS_Trapdoor(ZZX TD_w[2], int * I, vec_ZZ * w, size_t m, const SK_Data * const SKD);
unsigned long PECKS_Verify_Trapdoor(const ZZX TD_w[2], int * I, vec_ZZ * w, size_t m, const SK_Data * const SKD);
void PECKS_Peck(long SE[1+2*l0][N0], vec_ZZ * w, const PK_Data * const PKD);
bool PECKS_Test(const PK_Data * const PKD, long SE[1+2*l0][N0], int * I, size_t m, ZZX t_w);

float Trapdoor_Bench(const unsigned int nb_trap, SK_Data * SKD);
float Peck_Bench(const unsigned int nb_peck, PK_Data * PKD);
float Test_Bench(const unsigned int nb_test, PK_Data * PKD, SK_Data * SKD);


unsigned int Trapdoor_Test(const unsigned int nb_trap, SK_Data * SKD);
unsigned int Peck_Test(const unsigned int nb_peck, PK_Data * PKD, SK_Data * SKD);


#endif
