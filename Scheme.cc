#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <gmp.h>
#include <cmath>
#include "Sampling.h"
#include "params.h"
#include "FFT.h"
#include "Random.h"
#include "Algebra.h"


using namespace std;
using namespace NTL;

const ZZX phi = Cyclo();





//GPV Sampling algorithm. 
// result is stored in vector v
// ✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅
void GPV(RR_t * v, const RR_t * const c, const RR_t s, const SK_Data * const SKD)
{
    int i;
    unsigned j;
    RR_t ci[2*N0], zi, cip, sip, aux;

    for(j=0; j<2*N0;j++)
    {
        ci[j] = c[j];
    }
 

    for(i=2*N0-1; i>=0; i--)
    {
        aux = (SKD->GS_Norms)[i];
        cip = DotProduct(ci, SKD->Bstar[i])/(aux*aux);
        sip = s/aux;
        zi = Sample4(cip, sip*PiPrime);

        for(j=0; j<2*N0; j++)
        {
            ci[j] -= zi*(SKD->B)[i][j];
        }
    }

    for(j=0; j<2*N0; j++)
    {
        v[j] = c[j] - ci[j];
    }

}


//==============================================================================
//Generates from parameters N and q :
// - a public key : polynomial h
// - a private key : polynomials f,g,F,G
// ✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅
//==============================================================================
void Keygen(ZZ_pX& PublicKey, ZZX* PrivateKey)
{
    ZZ SqNorm;
    ZZX f,g,F,G;

    SqNorm = conv<ZZ>(1.36*q0/2);

    GenerateBasis(f, g, F, G, SqNorm);
    PrivateKey[0] = f;
    PrivateKey[1] = g;
    PrivateKey[2] = F;
    PrivateKey[3] = G;

    for(unsigned int i=0; i<4; i++)
    {
            PrivateKey[i].SetLength(N0);
    }

    PublicKey = Quotient(f, g);
}

//==============================================================================
//Computes the private basis B from private key PrivateKey and parameter N
// ✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅
//==============================================================================
void CompletePrivateKey(mat_ZZ& B, const ZZX * const PrivateKey)
{
    ZZX f,g,F,G;
    f = PrivateKey[0];
    g = PrivateKey[1];
    F = PrivateKey[2];
    G = PrivateKey[3];

    f = -f;
    F = -F;

    B = BasisFromPolynomials(g, f, G, F);
}





//==============================================================================
//==============================================================================
//                            MAIN PROGRAMS
//==============================================================================
//==============================================================================

// stores privateKey and it's Matrix(B and BStar), privateKeyFFT, norm of BStar(gramShcmidth of B), std_dev of Gaussian 
// result is stored in SKD
// ✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅
void CompleteSK(SK_Data * SKD, ZZX * SK)
{
    unsigned int i, j;
    mat_ZZ B0;

    for(i=0; i<4; i++)
    {
        SKD->PrK[i] = SK[i];
        ZZXToFFT(SKD->PrK_fft[i], SK[i]);
    }

    CompletePrivateKey(B0, SK);

    for(i=0; i<2*N0; i++)
    {
        for(j=0; j<2*N0; j++)
        {
            SKD->B[i][j] = ( (RR_t) conv<double>(B0[i][j]) );
        }
    }

    for(i=0; i<1; i++)
    {
        FastMGS(SKD->Bstar, SKD->B);
    }

    for(i=0; i<2*N0; i++)
    {
        SKD->GS_Norms[i] = sqrt( DotProduct(SKD->Bstar[i], SKD->Bstar[i]) );
    }

    SKD->sigma = 2*SKD->GS_Norms[0];

}


// stores h and h_FFT in PKD
// ✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅
void CompletePK(PK_Data * PKD, ZZ_pX PK)
{
    PKD->h = PK;
    ZZXToFFT(PKD->h_FFT, conv<ZZX>(PK));
}



// ✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅
void PECKS_Trapdoor(ZZX TD_w[2], int * I, vec_ZZ * w, size_t m, const SK_Data * const SKD)
{
    // looping variable: (for)
    unsigned int i, j;

    // c: center of Gaussian Distribution, td: trapdoor array, sigma: std deviation of gaussian distribution
    RR_t c[2*N0], td[2*N0], sigma;
    // f,g are components of private key (i.e., of Basis B)
    ZZX f,g;

    // load f,g from SKD
    f = SKD -> PrK[0];
    g = SKD -> PrK[1];

    // load std. deviation of distribution from SKD
    sigma = SKD->sigma;
    // set the length of Trapdoor vector: TD_w = (s, t_w)
    TD_w[0].SetLength(N0);
    TD_w[1].SetLength(N0);

    // set the center of Gaussian distribution for GPV sampling: (note that t = H(w) is implemented by w mod q0
    for(i=0;i<N0;i++)
    {
        c[i] = 0;
        for(j=0; j<m; j++) {
            
            c[i] += ((RR_t)conv<double>(w[I[j]][i]));
        }
        c[i] = fmodl(c[i], q0);
        c[i+N0] = 0;
    }
    // sample a gaussian vector (td) from lattice:
    GPV(td, c, sigma, SKD);
    
    // (s, t_w) = (t, 0) - GPV(B, sigma, (t, 0)):
    for(i=0; i<N0; i++)
    {
        td[i] = c[i] - td[i];
        td[i+N0] = - td[i+N0];
    }
    
    // save s, t_w in TD_w vector (the trapdoor is t_w but we save s too for testing)
    for(i=0; i<N0; i++)
    {
        TD_w[0][i] = td[i];
        TD_w[1][i] = td[i+N0];
    }
    
}

// verify if the GPV-Sampled vector has the property s + t_w*h = t
// ✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅
unsigned long PECKS_Verify_Trapdoor(const ZZX TD_w[2], int * I, vec_ZZ * w, size_t m, const SK_Data * const SKD)
{
    unsigned int i;
    ZZX f,g,t,aux;

    f = SKD -> PrK[0];
    g = SKD -> PrK[1];
    
    t = conv<ZZX>(w[I[0]]);
    for(i=1; i<m; i++)
    {
        t = t + conv<ZZX>(w[I[i]]);
    }
    // note that h = g*(f^-1) mod q
    // aux = ((s-t)*f + g*t_w) mod phi
    aux = ((TD_w[0] - t)*f + g*TD_w[1])%phi;

    for(i=0; i<N0; i++)
    {
        aux[i] %= q1;
    }

    if( IsZero(aux) != 0)
    {
        cout << "The signature (s, t_w) doesn't verify the required equality [ (s-t)*f + g*t_w = 0 ] !\nActually, (s-t)*f + g*t_w = " << aux << endl << endl;
    }
    return IsZero(aux);
}

// Compute the PEKS of the keyword.
// the result is saved in SE[2]
// ✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅
void PECKS_Peck(long SE[1+2*l0][N0], vec_ZZ * w, const PK_Data * const PKD)
{
    unsigned long i, j;
    
    // variables for small random vectors r, e1, e1
    long r[N0], e1[N0], e2[l0][N0];
    // keyword enkryption is done via KEM mechanism. key = k:
    // t = Hash of keyword w
    long k[N0], t[l0][N0];
    // variables for fft of r, t, and aux variables
    CC_t r_FFT[N0], t_FFT[l0][N0], aux1_FFT[N0], aux2_FFT[l0][N0];
    for(j=0;j<l0;j++)
    {
        for(i=0;i<N0;i++)
        {
            t[j][i] = ((long) conv<double>(w[j][i])) ;
        }
    }
    
    
    // sample r, e1, e2 from {-1, 0, 1}^N and k from {0,1}^N
    for(i=0; i<N0; i++)
    {
        e1[i] = (rand()%3) - 1;
        for(j=0;j<l0;j++)
        {
            e2[j][i] = (rand()%3) - 1;
        }
        k[i] = rand()%2;
        r[i] = (rand()%3) - 1;
    }

    // compute FFT of r and save it into r_FFT
    MyIntFFT(r_FFT, r);
    // compute FFT of t and save it into t_FFT
    for(j=0;j<l0;j++)
    {
        MyIntFFT(t_FFT[j], t[j]);
    }
    

    // calculate r_FTT * h_FFT and r_FFT * t_FFT
    for(i=0; i<N0; i++)
    {
        aux1_FFT[i] = r_FFT[i]*((PKD->h_FFT)[i]);
        for(j=0;j<l0;j++)
        {
            aux2_FFT[j][i] = r_FFT[i]*t_FFT[j][i];
        }
    }
    // find the IFFT(r_FTT * h_FFT) which is r*h . save it into SE[0]
    MyIntReverseFFT(SE[0], aux1_FFT);
    // find the IFFT(r_FTT * t_FFT) which is r*t . save it into SE[1]
    for(j=0;j<l0;j++) 
    {
        MyIntReverseFFT(SE[1+j], aux2_FFT[j]);
    }

    // ok now we have SE[0] = r*h, SE[1] = r*t
    
    for(i=0; i<N0; i++)
    {
        SE[0][i] = (SE[0][i] + e1[i]               + q0/2)%q0 - (q0/2);
        for(j=0;j<l0;j++)
        {
            SE[1+j][i] = (SE[1+j][i] + e2[j][i] + (q0/2)*k[i] + q0/2)%q0 - (q0/2);
            SE[1+l0+j][i] = (SE[1+j][i] + k[i]                + q0/2)%q0 - (q0/2);// I used a simple mod Hash function
        }
    } 

}



// ✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅
bool PECKS_Test(const PK_Data * const PKD, long SE[1+2*l0][N0], int * I, size_t m, ZZX t_w)
{
    unsigned int i, j;
    CC_t b_FFT[N0], tw_FFT[N0], aux_FFT[N0];
    long y[N0];
    long sum_of_ci[N0];
    long sum_of_di[N0];
    long sum_of_hashes[N0];

    MyIntFFT(b_FFT, SE[0]);
    ZZXToFFT(tw_FFT, t_w);
    for(i=0; i<N0; i++)
    {
        aux_FFT[i] = b_FFT[i]*tw_FFT[i];
        sum_of_ci[i] = 0;
        sum_of_di[i] = 0;
        sum_of_hashes[i] = 0;
    }

    MyIntReverseFFT(y, aux_FFT);
    // now y = b * t_w 
    

    for(i=0; i<N0; i++)
    {
        for (j = 0; j < m; j++)
        {
            sum_of_ci[i] += SE[1+I[j]][i];
            sum_of_di[i] += SE[1+l0+I[j]][i];
        }
        y[i] = ((unsigned long)(sum_of_ci[i] - y[i])) % q0; // now y = C[I_1] + ... + C[I_m] - b * t_w 
        y[i] = (y[i] + ((m*q0)/4)) / ((m*q0)/2);
        y[i] %= 2;

        for (j = 0; j < m; j++)
        {
            sum_of_hashes[i] = (sum_of_hashes[i]+ SE[1+I[j]][i] + y[i]                + q0/2)%q0 - (q0/2);
        }
        if (sum_of_di[i] != sum_of_hashes[i])
        {
            return false;
        }
    }
    
    return true;
}





//==============================================================================
//==============================================================================
//                             BENCHES AND TESTS
//                   FOR EXTRACTION AND ENCRYPTION/DECRYPTION
//==============================================================================
//==============================================================================


float Trapdoor_Bench(const unsigned int nb_trap, SK_Data * SKD)
{
    clock_t t1, t2;
    float diff;
    unsigned int i;
    unsigned int m = l0;
    int I[m];
    for (size_t j = 0; j < m; j++)
    {
        I[j] = rand() % l0;
    }
    vec_ZZ w[m];
    ZZX TD_w[2];

    
    t1 = clock();

    for(i=0; i<nb_trap; i++)
    {
        for (size_t j = 0; j < m; j++)
        {
            w[j] = RandomVector();
        }    
        PECKS_Trapdoor(TD_w, I, w, m, SKD);
    }

    t2 = clock();
    diff = ((float)t2 - (float)t1)/1000000.0F;
    cout << "\n\nIt took " << diff << " seconds to generate " << nb_trap << " trapdoors." << endl;
    cout << "That's " << (diff/nb_trap)*1000 << " milliseconds per trapdoor." << endl << endl;

    return (diff/nb_trap)*1000;
}


float Peck_Bench(const unsigned int nb_peck, PK_Data * PKD)
{
    clock_t t1, t2;
    double diff;
    unsigned int i;
    vec_ZZ w[l0];
    long SE[1+2*l0][N0];

    // ZZX TD_w;
    
    t1 = clock();

    for(i=0; i<nb_peck; i++)
    {
        for(size_t j = 0; j < l0; j++)
        {
            w[j] = RandomVector();
        }
        PECKS_Peck(SE, w, PKD);
    }

    t2 = clock();
    diff = ((double)t2 - (double)t1)/1000000.0l;
    cout << "\n\nIt took " << diff << " seconds to do " << nb_peck << " Peck operations." << endl;
    cout << "That's " << (diff/nb_peck)*1000 << " milliseconds per Peck." << endl;
    cout << "That's " << (diff/nb_peck)*1000*1024/N0 << " milliseconds per Peck per Kilobit." << endl << endl;
    return (diff/nb_peck)*1000;
}

float Test_Bench(const unsigned int nb_test, PK_Data * PKD, SK_Data * SKD)
{
    clock_t t1, t2;
    double diff;
    unsigned int i;
    
    // first generate some keywords, trapdoors and Searchable Encryptions
    size_t m = l0;
    int I[m];
    for (size_t j = 0; j < m; j++)
    {
        I[j] = rand() % l0;
    }
    vec_ZZ w[m];
    ZZX TD_w[2];

        
    long SE[1+2*l0][N0];
    
    
    // Now do the Tests:

    t1 = clock();

    for(i=0; i<nb_test; i++)
    {
        for (size_t j = 0; j < m; j++)
        {
            w[j] = RandomVector();
        }    
        PECKS_Peck(SE, w, PKD);
        PECKS_Trapdoor(TD_w, I, w, m, SKD);
        PECKS_Test(PKD, SE, I, m, TD_w[1]);
    }

    t2 = clock();
    diff = ((double)t2 - (double)t1)/1000000.0l;
    cout << "\n\nIt took " << diff << " seconds to do " << nb_test << " Test operations." << endl;
    cout << "That's " << (diff/nb_test)*1000 << " milliseconds per Test." << endl;
    cout << "That's " << (diff/nb_test)*1000*1024/N0 << " milliseconds per Test per Kilobit." << endl << endl;
    return (diff/nb_test)*1000;
}


unsigned int Trapdoor_Test(const unsigned int nb_trap, SK_Data * SKD)
{
    unsigned int i, rep;
    
    size_t m = l0;
    int I[m];
    for (size_t j = 0; j < m; j++)
    {
        I[j] = rand() % l0;
    }
    vec_ZZ w[m];
    ZZX TD_w[2];


    rep = 0;

    for(i=0; i<nb_trap; i++)
    {
        for (size_t j = 0; j < m; j++)
        {
            w[j] = RandomVector();
        }    
        PECKS_Trapdoor(TD_w, I, w, m, SKD);
        rep += PECKS_Verify_Trapdoor(TD_w, I, w, m, SKD);
    }

    cout << endl;
    if(rep == 0){
        cout << endl << nb_trap << " Trapdoor generations successfully performed!" << endl << endl;   
    }
    else {
        cout << endl << rep << " out of " << nb_trap << " Trapdoor generations failed miserabily!" << endl << endl;
    }
    return rep;
}


unsigned int Peck_Test(const unsigned int nb_peck, PK_Data * PKD, SK_Data * SKD)
{
    unsigned int i, rep;
    bool test_result = false;
    size_t m = l0;
    int I[m];
    for (size_t j = 0; j < m; j++)
    {
        I[j] = rand() % l0;
    }
    vec_ZZ w[m];
    ZZX TD_w[2];

    long SE[1+2*l0][N0];

    
    rep = 0;

    for(i=0; i<nb_peck; i++)
    {
        for (size_t j = 0; j < m; j++)
        {
            w[j] = RandomVector();
        }    
        PECKS_Trapdoor(TD_w, I, w, m, SKD);
        PECKS_Peck(SE, w, PKD);
        if (PECKS_Verify_Trapdoor(TD_w, I, w, m, SKD) > 0) {
            cout << "Bad Trapdoors generated!! \n";
        }
        test_result = PECKS_Test(PKD, SE, I, m, TD_w[1]);
        if(!test_result)
        {
            // cout << "ERROR : Test = false " << endl;
            rep++;
        }
        

    }

    cout << endl;
    if(rep == 0)
    {    cout << endl << nb_peck << " Tests successfully performed!" << endl << endl;    }
    else
    {    cout << endl << rep << " out of " << nb_peck << " Tests failed miserabily!" << endl << endl;    }
    return rep;
}
