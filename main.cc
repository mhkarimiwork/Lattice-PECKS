#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <gmp.h>
#include <fstream>

#include "params.h"
#include "io.h"
#include "FFT.h"
#include "Sampling.h"
#include "Random.h"
#include "Algebra.h"
#include "Scheme.h"


using namespace std;
using namespace NTL;




//==============================================================================
//==============================================================================
//                                  MAIN
//==============================================================================
//==============================================================================


int main()
{
    cout << "\n=======================================================================\n";
    cout << "This program is a proof-of concept for efficient PECKS over lattices.\n";
    cout << "It generates a NTRU lattice of dimension 2N and associated modulus q,\n";
    cout << "and perform benches and tests, for trapdoor generation, PECK and Test Operations.";
    cout << "\n=======================================================================\n\n";
    
    // open a file for writing the result:
    ofstream myfile;
    myfile.open ("results.txt");
    myfile <<  "===================================================================\n";
    myfile << "N = " << N0 <<endl;
    myfile << "q = " << q0 <<endl;
    myfile << "l = " << l0 <<endl;
    myfile <<  "===================================================================\n\n\n";



    ZZX SK[4];
    ZZ_pX phiq, PK;
    unsigned int i;
    float diff;
    SK_Data * SKD = new SK_Data;
    PK_Data * PKD = new PK_Data;
    clock_t t1, t2;
    const ZZX phi = Cyclo();

    srand(rdtsc()); // initialisation of rand

    cout << "N = " << N0 << endl;
    cout << "q = " << q0 << endl;
    cout << "l = " << l0 << endl;

    ZZ_p::init(q1);
    zz_p::init(q0);

    phiq = conv<ZZ_pX>(phi);
    ZZ_pXModulus PHI(phiq);


    cout << "\n===================================================================\n KEY GENERATION";
    cout << "\n===================================================================\n";
    t1 = clock();
    for(i=0; i<1; i++)
    {
        Keygen(PK, SK);
    }
    t2 = clock();

    CompleteSK(SKD, SK);
    CompletePK(PKD, PK);


    diff = ((float)t2 - (float)t1)/1000000.0F;
    cout << "It took " << diff << " seconds to generate the  key pair" << endl;
    myfile << "Key pair generation time:   " << diff << " seconds" << endl;


    unsigned int TRIES = 10;
    //==============================================================================
    // Trapdoor Generation bench and PECK/TEST bench
    //==============================================================================
    cout << "\n===================================================================\n RUNNING PECKS.Trapdoor BENCH FOR ";
    cout << TRIES << " DIFFERENT Keywords\n===================================================================\n";
    myfile << "Trapdoor(sk, W) avg time:   " << Trapdoor_Bench(TRIES, SKD) << " milliseconds" << endl;
    

    cout << "\n===================================================================\n RUNNING PECKS.PECK BENCH FOR ";
    cout << TRIES << " DIFFERENT KEYWORDS\n===================================================================\n";
    myfile << "PECKS(pk, Q) avg time:      " << Peck_Bench(TRIES, PKD) << " milliseconds" << endl;
    

    cout << "\n===================================================================\n RUNNING PECKS.Test BENCH FOR ";
    cout << TRIES << " DIFFERENT KEYWORDS\n===================================================================\n";
    myfile << "Test(pk, SE, T_Q) avg time: " << Test_Bench(TRIES, PKD, SKD) << " milliseconds" << endl;
    
    
    myfile <<  "===================================================================\n\n\n";
    //==============================================================================
    // Trapdoor generation test and PECK/Test test
    //==============================================================================
    cout << "\n===================================================================\n CHECKING TRAPDOOR VALIDITY FOR ";
    cout << TRIES << " DIFFERENT KEYWORDS\n===================================================================\n";
    myfile << "Invalid Trapdoors in " << TRIES << " tries: " << Trapdoor_Test(TRIES, SKD) << endl;

    cout << "\n===================================================================\n CHECKING Trapdoor/PECK Validity with Test Algorithm FOR ";
    cout << TRIES << " DIFFERENT KEYWORDS\n===================================================================\n";
    myfile << "Invalid Trapdoor/PECK in " << TRIES << " tries: " << Peck_Test(TRIES, PKD, SKD) << endl;
    

    free(SKD);
    free(PKD);
    myfile.close();
    return 0;
}
