#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <gmp.h>


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

    CompleteSK(SKD, SK);
    CompletePK(PKD, PK);

    t2 = clock();

    diff = ((float)t2 - (float)t1)/1000000.0F;
    cout << "It took " << diff << " seconds to generate the  key pair" << endl;



    //==============================================================================
    // Trapdoor Generation bench and PECK/TEST bench
    //==============================================================================
    const unsigned int nb_trap = 100;
    const unsigned int nb_peck = 1000;

    cout << "\n===================================================================\n RUNNING TRAPDOOR BENCH FOR ";
    cout << nb_trap << " DIFFERENT Keywords\n===================================================================\n";
    Trapdoor_Bench(nb_trap, SKD);

    cout << "\n===================================================================\n RUNNING PECK BENCH FOR ";
    cout << nb_peck << " DIFFERENT KEYWORDS\n===================================================================\n";
    Peck_Bench(nb_peck, PKD);



    ///==============================================================================
    // Trapdoor generation test and PECK/Test test
    //==============================================================================
    cout << "\n===================================================================\n CHECKING TRAPDOOR VALIDITY FOR ";
    cout << nb_trap << " DIFFERENT KEYWORDS\n===================================================================\n";
    Trapdoor_Test(nb_trap, SKD);

    cout << "\n===================================================================\n CHECKING PECK/TEST VALIDITY FOR ";
    cout << nb_trap << " DIFFERENT KEYWORDS\n===================================================================\n";
    Peck_Test(nb_peck, PKD, SKD);

    free(SKD);
    free(PKD);
    return 0;
}
