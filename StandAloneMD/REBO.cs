using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace StandAloneMD
{
    class REBO : Potential
    {
        private double[,] SPGC, SPGH;
        private int[] IGC, IGH;
        private int[,] IN3;
        private double[, ,] XH, XH1, XH2;
        private double ATT, XQM;
        private double[, , ,] CLM;
        private double[, , , ,] CLMN;

        public REBO()
        {
            setParameters();
        }

        public override void preCompute()
        {
            
        }

        //the function returns the LennarJones force on the atom given the list of the atoms that are within range of it
        private float calcAcceleration(float distance, Atom firstAtom, Atom secondAtom)
        {
            
            return 0.0f;
        }

        //the function returns the LennarJones force on the atom given the list of the atoms that are within range of it
        private float calcPotential(float distance, Atom firstAtom, Atom secondAtom)
        {
            
            return 0.0f;
        }

        //the function returns the Lennard-Jones force on the atom given the list of all the atoms in the simulation
        public override void getForce(Atom firstAtom, Atom secondAtom)
        {
            
        }

        //the function returns the Lennard-Jones force on the atom given the list of all the atoms in the simulation
        public override float getPotential(Atom firstAtom, Atom secondAtom)
        {
            return 0.0f;
        }

        public override void calculateVerletRadius()
        {
            
        }

        //This function creates a list of all neighbor list for each atom
        public override void calculateNeighborList()
        {
            
        }

        private void setParameters()
        {
            string line;

            //Parameters for hydrocarbons
            SPGC = new double[6,5] { { 0.2817216000000E+00, 0.2817216000000E+00, 0.6900668660000E+00, 0.3754490870000E+00, 0.2718558000000E+00 }, { 0.1062912000000E+01, 0.1062912000000E+01, 0.5460691360000E+01, 0.1407252749388E+01, 0.4892727456293E+00 }, { 0.2136736000000E+01, 0.2136736000000E+01, 0.2301345680000E+02, 0.2255103926323E+01, -0.4328199017473E+00 }, { 0.2533952000000E+01, 0.2533952000000E+01, 0.5491519344000E+02, 0.2028902219952E+01, -0.5616795197048E+00 }, { 0.1554736000000E+01, 0.1554736000000E+01, 0.6862037040000E+02, 0.1426981217906E+01, 0.1270874966906E+01 }, { 0.3863296000000E+00, 0.3863296000000E+00, 0.3470897779200E+02, 0.5063107994308E+00, -0.3750409108350E-01 } };
            IGC = new int[25] {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2, 2, 1, 1, 1, 1, 1};
            SPGH = new double[6, 3] { { 270.467795364007301, 16.956325544514659, 19.065031149937783 }, { 1549.701314596994564, -21.059084522755980, 2.017732531534021 }, { 3781.927258631323866, -102.394184748124742, -2.566444502991983 }, { 4582.337619544424228, -210.527926707779059, 3.291353893907436 }, { 2721.538161662818368, -229.759473570467513, -2.653536801884563 }, { 630.658598136730774, -94.968528666251945, 0.837650930130006 } };
            IGH = new int[25] { 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1};
            
            XH = new double[2, 10, 10];
            Array.Clear(XH,0,200);
            XH1 = new double[2, 10, 10];
            Array.Clear(XH1, 0, 200);
            XH2 = new double[2, 10, 10];
            Array.Clear(XH2, 0, 200);
            
            ATT = 3.20E0;
            XQM = 3.70E0;

            //Zero bicubic spline coefficients
            CLM = new double[2, 10, 10, 16];
            Array.Clear(CLM, 0, 2 * 10 * 10 * 16);

            // open Spline file inter2d_iv.d
            StreamReader inter2d_iv = new StreamReader("../../../Spline/inter2d_iv.d");;

            //bicubic spline
            line = inter2d_iv.ReadLine();
            int I2D = int.Parse(line);

            //read integer values of bicubic spline
            while (line != null)
            {
                line = inter2d_iv.ReadLine();
                if (line != null)
                {
                    string[] words = line.Split(" ".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);

                    int I = int.Parse(words[0]);
                    if (I <= 0)  break; 

                    int J = int.Parse(words[1]);
                    int K = int.Parse(words[2]);
                    XH[I-1, J-1, K-1] = double.Parse(words[3]);
                    //Console.WriteLine(XH[I - 1, J - 1, K - 1]);
                }
            }

            XH1[2 - 1, 3 - 1, 1 - 1] = (XH[2 - 1, 4 - 1, 1 - 1] - XH[2 - 1, 2 - 1, 1 - 1]) / 2.0;
            XH1[2 - 1, 2 - 1, 2 - 1] = (XH[2 - 1, 3 - 1, 2 - 1] - XH[2 - 1, 1 - 1, 2 - 1]) / 2.0;

            XH2[2 - 1, 2 - 1, 2 - 1] = (XH[2 - 1, 2 - 1, 3 - 1] - XH[2 - 1, 2 - 1, 1 - 1]) / 2.0;
            XH2[2 - 1, 1 - 1, 3 - 1] = (XH[2 - 1, 1 - 1, 4 - 1] - XH[2 - 1, 1 - 1, 2 - 1]) / 2.0;

            //read bicubic spline coefficients
            while (line != null)
            {
                line = inter2d_iv.ReadLine();
                if (line != null)
                {
                    string[] words = line.Split(" ".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);

                    int I = int.Parse(words[0]);
                    int L = int.Parse(words[1]);
                    int M = int.Parse(words[2]);

                    line = inter2d_iv.ReadLine();
                    words = line.Split(" ".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                    CLM[I - 1, L - 1, M - 1, 0] = double.Parse(words[0]);
                    CLM[I - 1, L - 1, M - 1, 1] = double.Parse(words[1]);
                    CLM[I - 1, L - 1, M - 1, 2] = double.Parse(words[2]);
                    CLM[I - 1, L - 1, M - 1, 3] = double.Parse(words[3]);

                    line = inter2d_iv.ReadLine();
                    words = line.Split(" ".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                    CLM[I - 1, L - 1, M - 1, 4] = double.Parse(words[0]);
                    CLM[I - 1, L - 1, M - 1, 5] = double.Parse(words[1]);
                    CLM[I - 1, L - 1, M - 1, 6] = double.Parse(words[2]);
                    CLM[I - 1, L - 1, M - 1, 7] = double.Parse(words[3]);

                    line = inter2d_iv.ReadLine();
                    words = line.Split(" ".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                    CLM[I - 1, L - 1, M - 1, 8] = double.Parse(words[0]);
                    CLM[I - 1, L - 1, M - 1, 9] = double.Parse(words[1]);
                    CLM[I - 1, L - 1, M - 1, 10] = double.Parse(words[2]);
                    CLM[I - 1, L - 1, M - 1, 11] = double.Parse(words[3]);

                    line = inter2d_iv.ReadLine();
                    words = line.Split(" ".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                    CLM[I - 1, L - 1, M - 1, 12] = double.Parse(words[0]);
                    CLM[I - 1, L - 1, M - 1, 13] = double.Parse(words[1]);
                    CLM[I - 1, L - 1, M - 1, 14] = double.Parse(words[2]);
                    CLM[I - 1, L - 1, M - 1, 15] = double.Parse(words[3]);

                }
            }

            // close bicubic spline file
            inter2d_iv.Close();

            //read tricubic spline coefficients
            IN3 = new int[64, 3];
            int ic = 0;
            for (int I=0; I<4; I++)
            {
                for (int J=0; J < 4; J++)
                {
                    for (int K=0; K < 4; K++)
                    {
                        ic++;
                        IN3[ic-1, 0] = I;
                        IN3[ic-1, 1] = J;
                        IN3[ic-1, 2] = K;
                    }
                }
            }


            //tricubic spline coefficient
            CLMN = new double[3, 10, 10, 10, 64];
            /******************************************************************************************************************************************************************/
            // open tricubic spline file inter3d_iv_new.d
            StreamReader inter3d_iv = new StreamReader("../../../Spline/inter3d_iv_new.d"); ;
            line = inter3d_iv.ReadLine();
            int I3D = int.Parse(line);
            //read tricubic spline coefficient CLMN[0,L,M,N,I]
            while (line != null)
            {
                line = inter3d_iv.ReadLine();
                if (line != null)
                {
                    string[] words = line.Split(" ".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                    int L = int.Parse(words[0]);
                    int M = int.Parse(words[1]);
                    int N = int.Parse(words[2]);
                    for (int I = 0; I < 21; I++ )
                    {
                        line = inter3d_iv.ReadLine();
                        words = line.Split(" ".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                        CLMN[0, L - 1, M - 1, N - 1, 3 * I + 0] = double.Parse(words[0]);
                        CLMN[0, L - 1, M - 1, N - 1, 3 * I + 1] = double.Parse(words[1]);
                        CLMN[0, L - 1, M - 1, N - 1, 3 * I + 2] = double.Parse(words[2]);
                    }
                    line = inter3d_iv.ReadLine();
                    words = line.Split(" ".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                    CLMN[0, L - 1, M - 1, N - 1, 63] = double.Parse(words[0]);
                }
            }
            inter3d_iv.Close();
            /******************************************************************************************************************************************************************/
            // open tricubic spline file inter3d_ch.d
            StreamReader inter3d_ch = new StreamReader("../../../Spline/inter3d_ch.d"); ;
            line = inter3d_ch.ReadLine();
            I3D = int.Parse(line);
            //read tricubic spline coefficient CLMN[1,L,M,N,I]
            while (line != null)
            {
                line = inter3d_ch.ReadLine();
                if (line != null)
                {
                    string[] words = line.Split(" ".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                    int L = int.Parse(words[0]);
                    int M = int.Parse(words[1]);
                    int N = int.Parse(words[2]);
                    for (int I = 0; I < 21; I++)
                    {
                        line = inter3d_ch.ReadLine();
                        words = line.Split(" ".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                        CLMN[1, L - 1, M - 1, N - 1, 3 * I + 0] = double.Parse(words[0]);
                        CLMN[1, L - 1, M - 1, N - 1, 3 * I + 1] = double.Parse(words[1]);
                        CLMN[1, L - 1, M - 1, N - 1, 3 * I + 2] = double.Parse(words[2]);
                    }
                    line = inter3d_ch.ReadLine();
                    words = line.Split(" ".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                    CLMN[1, L - 1, M - 1, N - 1, 63] = double.Parse(words[0]);
                }
            }
            inter3d_ch.Close();
            /******************************************************************************************************************************************************************/
            // open tricubic spline file inter3d_h.d
            StreamReader inter3d_h = new StreamReader("../../../Spline/inter3d_h.d"); ;
            line = inter3d_h.ReadLine();
            I3D = int.Parse(line);
            //read tricubic spline coefficient CLMN[2,L,M,N,I]
            while (line != null)
            {
                line = inter3d_h.ReadLine();
                if (line != null)
                {
                    string[] words = line.Split(" ".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                    int L = int.Parse(words[0]);
                    int M = int.Parse(words[1]);
                    int N = int.Parse(words[2]);
                    for (int I = 0; I < 21; I++)
                    {
                        line = inter3d_h.ReadLine();
                        words = line.Split(" ".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                        CLMN[2, L - 1, M - 1, N - 1, 3 * I + 0] = double.Parse(words[0]);
                        CLMN[2, L - 1, M - 1, N - 1, 3 * I + 1] = double.Parse(words[1]);
                        CLMN[2, L - 1, M - 1, N - 1, 3 * I + 2] = double.Parse(words[2]);
                    }
                    line = inter3d_h.ReadLine();
                    words = line.Split(" ".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                    CLMN[2, L - 1, M - 1, N - 1, 63] = double.Parse(words[0]);
                }
            }
            inter3d_h.Close();

            Console.WriteLine("end of setParameters function!!!");

        }
    }
}
