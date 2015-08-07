using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Diagnostics;

namespace StandAloneMD
{
    class REBO : Potential
    {
		private int NDIHED;
		private int NTAB = 10000; // this is the array size of potential lookup table.
        private int[] IGC, IGH;
		private int[] NABORS;
		private List<int> IVCT2B = new List<int>();
		private List<int> JVCT2B = new List<int>();
		private int[] LCHECK;
        private int[,] IN2, IN3;
		private double RLL = 0.5; // this is most probably the shell thickness used for neighbor list calculation
        private double SIGMA, EPSI, ATT, XQM, PQ, PIDT, XXDB;
		private double tote; // total potential energy
		private double[] XN1, XTN2, XTN1, ADB, CDB, CDB2, DDB, DDB2, HDB;
		private double[] eatom;
		private double[] RCOR, WW, DWW, EXX1, DEXX1;
		private double[,] AD, AXL, BD, BXL, CD, CXL, DD, DXL, ED, RB1, RB2, PID, RMAX, RLIST, SPGC, SPGH, CHI, DDTAB;
		private double[,] COR;
		private double[,] RNP; // this is the array that holds the value of forces acting on atoms
		private double[, ,] XH, XH1, XH2, XDB, REG;
		private double[, ,] TABFC, TABDFC, ATABLE, DATABLE, RTABLE, DRTABLE;
        private double[, , ,] CLM, TLMN;
        private double[, , , ,] CLMN;
		private bool neighborListFlag = true;

        public REBO()
        {
            setParameters();
			mtable();
			
        }

        public override void preCompute()
        {
			setParameters();
			mtable();
			neighborListFlag = true;
        }

        //the function returns the force on the atom given the list of the atoms that are within range of it
        private float calcAcceleration(float distance, Atom firstAtom, Atom secondAtom)
        {
            
            return 0.0f;
        }

        //the function returns the force on the atom given the list of the atoms that are within range of it
        private float calcPotential(float distance, Atom firstAtom, Atom secondAtom)
        {
            
            return 0.0f;
        }

        //the function returns the force on the atom given the list of all the atoms in the simulation
		public override void getForce(Atom firstAtom, Atom secondAtom)
        {
			caguts();
			pibond();

			for (int i = 0; i < Atom.AllAtoms.Count; i++)
			{
				Atom currAtom = Atom.AllAtoms[i];
				for (int j = 0; j < 3; j++)
				{
					if (double.IsNaN(RNP[i, j]))
					{
						Console.WriteLine("Calculated force is NaN!!!");
						Console.ReadLine();

					}
					currAtom.accelerationNew[j] = (float)(RNP[i, j] * 1.0e15) / currAtom.massamu;
				}
			}
        }

        //the function returns the force on the atom given the list of all the atoms in the simulation
        public override float getPotential(Atom firstAtom, Atom secondAtom)
        {
			tote = 0.0;
			return 0.0f;
			
        }

		// the function calculates the verlet radius that is used to make neighbor list
        public override void calculateVerletRadius()
        {
            
        }

        //This function creates a list of all neighbor list for each atom
        public override void calculateNeighborList()
        {
			neighborListFlag = true;
        }

        private void setParameters() //Set the potential coefficients and parameters for the potentials
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
            
			int ic = 0;
            IN2 = new int[16, 2];
            for (int I = 0; I < 4; I++ )
            {
                for (int J = 0; J<4; J++)
                {
                    ic++;
                    IN2[ic - 1, 0] = I;
                    IN2[ic - 1, 1] = J;
                }
            }

            //read tricubic spline coefficients
            IN3 = new int[64, 3];
            ic = 0;
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
			Array.Clear(CLMN, 0, 3 * 10 * 10 * 64);
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

            for (int L = 0; L < 10; L++ )
            {
				for (int M = 0; M < 10; M++)
				{
					for (int I = 0; I<64; I++)
					{
						CLMN[1 - 1, L, M, 10 - 1, I] = CLMN[1 - 1, L, M, 9 - 1, I];
						CLMN[2 - 1, L, M, 10 - 1, I] = CLMN[2 - 1, L, M, 9 - 1, I];
						for (int N = 5; N < 10; N++)
						{
							CLMN[3 - 1, L, M, N, I] = CLMN[3 - 1, L, M, 5 - 1, I];
						}
					}
				}
            }

			// read tricubic spline coefficients for torsional potential
			TLMN = new double[10, 10, 10, 64];
			Array.Clear(TLMN, 0, 10 * 10 * 10 * 64);
			/******************************************************************************************************************************************************************/
			// open tricubic spline file inter3d_iv_new.d
			StreamReader inter3dtors = new StreamReader("../../../Spline/inter3dtors.d"); ;
			line = inter3dtors.ReadLine();
			int ITD = int.Parse(line);
			while (line != null)
			{
				line = inter3dtors.ReadLine();
				if (line != null)
				{
					string[] words = line.Split(" ".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
					int L = int.Parse(words[0]);
					int M = int.Parse(words[1]);
					int N = int.Parse(words[2]);
					for (int I = 0; I < 21; I++)
					{
						line = inter3dtors.ReadLine();
						words = line.Split(" ".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
						TLMN[L - 1, M - 1, N - 1, 3 * I + 0] = double.Parse(words[0]);
						TLMN[L - 1, M - 1, N - 1, 3 * I + 1] = double.Parse(words[1]);
						TLMN[L - 1, M - 1, N - 1, 3 * I + 2] = double.Parse(words[2]);
					}
					line = inter3dtors.ReadLine();
					words = line.Split(" ".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
					TLMN[L - 1, M - 1, N - 1, 63] = double.Parse(words[0]);
				}
			}
			inter3dtors.Close();

			for (int L = 0; L < 10; L++ )
			{
				for (int M = 0; M<10; M++)
				{
					for (int N = 3; N<10; N++)
					{
						for (int I=0; I<64; I++)
						{
							TLMN[L, M, N, I] = TLMN[L, M, 3 - 1, I];
						}
					}
				}
			}

			if ((ITD != I2D) & (ITD != I3D))
			{
				Console.WriteLine("Incompatible potential types!!!");
				Console.ReadLine();
			}

			PQ = Math.PI / (XQM - ATT);

			ADB = new double[4];
			CDB = new double[4];
			CDB2 = new double[4];
			DDB = new double[4];
			DDB2 = new double[4];
			HDB = new double[4];
			XTN1 = new double[4];
			XTN2 = new double[4];
			AD = new double[4, 4];
			AXL = new double[4, 4];
			BD = new double[4, 4];
			BXL = new double[4, 4];
			CD = new double[4, 4];
			CXL = new double[4, 4];
			DD = new double[4, 4];
			DXL = new double[4, 4];
			ED = new double[4, 4];
			RB1 = new double[4, 4];
			RB2 = new double[4, 4];
			RMAX = new double[4, 4];
			PID = new double[4, 4];
			CHI = new double[4, 4];
			XDB = new double[4, 4, 4];
			REG = new double[4, 4, 4];

			Array.Clear(ADB, 0, 4);
			Array.Clear(CDB, 0, 4);
			Array.Clear(CDB2, 0, 4);
			Array.Clear(DDB, 0, 4);
			Array.Clear(DDB2, 0, 4);
			Array.Clear(HDB, 0, 4);
			Array.Clear(XTN1, 0, 4);
			Array.Clear(XTN2, 0, 4);
			Array.Clear(AD, 0, 16);
			Array.Clear(AXL, 0, 16);
			Array.Clear(BD, 0, 16);
			Array.Clear(BXL, 0, 16);
			Array.Clear(CD, 0, 16);
			Array.Clear(CXL, 0, 16);
			Array.Clear(DD, 0, 16);
			Array.Clear(DXL, 0, 16);
			Array.Clear(ED, 0, 16);
			Array.Clear(RB1, 0, 16);
			Array.Clear(RB2, 0, 16);
			Array.Clear(RMAX, 0, 16);
			Array.Clear(PID, 0, 16);
			Array.Clear(CHI, 0, 16);
			Array.Clear(XDB, 0, 4*4*4);
			Array.Clear(REG, 0, 4 * 4 * 4);

			//******IMPORTANT*******
			// TO INCLUDE DIHEDRAL TERMS
			// SET NDIHED=2, OTHERWISE
			// SET NDIHED=10
			NDIHED = 2;

			// carbon
			AD[0,0]=12388.79197798375;
			AXL[0,0]=4.720452312717397;
			BD[0,0]=17.56740646508968;
			BXL[0,0]=1.433213249951261;
			CD[0,0]=30.71493208065162;
			CXL[0,0]=1.382691250599169;
			DD[0,0]=10953.54416216992;
			DXL[0,0]=4.746539060659529;
			ED[0,0]=0.3134602960832605;
			RB1[0,0]=1.7;
			RB2[0,0]=2.0;
			RMAX[0, 0] = RB2[0, 0];
			PID[0,0]=Math.PI/(RB2[0,0]-RB1[0,0]);

			// hydrogen
			AD[1,1]=29.6325931;
			AXL[1, 1] = 1.715892169856421;
			BD[1, 1] = 0.0;
			BXL[1, 1] = 1.0;
			CD[1, 1] = 0.0;
			CXL[1, 1] = 1.0;
			DD[1, 1] = 32.81735574722296;
			DXL[1, 1] = 3.536298648376465;
			ED[1, 1] = 0.3704714870452888;
			RB1[1, 1] = 1.10;
			RB2[1, 1] = 1.70;
			RMAX[1, 1] = RB2[1, 1];
			PID[1, 1] = Math.PI / (RB2[1, 1] - RB1[1, 1]);

			//carbon-hydrogen
			AD[1,0]=32.35518665873256;
			AXL[1,0]=1.434458059249837;
			DD[1,0]=149.9409872288120;
			DXL[1,0]= 4.102549828548784;
			ED[1,0]=0.3407757282257080;
			BD[1,0]=0.0;
			BXL[1,0]=1.0;
			CD[1,0]=0.0;
			CXL[1,0]=1.0;
			AD[0,1]=AD[1,0];
			AXL[0,1]=AXL[1,0];
			BD[0,1]=BD[1,0];
			BXL[0,1]=BXL[1,0];
			CD[0,1]=CD[1,0];
			CXL[0,1]=CXL[1,0];
			DD[0,1]=DD[1,0];
			DXL[0,1]=DXL[1,0];
			ED[0,1]=ED[1,0];
			RB1[1,0]=1.3;
			RB2[1,0]=1.8;
			RMAX[1,0]=RB2[1,0];
			PID[1,0]=Math.PI/(RB2[1,0]-RB1[1,0]);
			PIDT=Math.PI/0.30;
			RB1[0,1]=RB1[1,0];
			RB2[0,1]=RB2[1,0];
			RMAX[0,1]=RB2[0,1];
			PID[0,1]=Math.PI/(RB2[0,1]-RB1[0,1]);

			for (int I=0; I<2; I++)
			{
				for (int J=0; J<2; J++)
				{
					for (int K=0; K<2; K++)
					{
						XDB[I,J,K]=0.0;
						REG[I,J,K]=1.0;
					}
				}
			}

			XXDB=4.0;
			double RHH=0.7415886997;
			double RCH=1.09;
			XDB[1, 1, 1]=4.0;
			XDB[1, 0, 1] = 4.0;
			XDB[1, 1, 0] = 4.0;
			XDB[1, 0, 0] = 4.0;
			XDB[0, 1, 0] = 0.0;
			XDB[0, 1, 1] = 0.0;
			REG[1, 0, 1] = Math.Exp(XDB[1, 0, 1] * (RHH - RCH));
			REG[1, 1, 0] = Math.Exp(XDB[1, 1, 0] * (RCH - RHH));

			//Tersoff-III silicon
			DXL[2,2]=2.4799;
			AXL[2,2]=1.7322;
			DD[2,2]=1830.8;
			AD[2,2]=471.18;
			XTN2[2]=0.78734;
			XTN1[2] = 1.0 / (2.0 * XTN2[2]);
			ADB[2]=1.0999E-6;
			CDB[2]=1.0039E+5;
			CDB2[2] = CDB[2]*CDB[2];
			DDB[2]=16.218;
			DDB2[2] = DDB[2]*DDB[2];
			HDB[2]=0.59826;
			RB1[2,2]=2.7;
			RB2[2,2]=3.0;
			RMAX[2,2]=RB2[2,2];
			PID[2,2]=Math.PI/(RB2[2,2]-RB1[2,2]);
			CHI[2,2] = 1.0;

			//Tersoff germanium
			DXL[3,3] = 2.4451E0;
			AXL[3, 3] = 1.7047E0;
			DD[3, 3] = 1769.0E0;
			AD[3, 3] = 0.50E0 * 419.23E0;
			XTN2[3] = 0.75627E0;
			XTN1[3] = 1 / (2.0 * XTN2[3]);
			ADB[3] = 9.0166E-07;
			CDB[3] = 1.0643E+5;
			CDB2[3] = CDB[3] * CDB[3];
			DDB[3] = 15.652E0;
			DDB2[3] = DDB[3] * DDB[3];
			HDB[3] = -0.43884E0;
			RB1[3, 3] = 2.7E0;
			RB2[3, 3] = 3.0E0;
			RMAX[3, 3] = RB2[3, 3];
			PID[3, 3] = Math.PI / (RB2[3, 3] - RB1[3, 3]);
			CHI[3, 3] = 1.0E0;

			// silicon-germanium
			DXL[3,2] = (DXL[3,3]+DXL[2,2])/2.0E0;
			DXL[2,3] = DXL[3,2];
			AXL[3,2] = (AXL[3,3]+AXL[2,2])/2.0E0;
			AXL[2,3] = AXL[3,2];
			DD[3,2] = Math.Sqrt(DD[3,3]*DD[2,2]);
			DD[2,3] = DD[3,2];
			AD[3,2] = Math.Sqrt(AD[2,2]*AD[3,3]);
			AD[2,3] = AD[3,2];
			RB1[3,2] = Math.Sqrt(RB1[2,2]*RB1[3,3]);
			RB1[2,3] = RB1[3,2];
			RB2[3,2] = Math.Sqrt(RB2[2,2]*RB2[3,3]);
			RB2[2,3] = RB2[3,2];
			RMAX[3,2] = RB2[3,2];
			RMAX[2,3] = RMAX[3,2];
			PID[3,2] = Math.PI/(RB2[3,2]-RB1[3,2]);
			PID[2,3] = PID[3,2];
			CHI[3,2] = 1.00061E0;
			CHI[2,3] = CHI[3,2];
			SIGMA=1.0E0;
			EPSI=11605.0E0;

			RLIST = new double[4, 4];
			Array.Clear(RLIST, 0, 16);

			for (int I=0; I<4; I++)
			{
				for (int J=0; J<4; J++)
				{
					RLIST[I,J]=Math.Pow((RMAX[I,J]+RLL),2.0); // this is the cut-off radius^2 for creating neighbor list.
					RMAX[I,J]=Math.Pow(RMAX[I,J],2.0);
				}
			}
        }

		private void mtable() //Creates tables of the potentials VA and VR from the Brenner paper
		{
			DDTAB = new double[4, 4];
			Array.Clear(DDTAB, 0, 16);
			TABFC = new double[4,4,NTAB];
			Array.Clear(TABFC,0,4*4*NTAB);
			TABDFC = new double[4,4,NTAB];
			Array.Clear(TABDFC,0,4*4*NTAB);
			ATABLE = new double[4,4,NTAB];
			Array.Clear(ATABLE,0,4*4*NTAB);
			DATABLE = new double[4,4,NTAB];
			Array.Clear(DATABLE,0,4*4*NTAB);
			RTABLE = new double[4,4,NTAB];
			Array.Clear(RTABLE,0,4*4*NTAB);
			DRTABLE = new double[4,4,NTAB];
			Array.Clear(DRTABLE,0,4*4*NTAB);

			// generate lookup tables for bond-order potentials
			for (int ki=0; ki<4; ki++)
			{
				for (int kj=0; kj<4; kj++)
				{
					DDTAB[ki, kj] = RB2[ki, kj] / ((double)(NTAB - 2));
					DDTAB[kj, ki] = DDTAB[ki, kj];
					double rc = 0.0;
					for (int i=1; i<NTAB-1; i++)
					{
						if (DDTAB[ki,kj] != 0.0)
						{
							rc = rc + DDTAB[ki, kj];
							double rsq = rc * rc;
							// cut-off function
							double FC = 0.0;
							double DFC = 0.0;

							if (rc < RB2[ki,kj])
							{
								double DTEMP = PID[ki, kj] * (rc - RB1[ki, kj]);
								FC = (1.0 + Math.Cos(DTEMP)) / 2.0;
								DFC = -PID[ki, kj] / 2.0 * Math.Sin(DTEMP);
							}

							if (rc <= RB1[ki,kj])
							{
								FC = 1.0;
								DFC = 0.0;
							}

							TABFC[ki,kj,i] = FC;
							TABFC[kj,ki,i] = TABFC[ki,kj,i];
							TABDFC[ki,kj,i] = DFC;
							TABDFC[kj,ki,i] = TABDFC[ki,kj,i];
							// attractive pair terms
							double VA=AD[ki,kj]*Math.Exp(-AXL[ki,kj]*rc);
							double DVA=-AXL[ki,kj]*VA;

							double VB=BD[ki,kj]*Math.Exp(-BXL[ki,kj]*rc);
							double DVB=-BXL[ki,kj]*VB;

							double VC=CD[ki,kj]*Math.Exp(-CXL[ki,kj]*rc);
							double DVC=-CXL[ki,kj]*VC;

							double VV=(VA+VB+VC)/2.0;
							double DVV=(DVA+DVB+DVC)/2.0;
							ATABLE[ki,kj,i] = FC*VV;
							ATABLE[kj,ki,i] = ATABLE[ki,kj,i];
							DATABLE[ki,kj,i] = (FC*DVV+DFC*VV)/rc;
							DATABLE[kj,ki,i] = DATABLE[ki,kj,i];
							// repulsive pair terms
							double FF1=DD[ki,kj]*Math.Exp(-DXL[ki,kj]*rc);
							double DF1=-DXL[ki,kj]*FF1;

							double FF2=(1.0+ED[ki,kj]/rc);
							double DF2=-ED[ki,kj]/rsq;

							VV = FF1 * FF2;
							double DVM = (DF1 * FF2 + FF1 * DF2);
							RTABLE[ki, kj, i] = VV * FC;
							RTABLE[kj, ki, i] = RTABLE[ki, kj, i];
							DRTABLE[ki, kj, i] = -(FC * DVM + DFC * VV) / rc;
							DRTABLE[kj, ki, i] = DRTABLE[ki, kj, i];
						}
						else
						{
							TABFC[ki,kj,i]=0.0;
							TABFC[kj,ki,i] = TABFC[ki,kj,i];
							TABDFC[ki,kj,i] = 0.0;
							TABDFC[kj,ki,i] = TABDFC[ki,kj,i];
							ATABLE[ki,kj,i] = 0.0;
							ATABLE[kj,ki,i] = ATABLE[ki,kj,i];
							DATABLE[ki,kj,i] = 0.0;
							DATABLE[kj,ki,i] = DATABLE[ki,kj,i];
							RTABLE[ki,kj,i] = 0.0;
							RTABLE[kj,ki,i] = RTABLE[ki,kj,i];
							DRTABLE[ki,kj,i] = 0.0;
							DRTABLE[kj,ki,i] = DRTABLE[ki,kj,i];
						}
					}
					ATABLE[ki, kj, 0] = ATABLE[ki, kj, 1];
					ATABLE[kj, ki, 0] = ATABLE[ki, kj, 0];
					DATABLE[ki, kj, 0] = DATABLE[ki, kj, 1];
					DATABLE[kj, ki, 0] = DATABLE[ki, kj, 0];
					RTABLE[ki, kj, 0] = RTABLE[ki, kj, 1];
					RTABLE[kj, ki, 0] = RTABLE[ki, kj, 0];
					DRTABLE[ki, kj, 0] = DRTABLE[ki, kj, 1];
					DRTABLE[kj, ki, 0] = DRTABLE[ki, kj, 0];
					TABFC[ki, kj, 0] = 0.0;
					TABFC[kj, ki, 0] = 0.0;
					TABDFC[ki, kj, 0] = 0.0;
					TABDFC[kj, ki, 0] = 0.0;

					ATABLE[ki, kj, NTAB - 1] = 0.0;
					ATABLE[kj, ki, NTAB - 1] = ATABLE[ki, kj, NTAB - 1];
					DATABLE[ki, kj, NTAB - 1] = 0.0;
					DATABLE[kj, ki, NTAB - 1] = DATABLE[ki, kj, NTAB - 1];
					RTABLE[ki, kj, NTAB - 1] = 0.0;
					RTABLE[kj, ki, NTAB - 1] = RTABLE[ki, kj, NTAB - 1];
					DRTABLE[ki, kj, NTAB - 1] = 0.0;
					DRTABLE[kj, ki, NTAB - 1] = DRTABLE[ki, kj, NTAB - 1];
					TABFC[ki, kj, NTAB - 1] = 0.0;
					TABFC[kj, ki, NTAB - 1] = TABFC[ki, kj, NTAB - 1];
					TABDFC[ki, kj, NTAB - 1] = 0.0;
					TABDFC[kj, ki, NTAB - 1] = TABDFC[ki, kj, NTAB - 1];
				}
			}
		}

		public void caguts() // Finds forces on atoms from pair potentials.
		{
			// calculate two-body forces and neighbor list for hydrocarbons
			double[] RR = new double[3];
			double[] RI = new double[3];
			NABORS = new int[Atom.AllAtoms.Count + 1];

			RNP = new double[Atom.AllAtoms.Count, 3];
			Array.Clear(RNP, 0, Atom.AllAtoms.Count * 3);

			{
				bool flg = checkNaN(RNP);
				if (flg)
				{
					Debugger.Break();
				}
			}

			eatom = new double[Atom.AllAtoms.Count];
			Array.Clear(eatom, 0, eatom.Length);
			int KEND = 0; // this is the total number of neighbor pairs that are calculated and stored in LIST, IVCT2B, and JVCT2B
			if (neighborListFlag)
			{
				// set up neighbor list
				int K = 0;
				IVCT2B.Clear();
				JVCT2B.Clear();
				for (int I = 0; I < Atom.AllAtoms.Count; I++)
				{
					NABORS[I] = K;
					RI[0] = Atom.AllAtoms[I].position[0];
					RI[1] = Atom.AllAtoms[I].position[1];
					RI[2] = Atom.AllAtoms[I].position[2];
					int ki = 5;
					if (Atom.AllAtoms[I].atomicNumber == 6)
					{
						ki = 0;
					}
					else if (Atom.AllAtoms[I].atomicNumber == 1)
					{
						ki = 1;
					}
					else
					{
						Console.WriteLine("Atom types should be either carbon or hydrogen to use REBO potential!");
						Console.ReadLine();
						return;
					}

					for (int J = 0; J < Atom.AllAtoms.Count; J++)
					{
						if (I != J)
						{
							int kj = 5;
							if (Atom.AllAtoms[J].atomicNumber == 6)
							{
								kj = 0;
							}
							else if (Atom.AllAtoms[J].atomicNumber == 1)
							{
								kj = 1;
							}
							else
							{
								Console.WriteLine("Atom types should be either carbon or hydrogen to use REBO potential!");
								Console.ReadLine();
								return;
							}
							double RLIS = RLIST[ki, kj];
							double rsq = 0.0;
							for (int L = 0; L < 3; L++)
							{
								RR[L] = RI[L] - Atom.AllAtoms[J].position[L];
								rsq = rsq + RR[L] * RR[L];
							}

							if (rsq <= RLIS)
							{
								K++;
								IVCT2B.Add(I);
								JVCT2B.Add(J);
							}


						}
					}
				}
				KEND = K; // KEND is number of interacting atom pairs
				NABORS[Atom.AllAtoms.Count] = K;
			}

			//debug
			if (KEND != IVCT2B.Count)
			{
				Console.WriteLine("There is a problem with index number of IVCT2B versus KEND! Please investigate!");
				Console.ReadLine();
			}
			//debug


			LCHECK = new int[KEND];
			Array.Clear(LCHECK, 0, KEND);
			COR = new double[KEND, 3];
			Array.Clear(COR, 0, 3 * KEND);
			RCOR = new double[KEND];
			Array.Clear(RCOR, 0, KEND);
			WW = new double[KEND];
			Array.Clear(WW, 0, KEND);
			DWW = new double[KEND];
			Array.Clear(DWW, 0, KEND);
			EXX1 = new double[KEND];
			Array.Clear(EXX1, 0, KEND);
			DEXX1 = new double[KEND];
			Array.Clear(DEXX1, 0, KEND);
			double[,] RPP = new double[KEND, 3];
			Array.Clear(RPP, 0, 3 * KEND);

			for (int K = 0; K < KEND; K++)
			{
				int I = IVCT2B[K];
				int J = JVCT2B[K];
				int KI = 5;
				if (Atom.AllAtoms[I].atomicNumber == 6)
				{
					KI = 0;
				}
				else if (Atom.AllAtoms[I].atomicNumber == 1)
				{
					KI = 1;
				}
				else
				{
					Console.WriteLine("Atom types should be either carbon or hydrogen to use REBO potential!");
					Console.ReadLine();
					return;
				}
				int KJ = 5;
				if (Atom.AllAtoms[J].atomicNumber == 6)
				{
					KJ = 0;
				}
				else if (Atom.AllAtoms[J].atomicNumber == 1)
				{
					KJ = 1;
				}
				else
				{
					Console.WriteLine("Atom types should be either carbon or hydrogen to use REBO potential!");
					Console.ReadLine();
					return;
				}

				LCHECK[K] = 0;
				double rsq = 0.0;
				for (int L = 0; L < 3; L++)
				{
					RR[L] = Atom.AllAtoms[I].position[L] - Atom.AllAtoms[J].position[L];
					rsq = rsq + RR[L] * RR[L];
					COR[K, L] = RR[L];
				}

				if (rsq <= RMAX[KI, KJ])
				{

					if ((KJ <= 1) && (KI <= 1)) LCHECK[K] = 1;
					if ((KJ >= 2) && (KI >= 2)) LCHECK[K] = 2;

					double rc = Math.Sqrt(rsq);
					double rt = rc / DDTAB[KI, KJ];
					int it = Math.Min((int)rt, NTAB - 2);

					RCOR[K] = rc;
					WW[K] = TABFC[KI, KJ, it] + (TABFC[KI, KJ, it + 1] - TABFC[KI, KJ, it]) * (rt - (double)it + 1.0);
					DWW[K] = TABDFC[KI, KJ, it] + (TABDFC[KI, KJ, it + 1] - TABDFC[KI, KJ, it]) * (rt - (double)it + 1.0);
					EXX1[K] = ATABLE[KI, KJ, it] + (ATABLE[KI, KJ, it + 1] - ATABLE[KI, KJ, it]) * (rt - (double)it + 1.0);
					DEXX1[K] = DATABLE[KI, KJ, it] + (DATABLE[KI, KJ, it + 1] - DATABLE[KI, KJ, it]) * (rt - (double)it + 1.0);

					if (I < J)
					{
						double vv = RTABLE[KI, KJ, it] + (RTABLE[KI, KJ, it + 1] - RTABLE[KI, KJ, it]) * (rt - (double)it + 1.0); // this is the potential energy between atom pair I and J
						double rp = DRTABLE[KI, KJ, it] + (DRTABLE[KI, KJ, it + 1] - DRTABLE[KI, KJ, it]) * (rt - (double)it + 1.0);
						tote = tote + vv;
						eatom[I] = eatom[I] + vv / 2.0;
						eatom[J] = eatom[J] + vv / 2.0;

						for (int L = 0; L < 3; L++)
						{
							RPP[K, L] = rp * RR[L];
						}
					}
				}
			}

			for (int K = 0; K < KEND; K++)
			{
				if (LCHECK[K] != 0)
				{
					int I = IVCT2B[K];
					int J = JVCT2B[K];
					if (I < J)
					{
						for (int L = 0; L < 3; L++)
						{
							RNP[I, L] = RNP[I, L] + RPP[K, L];
							RNP[J, L] = RNP[J, L] + RPP[K, L];
						}
					}
				}
			}

			{
				bool flg = checkNaN(RNP);
				if (flg)
				{
					Debugger.Break();
				}
			}

			// add a check to see if the atoms are carbohydrates or si-germanium and perform the corresponding routines
			// call pibond for carbohydrates
			// call sili_germ for silicon germanium atoms.

			for (int i = 0; i < Atom.AllAtoms.Count; i++)
			{
				if ((Atom.AllAtoms[i].atomicNumber != 1) && (Atom.AllAtoms[i].atomicNumber != 6))
				{
					Console.WriteLine("Atoms are not hydrogen or carbon! Correct potential type should be used!!!!");
					Console.ReadLine();
				}
			}
		}

		private void pibond() // This function calculates the forces due to pi-bonds between atoms
		{
			double[,] XHC = new double[Atom.AllAtoms.Count,2];



			// Find number of hydrogens and carbons connected to each atom
			for (int i = 0; i < Atom.AllAtoms.Count; i++ )
			{
				int JBEGIN = NABORS[i];
				int JEND = NABORS[i + 1];

				XHC[i, 0] = 1;
				XHC[i, 1] = 1;
				for (int j = JBEGIN; j < JEND; j++ )
				{
					int JN = JVCT2B[j];
					int kj = 5;

					if (Atom.AllAtoms[JN].atomicNumber == 6)
					{
						kj = 0;
					}
					else if (Atom.AllAtoms[JN].atomicNumber == 1)
					{
						kj = 1;
					}
					else
					{
						Console.WriteLine("Atom types should be either carbon or hydrogen to use REBO potential!");
						Console.ReadLine();
						return;
					}

					XHC[i,kj] = XHC[i,kj] + WW[j];
					

				}

			}

			// Sum over bonds between atoms I and J
			for (int i = 0; i < Atom.AllAtoms.Count; i++ )
			{
				int JBEGIN = NABORS[i];
				int JEND = NABORS[i + 1];

				int ki = 5;
				if (Atom.AllAtoms[i].atomicNumber == 6)
				{
					ki = 0;
				}
				else if (Atom.AllAtoms[i].atomicNumber == 1)
				{
					ki = 1;
				}
				else
				{
					Console.WriteLine("Atom types should be either carbon or hydrogen to use REBO potential!");
					Console.ReadLine();
					return;
				}

				for (int j=JBEGIN; j<(JEND); j++)
				{
					int jn = JVCT2B[j];
					if (i < jn)
					{
						double[] CJ = new double[3];
						for (int mm = 0; mm < 3; mm++)
						{
							CJ[mm] = COR[j, mm];
						}
						double SIJ = RCOR[j];
						double RSQIJ = SIJ * SIJ;

						int kj = 5;
						if (Atom.AllAtoms[jn].atomicNumber == 6)
						{
							kj = 0;
						}
						else if (Atom.AllAtoms[jn].atomicNumber == 1)
						{
							kj = 1;
						}
						else
						{
							Console.WriteLine("Atom types should be either carbon or hydrogen to use REBO potential!");
							Console.ReadLine();
							return;
						}

						int kikj = ki + kj;

						// I side of bond
						int nk = 0;
						List<double[]> xk = new List<double[]>();
						List<double> cosk = new List<double>();
						List<double> sink = new List<double>();
						List<double> cfuni = new List<double>();
						List<double> dcfuni = new List<double>();
						List<double> dctjk = new List<double>();
						List<double> dctij = new List<double>();
						List<double> dctik = new List<double>();
						List<double> xsik = new List<double>();
						List<double> xsjk = new List<double>();


						double xsij = 0;
						double ssumk = 0;
						double conk = 0;
						double[] xni = new double[2] { XHC[i, 0], XHC[i, 1] };
						xni[kj] = xni[kj] - WW[j];
						double qi = xni[0] + xni[1] - 2.0e0;
						double sdalik = 0;

						for (int k = JBEGIN; k < JEND; k++)
						{
							double ali = 0;
							double dali = 0;
							double daldik = 0;
							double gangle = 0;
							double dgdthet = 0;

							if ((k != j))
							{
								int KN = JVCT2B[k];

								int kk = 5;
								if (Atom.AllAtoms[KN].atomicNumber == 6)
								{
									kk = 0;
								}
								else if (Atom.AllAtoms[KN].atomicNumber == 1)
								{
									kk = 1;
								}
								else
								{
									Console.WriteLine("Atom types should be either carbon or hydrogen to use REBO potential!");
									Console.ReadLine();
									return;
								}

								nk = nk + 1;
								double s3 = RCOR[k];
								double rsq3 = s3 * s3;
								double rsq2 = 0;

								double[] tempArray = new double[3];
								for (int mm = 0; mm < 3; mm++)
								{
									tempArray[mm] = COR[k, mm] - CJ[mm];
									rsq2 = rsq2 + tempArray[mm] * tempArray[mm];
								}
								xk.Add(tempArray);

								double ss = 2.0 * SIJ * s3;
								double rr = RSQIJ - rsq3;
								double costh = (RSQIJ + rsq3 - rsq2) / ss;
								if (costh > 1.0) costh = 1.0;
								if (costh < -1.0) costh = -1.0;
								cosk.Add(costh);
								sink.Add(Math.Sqrt(1.0 - costh * costh));
								if (Math.Acos(costh) > Math.PI) sink[nk - 1] = -sink[nk - 1];


								if (ki == 0)
								{
									int ig = IGC[(int)(-costh * 12.0) + 13 - 1] - 1;
									if (ig != 3)
									{
										gangle = SPGC[0, ig] + SPGC[1, ig] * costh;
										dgdthet = SPGC[1, ig];
										for (int jj = 2; jj < 6; jj++)
										{
											gangle = gangle + SPGC[jj, ig] * (Math.Pow(costh, jj));
											dgdthet = dgdthet + SPGC[jj, ig] * ((double)jj) * Math.Pow(costh, jj - 1);
										}
									}
									else
									{
										ali = 0;
										dali = 0;
										if (qi < XQM)
										{
											ali = 1.0;
											if (qi > ATT)
											{
												double dtemp = PQ * (qi - ATT);
												ali = (1.0 + Math.Cos(dtemp)) / 2.0;
												dali = -PQ / 2.0 * Math.Sin(dtemp);
											}
										}
										gangle = SPGC[0, ig] + SPGC[1, ig] * costh;
										dgdthet = SPGC[1, ig];
										int ig1 = ig + 1;
										double gangle1 = SPGC[0, ig1] + SPGC[1, ig1] * costh;
										double dgdthet1 = SPGC[1, ig1];
										for (int jj = 2; jj < 6; jj++)
										{
											gangle = gangle + SPGC[jj, ig] * Math.Pow(costh, jj);
											dgdthet = dgdthet + SPGC[jj, ig] * ((double)jj) * Math.Pow(costh, jj - 1);
											gangle1 = gangle1 + SPGC[jj, ig1] * Math.Pow(costh, jj);
											dgdthet1 = dgdthet1 + SPGC[jj, ig1] * ((double)jj) * Math.Pow(costh, jj - 1);
										}
										daldik = dali * (gangle1 - gangle);
										gangle = gangle + ali * (gangle1 - gangle);
										dgdthet = dgdthet + ali * (dgdthet1 - dgdthet);
									}
								}
								else
								{
									int ig = IGH[(int)(-costh * 12.0) + 13 - 1] - 1;
									gangle = SPGH[0, ig] + SPGH[1, ig] * costh;
									dgdthet = SPGH[1, ig];
									for (int jj = 2; jj < 6; jj++)
									{
										gangle = gangle + SPGH[jj, ig] * Math.Pow(costh, jj);
										dgdthet = dgdthet + SPGH[jj, ig] * (double)(jj) * Math.Pow(costh, jj - 1);
									}
								}

								double fc = WW[k];
								double dfc = DWW[k];
								cfuni.Add(0);
								dcfuni.Add(0);

								if (kk == 0)
								{
									double xx = XHC[KN, 0] + XHC[KN, 1] - fc - 2.0;
									if (xx < 3.0)
									{
										if (xx <= 2.0)
										{
											cfuni[nk - 1] = 1.0;
										}
										else
										{
											double px = Math.PI * (xx - 2.0);
											cfuni[nk - 1] = (1.0 + Math.Cos(px)) / 2.0;
											dcfuni[nk - 1] = -fc * Math.Sin(px) * Math.PI / 2.0;
										}
									}
								}
								conk = conk + fc * cfuni[nk - 1];

								double exx = 0;
								if (XDB[ki, kj, kk] != 0)
								{
									exx = REG[ki, kj, kk] * Math.Exp(XDB[ki, kj, kk] * (SIJ - s3));
								}
								else
								{
									exx = 1.0;
								}

								double dctdjk = -2.0 / ss;
								double dctdij = (rr + rsq2) / (ss * RSQIJ);
								double dctdik = (-rr + rsq2) / (ss * rsq3);
								dctjk.Add(dctdjk);
								dctij.Add(dctdij);
								dctik.Add(dctdik);
								double gs = gangle * exx;
								ssumk = ssumk + fc * gs;
								double xtemp = fc * exx * dgdthet;
								double gfx = gs * fc * XDB[ki, kj, kk];
								xsij = xsij + xtemp * dctdij + gfx / SIJ;
								xsik.Add((gs * dfc - gfx) / s3 + xtemp * dctdik);
								sdalik = sdalik + exx * fc * daldik;
								xsjk.Add(xtemp * dctdjk);
							}
						}

						// J side of bond
						int nl = 0;
						List<double[]> xl = new List<double[]>();
						List<double> cosl = new List<double>();
						List<double> sinl = new List<double>();
						List<double> cfunj = new List<double>();
						List<double> dcfunj = new List<double>();

						List<double> dctil = new List<double>();
						List<double> dctji = new List<double>();
						List<double> dctjl = new List<double>();
						List<double> xsjl = new List<double>();
						List<double> xsil = new List<double>();


						double xsji = 0;
						double ssuml = 0;
						double conl = 0;
						int LBEGIN = NABORS[jn];
						int LEND = NABORS[jn + 1];
						double[] xnj = new double[2] { XHC[jn, 0], XHC[jn, 1] };
						xnj[ki] = xnj[ki] - WW[j];
						double qj = xnj[0] + xnj[1] - 2.0e0;
						double sdaljl = 0;

						for (int l = LBEGIN; l < LEND; l++)
						{
							double alj = 0;
							double dalj = 0;
							double daldjl = 0;
							int ln = JVCT2B[l];
							double gangle = 0;
							double dgdthet = 0;

							if ((ln != i))
							{

								int kl = 5;
								if (Atom.AllAtoms[ln].atomicNumber == 6)
								{
									kl = 0;
								}
								else if (Atom.AllAtoms[ln].atomicNumber == 1)
								{
									kl = 1;
								}
								else
								{
									Console.WriteLine("Atom types should be either carbon or hydrogen to use REBO potential!");
									Console.ReadLine();
									return;
								}

								nl = nl + 1;
								double s3 = RCOR[l];
								double rsq3 = s3 * s3;
								double rsq2 = 0;

								double[] tempArray = new double[3];
								for (int mm = 0; mm < 3; mm++)
								{
									tempArray[mm] = COR[l, mm] - CJ[mm];
									rsq2 = rsq2 + tempArray[mm] * tempArray[mm];
								}
								xl.Add(tempArray);

								double ss = 2.0 * SIJ * s3;
								double rr = RSQIJ - rsq3;
								double costh = (RSQIJ + rsq3 - rsq2) / ss;
								if (costh > 1.0) costh = 1.0;
								if (costh < -1.0) costh = -1.0;
								cosl.Add(costh);
								sinl.Add(Math.Sqrt(1.0 - costh * costh));
								if (Math.Acos(costh) > Math.PI) sinl[nl - 1] = -sinl[nl - 1];

								if (kj == 0)
								{
									int ig = IGC[(int)(-costh * 12.0) + 13 - 1] - 1;
									if (ig != 3)
									{
										gangle = SPGC[0, ig] + SPGC[1, ig] * costh;
										dgdthet = SPGC[1, ig];
										for (int jj = 2; jj < 6; jj++)
										{
											gangle = gangle + SPGC[jj, ig] * (Math.Pow(costh, jj));
											dgdthet = dgdthet + SPGC[jj, ig] * ((double)jj) * Math.Pow(costh, jj - 1);
										}
									}
									else
									{
										alj = 0;
										dalj = 0;
										if (qj < XQM)
										{
											alj = 1.0;
											if (qj > ATT)
											{
												double dtemp = PQ * (qj - ATT);
												alj = (1.0 + Math.Cos(dtemp)) / 2.0;
												dalj = -PQ / 2.0 * Math.Sin(dtemp);
											}
										}
										gangle = SPGC[0, ig] + SPGC[1, ig] * costh;
										dgdthet = SPGC[1, ig];
										int ig1 = ig + 1;
										double gangle1 = SPGC[0, ig1] + SPGC[1, ig1] * costh;
										double dgdthet1 = SPGC[1, ig1];
										for (int jj = 2; jj < 6; jj++)
										{
											gangle = gangle + SPGC[jj, ig] * Math.Pow(costh, jj);
											dgdthet = dgdthet + SPGC[jj, ig] * ((double)jj) * Math.Pow(costh, jj - 1);
											gangle1 = gangle1 + SPGC[jj, ig1] * Math.Pow(costh, jj);
											dgdthet1 = dgdthet1 + SPGC[jj, ig1] * ((double)jj) * Math.Pow(costh, jj - 1);
										}
										daldjl = dalj * (gangle1 - gangle);
										gangle = gangle + alj * (gangle1 - gangle);
										dgdthet = dgdthet + alj * (dgdthet1 - dgdthet);
									}
								}
								else
								{
									int ig = IGH[(int)(-costh * 12.0) + 13 - 1] - 1;
									gangle = SPGH[0, ig] + SPGH[1, ig] * costh;
									dgdthet = SPGH[1, ig];
									for (int jj = 2; jj < 6; jj++)
									{
										gangle = gangle + SPGH[jj, ig] * Math.Pow(costh, jj);
										dgdthet = dgdthet + SPGH[jj, ig] * (double)(jj) * Math.Pow(costh, jj - 1);
									}
								}

								double fc = WW[l];
								double dfc = DWW[l];
								cfunj.Add(0);
								dcfunj.Add(0);

								if (kl == 0)
								{
									double xx = XHC[ln, 0] + XHC[ln, 1] - fc - 2.0;
									if (xx < 3.0)
									{
										if (xx <= 2.0)
										{
											cfunj[nl - 1] = 1.0;
										}
										else
										{
											double px = Math.PI * (xx - 2.0);
											cfunj[nl - 1] = (1.0 + Math.Cos(px)) / 2.0;
											dcfunj[nl - 1] = -fc * Math.Sin(px) * Math.PI / 2.0;
										}
									}
								}
								conl = conl + fc * cfunj[nl - 1];

								double exx = 0;
								if (XDB[kj, ki, kl] != 0)
								{
									exx = REG[kj, ki, kl] * Math.Exp(XDB[kj, ki, kl] * (SIJ - s3));
								}
								else
								{
									exx = 1.0;
								}

								double dctdil = -2.0 / ss;
								double dctdji = (rr + rsq2) / (ss * RSQIJ);
								double dctdjl = (-rr + rsq2) / (ss * rsq3);
								dctil.Add(dctdil);
								dctji.Add(dctdji);
								dctjl.Add(dctdjl);
								double gs = gangle * exx;
								ssuml = ssuml + fc * gs;
								double xtemp = fc * exx * dgdthet;
								double gfx = gs * fc * XDB[kj, ki, kl];
								xsji = xsji + xtemp * dctdji + gfx / SIJ;
								xsjl.Add((gs * dfc - gfx) / s3 + xtemp * dctdjl);
								sdaljl = sdaljl + exx * fc * daldjl;
								xsil.Add(xtemp * dctdil);
							}
						}

						double exnij = 0;
						double[] dexni = new double[2] { 0, 0 };

						if (ki == 0)
						{
							int nh = (int)(xni[1] + Math.Pow(1, -12));
							int nc = (int)(xni[0] + Math.Pow(1, -12));
							if ((Math.Abs((double)nh - xni[1]) > Math.Pow(1, -8)) || (Math.Abs((double)nc - xni[0]) > Math.Pow(1, -8)))
							{
								bcuint(ki, kj, xni[1], xni[0], nh - 1, nc - 1, exnij, dexni[1], dexni[0]);
							}
							else
							{
								exnij = XH[kj, nh - 1, nc - 1];
								dexni[1] = XH1[kj, nh - 1, nc - 1];
								dexni[0] = XH2[kj, nh - 1, nc - 1];
							}

						}

						double exnji = 0;
						double[] dexnj = new double[2] { 0, 0 };
						if (kj == 0)
						{
							int nh = (int)(xnj[1] + Math.Pow(1, -12));
							int nc = (int)(xnj[0] + Math.Pow(1, -12));
							if ((Math.Abs((double)nh - xnj[1]) > Math.Pow(1, -8)) || (Math.Abs((double)nc - xnj[0]) > Math.Pow(1, -8)))
							{
								bcuint(kj, ki, xnj[1], xnj[0], nh - 1, nc - 1, exnji, dexnj[1], dexnj[0]);
							}
							else
							{
								exnji = XH[ki, nh - 1, nc - 1];
								dexnj[1] = XH1[ki, nh - 1, nc - 1];
								dexnj[0] = XH2[ki, nh - 1, nc - 1];
							}
						}

						double dij = 1.0 + exnij + ssumk;
						double bij = 1.0 / Math.Sqrt(dij);
						double dji = 1.0 + exnji + ssuml;
						double bji = 1.0 / Math.Sqrt(dji);
						double dbdzi = -0.5 * bij / dij;
						double dbdzj = -0.5 * bji / dji;
						double vatt = EXX1[j];
						double dradi = 0;
						double dradj = 0;
						double drdc = 0;
						double conjug = 1.0 + (conk * conk) + (conl * conl);
						double xnt1 = xni[0] + xni[1] - 1.0;
						double xnt2 = xnj[0] + xnj[1] - 1.0;

						double rad = 0;
						radic(ki, kj, xnt1, xnt2, conjug, rad, dradi, dradj, drdc);

						double btot = bji + bij + rad;

						//dihedral terms
						if (kikj + 2 == NDIHED)
						{

							double btor = 0;

							double ator = 0;
							double datori = 0;
							double datorj = 0;
							double datorc = 0;
							tor(xnt1, xnt2, conjug, ator, datori, datorj, datorc);

							if (Math.Abs(ator) > Math.Pow(1.0, -8))
							{
								nk = 0;
								for (int k = JBEGIN; k < JEND; k++)
								{
									if (k != j)
									{
										nk = nk + 1;
										if (Math.Abs(sink[nk - 1]) >= 0.1)
										{
											double sink2 = sink[nk - 1] * sink[nk - 1];
											int kn = JVCT2B[k];

											double[] ck = new double[3] { COR[k, 0], COR[k, 1], COR[k, 2] };
											double rck = RCOR[k];

											double fck = 0;
											double dfck = 0;
											if (Atom.AllAtoms[kn].atomicNumber == 1)
											{
												fck = 1.0;
												dfck = 0;
												if (rck < 1.6)
												{
													if (rck >= 1.3)
													{
														double dtemp = PIDT * (rck - 1.3);
														fck = (1.0 + Math.Cos(dtemp)) / 2.0;
														dfck = -PIDT / 2.0 * Math.Sin(dtemp);
													}
												}
											}
											else
											{
												fck = WW[k];
												dfck = DWW[k];
											}

											nl = 0;
											for (int l = LBEGIN; l < LEND; l++)
											{
												int ln = JVCT2B[l];
												if (ln != i)
												{
													nl = nl + 1;
													if (Math.Abs(sinl[nl - 1]) >= 0.1)
													{
														double sinl2 = sinl[nl - 1] * sinl[nl - 1];
														double[] cl = new double[3] { COR[l, 0], COR[l, 1], COR[l, 2] };
														double rcl = RCOR[l];
														double fcl = 0;
														double dfcl = 0;
														if (Atom.AllAtoms[ln].atomicNumber == 1)
														{
															fcl = 1.0;
															dfcl = 0;
															if (rcl < 1.6)
															{
																if (rcl >= 1.3)
																{
																	double dtemp = PIDT * (rcl - 1.3);
																	fcl = (1.0 + Math.Cos(dtemp)) / 2.0;
																	dfcl = -PIDT / 2.0 * Math.Sin(dtemp);
																}
															}
														}
														else
														{
															fcl = WW[l];
															dfcl = DWW[l];
														}
														double t1 = rck * rcl * SIJ * SIJ * sink[nk - 1] * sinl[nl - 1];
														double dt1dik = 1.0 / rck / rck - dctik[nk - 1] / sink2 * cosk[nk - 1];
														double dt1djk = -dctjk[nk - 1] / sink2 * cosk[nk - 1];
														double dt1djl = 1.0 / rcl / rcl - dctjl[nl - 1] / sinl2 * cosl[nl - 1];
														double dt1dil = -dctil[nl - 1] / sinl2 * cosl[nl];
														double dt1dij = 2.0 / SIJ / SIJ - dctij[nk - 1] / sink2 * cosk[nk - 1] - dctji[nl - 1] / sinl2 * cosl[nl - 1];

														double crkx = ck[1] * CJ[2] - CJ[1] * ck[2];
														double crlx = CJ[1] * cl[2] - cl[1] * CJ[2];
														double crky = ck[2] * CJ[0] - CJ[2] * ck[0];
														double crly = CJ[2] * cl[0] - cl[2] * CJ[0];
														double crkz = ck[0] * CJ[1] - CJ[0] * ck[1];
														double crlz = CJ[0] * cl[1] - cl[0] * CJ[1];

														double t2 = crkx * crlx + crky * crly + crkz * crlz;

														double cw = t2 / t1;
														double bt = 1.0 - cw * cw;
														btor = btor + bt * fck * fcl;

														double[] dt2dik = new double[3] { -CJ[2] * crly + CJ[1] * crlz, -CJ[0] * crlz + CJ[2] * crlx, -CJ[1] * crlx + CJ[0] * crly };
														double[] dt2djl = new double[3] { -CJ[1] * crkz + CJ[2] * crky, -CJ[2] * crkx + CJ[0] * crkz, -CJ[0] * crky + CJ[1] * crkx };
														double[] dt2dij = new double[3] { ck[2] * crly - cl[2] * crky - ck[1] * crlz + cl[2] * crkz, ck[0] * crlz - cl[0] * crkz - ck[2] * crlx + cl[2] * crkx, ck[1] * crlx - cl[1] * crkx - ck[0] * crly + cl[0] * crky };

														double aa = -vatt * 2.0 * cw / t1 * ator * fcl * fck;
														double aaa1 = vatt * bt * ator;
														double at2 = aa * t2;

														double rp1 = -dt1dij * at2;
														double rp2 = -dt1dik * at2 + aaa1 * fcl * dfck / rck;
														double rp3 = -dt1djl * at2 + aaa1 * fck * dfcl / rcl;
														double rp4 = -dt1djk * at2;
														double rp5 = -dt1dil * at2;

														for (int mm = 0; mm < 3; mm++)
														{
															double rep = rp1 * CJ[mm] + aa * dt2dij[mm];
															RNP[i, mm] = RNP[i, mm] + rep;
															RNP[jn, mm] = RNP[jn, mm] - rep;

															rep = rp2 * ck[mm] + aa * dt2dik[mm];
															RNP[i, mm] = RNP[i, mm] + rep;
															RNP[kn, mm] = RNP[kn, mm] - rep;

															rep = rp3 * cl[mm] + aa * dt2djl[mm];
															RNP[jn, mm] = RNP[jn, mm] + rep;
															RNP[ln, mm] = RNP[ln, mm] - rep;

															rep = rp4 * xk[nk - 1][mm];
															RNP[jn, mm] = RNP[jn, mm] + rep;
															RNP[kn, mm] = RNP[kn, mm] - rep;

															rep = rp5 * xl[nl - 1][mm];
															RNP[i, mm] = RNP[i, mm] + rep;
															RNP[ln, mm] = RNP[ln, mm] - rep;
														}

														{
															bool flg = checkNaN(RNP);
															if (flg)
															{
																Debugger.Break();
															}
														}
													}
												}
											}
										}
									}
								}
							}

							btot = btot + btor * ator;
							dradi = dradi + datori * btor;
							dradj = dradj + datorj * btor;
							drdc = drdc + datorc * btor;
						} // end dihedral forces

						tote = tote - btot * vatt;

						double vdbdi = vatt * dbdzi;
						double vdbdj = vatt * dbdzj;
						double vdrdc = vatt * drdc;
						double vdrdi = vatt * dradi;
						double vdrdj = vatt * dradj;

						double rp = vdbdi * xsij + vdbdj * xsji + btot * DEXX1[j];
						for (int mm = 0; mm < 3; mm++)
						{
							double rep = rp * CJ[mm];
							RNP[i, mm] = RNP[i, mm] + rep;
							RNP[jn, mm] = RNP[jn, mm] - rep;
						}

						{
							bool flg = checkNaN(RNP);
							if (flg)
							{
								Debugger.Break();
							}
						}

						// add many-body forces /////////////////////////////////////////////////////////////////////////////
						// I side of bond
						nk = 0;
						for (int k = JBEGIN; k < JEND; k++)
						{
							if (k != j)
							{
								int kn = JVCT2B[k];
								int kk = 5;
								if (Atom.AllAtoms[kn].atomicNumber == 6)
								{
									kk = 0;
								}
								else if (Atom.AllAtoms[kn].atomicNumber == 1)
								{
									kk = 1;
								}
								else
								{
									Console.WriteLine("Atom types should be either carbon or hydrogen to use REBO potential!");
									Console.ReadLine();
									return;
								}

								double dwr = DWW[k] / RCOR[k];
								nk = nk + 1;

								// first neighbors
								double rp1 = vdbdi * (xsik[nk - 1] + dwr * dexni[kk]) + dwr * (vdrdi + vdrdc * cfuni[nk - 1]) + vdbdi * dwr * sdalik;
								double rp2 = vdbdi * xsjk[nk - 1];
								for (int mm = 0; mm < 3; mm++)
								{
									double rep = rp1 * COR[k, mm];
									RNP[i, mm] = RNP[i, mm] + rep;
									RNP[kn, mm] = RNP[kn, mm] - rep;

									// angular forces
									rep = rp2 * xk[nk - 1][mm];
									RNP[jn, mm] = RNP[jn, mm] + rep;
									RNP[kn, mm] = RNP[kn, mm] - rep;
								}

								{
									bool flg = checkNaN(RNP);
									if (flg)
									{
										Debugger.Break();
									}
								}

								// second neighbors via RADIC
								double ddr = vdrdc * dcfuni[nk - 1] * 2.0 * conk;
								if (ddr != 0)
								{
									int MBEGIN = NABORS[kn];
									int MEND = NABORS[kn + 1];
									for (int m = MBEGIN; m < MEND; m++)
									{
										int mn = JVCT2B[m];
										if (mn != kn)
										{
											rp = ddr * DWW[m] / RCOR[m];
											for (int mm = 0; mm < 3; mm++)
											{
												double rep = rp * COR[m, mm];
												RNP[kn, mm] = RNP[kn, mm] + rep;
												RNP[mn, mm] = RNP[mn, mm] - rep;
											}

											{
												bool flg = checkNaN(RNP);
												if (flg)
												{
													Debugger.Break();
												}
											}
										}
									}

								}
							}
						}

						// J side of bond
						nl = 0;
						for (int l = LBEGIN; l < LEND; l++)
						{
							int ln = JVCT2B[l];
							if (ln != i)
							{
								int kl = 5;
								if (Atom.AllAtoms[ln].atomicNumber == 6)
								{
									kl = 0;
								}
								else if (Atom.AllAtoms[ln].atomicNumber == 1)
								{
									kl = 1;
								}
								else
								{
									Console.WriteLine("Atom types should be either carbon or hydrogen to use REBO potential!");
									Console.ReadLine();
									return;
								}

								double dwr = DWW[l] / RCOR[l];
								nl = nl + 1;

								// first neighbors
								double rp1 = vdbdj * (xsjl[nl - 1] + dwr * dexnj[kl]) + dwr * (vdrdj + vdrdc * cfunj[nl - 1]) + vdbdj * dwr * sdaljl;
								double rp2 = vdbdj * xsil[nl - 1];
								for (int mm = 0; mm < 3; mm++)
								{
									double rep = rp1 * COR[l, mm];
									RNP[jn, mm] = RNP[jn, mm] + rep;
									RNP[ln, mm] = RNP[ln, mm] - rep;

									// angular forces
									rep = rp2 * xl[nl - 1][mm];
									RNP[i, mm] = RNP[i, mm] + rep;
									RNP[ln, mm] = RNP[ln, mm] - rep;
								}

								// second neighbors via RADIC
								double ddr = vdrdc * dcfunj[nl - 1] * 2.0 * conl;
								if (ddr != 0)
								{
									int NBEGIN = NABORS[ln];
									int NEND = NABORS[ln + 1];
									for (int n = NBEGIN; n < NBEGIN; n++)
									{
										int nn = JVCT2B[n];
										if (nn != ln)
										{
											rp = ddr * DWW[n] / RCOR[n];
											for (int mm = 0; mm < 3; mm++)
											{
												double rep = rp * COR[n, mm];
												RNP[ln, mm] = RNP[ln, mm] + rep;
												RNP[nn, mm] = RNP[nn, mm] - rep;
											}
										}
									}
								}
							}
						}
					}
				}
			}



			return;
		}

		private void bcuint(int kl, int ki, double xx1, double xx2, int nh, int nc, double ansy, double ansy1, double ansy2) //Bicubic interpolation of CLM; also returns derivatives. 
		{
			ansy = 0;
			ansy1 = 0;
			ansy2 = 0;

			for (int j = 0; j < 16; j++ )
			{
				double x = CLM[ki, nh, nc, j] * (Math.Pow(xx1, IN2[j, 0])) * (Math.Pow(xx2, IN2[j, 1]));
				ansy = ansy + x;
				ansy1 = ansy1 + x * IN2[j, 0] / xx1;
				ansy2 = ansy2 + x * IN2[j, 1] / xx2;
			}

			return;
		}

		private void radic(int ki, int kj, double xnt1, double xnt2, double conjug, double rad, double drdl, double drdm, double drdn) //Tricubic interpolation of CLMN 
		{

			int L = (int)xnt1;
			int M = (int)xnt2;
			int N = (int)conjug;
			rad = 0;
			drdl = 0;
			drdm = 0;
			drdn = 0;
			int kikj = ki + kj;

			if (L > 4)
			{
				L = 4;
				xnt1 = 4.0;
			}

			if (M > 4)
			{
				M = 4;
				xnt2 = 4.0;
			}

			if (N > 9)
			{
				N = 9;
				conjug = 9.0;
			}

			// debug purposes
			// WARNING: this section is not part of the original code and has been added to enable the code to run without IndexOutOfRange error. Serious debuging is needed
			if (L < 1)
			{
				L = 1;
			}

			if (M < 1)
			{
				M = 1;
			}

			if (N < 1)
			{
				N = 1;
			}
			// end of debug code

			for (int j = 0; j < 64; j++ )
			{
				double x = CLMN[kikj, L - 1, M - 1, N - 1, j] * (Math.Pow(xnt1, IN3[j, 0])) * (Math.Pow(xnt2, IN3[j, 1])) * (Math.Pow(conjug, IN3[j, 2]));
				rad = rad + x;
				drdl = drdl + x * IN3[j, 0] / xnt1;
				drdm = drdm + x * IN3[j, 1] / xnt2;
				drdn = drdn + x * IN3[j, 2] / conjug;
			}

			return;
		}

		private void tor(double xnt1, double xnt2, double conjug, double ator, double drdl, double drdm, double drdn) //Torsional interaction TLMN tricubic interpolation
		{

			ator = 0;
			drdl = 0;
			drdm = 0;
			drdn = 0;

			if ((xnt1<=4) && (xnt2<=4))
			{
				int L = (int)xnt1;
				int M = (int)xnt2;
				int N = (int)conjug;

				// debug purposes
				// WARNING: this section is not part of the original code and has been added to enable the code to run without IndexOutOfRange error. Serious debuging is needed
				if (L < 1)
				{
					L = 1;
				}

				if (M < 1)
				{
					M = 1;
				}

				if (N < 1)
				{
					N = 1;
				}

				if (N > 10)
				{
					N = 10;
				}
				// end of debug code

				for (int j = 0; j < 64; j++)
				{
					double x = TLMN[L-1, M-1, N-1, j] * (Math.Pow(xnt1, IN3[j, 0])) * (Math.Pow(xnt2, IN3[j, 1])) * (Math.Pow(conjug, IN3[j, 2]));
					ator = ator + x;
					drdl = drdl + x * IN3[j, 0] / xnt1;
					drdm = drdm + x * IN3[j, 1] / xnt2;
					drdn = drdn + x * IN3[j, 2] / conjug;
				}
			}
			return;
		}
    
		private bool checkNaN (double [,] myArray)
		{
			bool flg = false;
			for (int i = 0; i < myArray.GetLength(0); i++)
			{
				for (int j = 0; j < myArray.GetLength(1); j++)
				{
					if (double.IsNaN(myArray[i, j]))
					{
						Console.WriteLine("Calculated force is NaN!!!");
						flg = true;
					}
				}
			}
			return flg;
		}
	}
}
