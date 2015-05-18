using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace StandAloneMD
{
    class REBO : Potential
    {
		private int NDIHED;
		private int NTAB = 10000; // this is the array size of potential lookup table.
		private int KEND; // this is the total number of neighbor pairs that are calculated and stored in LIST, IVCT2B, and JVCT2B
        private int[] IGC, IGH;
		private int[] NABORS;
		private List<int> LIST = new List<int>();
		private List<int> IVCT2B = new List<int>();
		private List<int> JVCT2B = new List<int>();
		private int[] LCHECK;
        private int[,] IN2, IN3;
		private double RLL = 0.5; // this is most probably the shell thickness used for neighbor list calculation
        private double SIGMA, EPSI, XQ, ATT, XQM, PQ, PIDT, XXDB;
		private double tote;
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
			tote = 0.0;
			return 0.0f;
			
        }

        public override void calculateVerletRadius()
        {
            
        }

        //This function creates a list of all neighbor list for each atom
        public override void calculateNeighborList()
        {
			neighborListFlag = true;
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

		private void mtable()
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

		public void caguts()
		{
			// calculate two-body forces and neighbor list for hydrocarbons
			double[] RR = new double[3];
			double[] RI = new double[3];
			int[] NABORS = new int[Atom.AllAtoms.Count + 1];

			RNP = new double[Atom.AllAtoms.Count, 3];
			Array.Clear(RNP, 0, Atom.AllAtoms.Count * 3);



			eatom = new double[Atom.AllAtoms.Count];
			Array.Clear(eatom, 0, eatom.Length);

			if (neighborListFlag)
			{
				// set up neighbor list
				int K = 0;
				LIST.Clear();
				IVCT2B.Clear();
				JVCT2B.Clear();
				for (int I = 0; I<Atom.AllAtoms.Count; I++)
				{
					NABORS[I] = K + 1;
					RI[0] = Atom.AllAtoms[I].position[0];
					RI[1] = Atom.AllAtoms[I].position[1];
					RI[2] = Atom.AllAtoms[I].position[2];
					int ki = 5;
					if (Atom.AllAtoms[I].atomicNumber == 6)
					{
						ki = 0;
					}
					else if(Atom.AllAtoms[I].atomicNumber == 1)
					{
						ki = 1;
					}
					else
					{
						Console.WriteLine("Atom types should be either carbon or hydrogen to use REBO potential!");
						Console.ReadLine();
						return;
					}

					for (int J = 0; J<Atom.AllAtoms.Count; J++)
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
							for (int L = 0; L < 3; L++ )
							{
								RR[L] = RI[L] - Atom.AllAtoms[J].position[L];
								rsq = rsq + RR[L] * RR[L];
							}
							
							if (rsq <= RLIS)
							{
								K++;
								LIST.Add(J);
								IVCT2B.Add(I);
								JVCT2B.Add(J);
							}


						}
					}
				}
				NABORS[Atom.AllAtoms.Count] = K+1;
				KEND = K;
			}

			//debug
			if (KEND != IVCT2B.Count)
			{
				Console.WriteLine("There is a problem with index number of IVCT2B versus KEND! Please investigate!");
				Console.ReadLine();
			}
			//debug


			LCHECK = new int[KEND];
			Array.Clear(LCHECK,0, KEND);
			COR = new double [KEND, 3];
			Array.Clear(COR,0, 3*KEND);
			RCOR = new double [KEND];
			Array.Clear(RCOR,0, KEND);
			WW = new double [KEND];
			Array.Clear(WW,0, KEND);
			DWW = new double [KEND];
			Array.Clear(DWW,0, KEND);
			EXX1 = new double [KEND];
			Array.Clear(EXX1,0, KEND);
			DEXX1 = new double [KEND];
			Array.Clear(DEXX1,0, KEND);
			double[,] RPP = new double[KEND,3];
			Array.Clear(RPP,0, 3*KEND);

			for (int K=0; K<KEND; K++)
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
				
				LCHECK[K]=0;
				double rsq = 0.0;
				for (int L=0; L<3; L++)
				{
					RR[L] = Atom.AllAtoms[I].position[L] - Atom.AllAtoms[J].position[L];
					rsq = rsq + RR[L] * RR[L];
					COR[K,L] = RR[L];
				}
				
				if (rsq <= RMAX[KI,KJ])
				{

					if ((KJ <= 1) && (KI <= 1)) LCHECK[K] = 1;
					if ((KJ >= 2) && (KI >= 2)) LCHECK[K] = 2;

					double rc = Math.Sqrt(rsq);
					double rt = rc / DDTAB[KI, KJ];
					int it = Math.Min((int)rt, NTAB - 2);

					RCOR[K]=rc;
					WW[K] = TABFC[KI, KJ, it] + (TABFC[KI,KJ,it+1]-TABFC[KI,KJ,it])*(rt-(double)it+1.0);
					DWW[K] = TABDFC[KI, KJ, it] + (TABDFC[KI, KJ, it + 1] - TABDFC[KI, KJ, it]) * (rt - (double)it + 1.0);
					EXX1[K] = ATABLE[KI, KJ, it] + (ATABLE[KI,KJ,it+1]-ATABLE[KI,KJ,it])*(rt-(double)it+1.0);
					DEXX1[K] = DATABLE[KI, KJ, it] + (DATABLE[KI, KJ, it + 1] - DATABLE[KI, KJ, it]) * (rt - (double)it + 1.0);
					
					if (I < J)
					{
						double vv = RTABLE[KI, KJ, it] + (RTABLE[KI, KJ, it + 1] - RTABLE[KI, KJ, it]) * (rt - (double)it + 1.0); // this is the potential energy between atom pair I and J
						double rp = DRTABLE[KI, KJ, it] + (DRTABLE[KI, KJ, it + 1] - DRTABLE[KI, KJ, it]) * (rt - (double)it + 1.0);
						tote = tote + vv;
						eatom[I] = eatom[I] + vv / 2.0;
						eatom[J] = eatom[J] + vv / 2.0;

						for (int L=0; L<3; L++)
						{
							RPP[K, L] = rp * RR[L];
						}
					}
				}
			}

			for (int K=0;K<KEND; K++)
			{
				if (LCHECK[K] != 0)
				{
					int I = IVCT2B[K];
					int J = JVCT2B[K];
					if (I < J)
					{
						for (int L=0; L<3; L++)
						{
							RNP[I, L] = RNP[I, L] + RPP[K, L];
							RNP[J, L] = RNP[J, L] + RPP[K, L];
						}
					}
				}
			}

			// add a check to see if the atoms are carbohydrates or si-germanium and perform the corresponding routines
			// call pibond for carbohydrates
			// call sili_germ for silicon germanium atoms.
		}
    }
}
