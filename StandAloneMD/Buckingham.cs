using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace StandAloneMD
{
    class Buckingham
    {
        //Cutoff distance for calculating LennarJones force. This quantity is unit less and normalized to sigmaValue for atom pair
        public static float cutoff = 4.0f; //[Angstroms]
        public static float cutoffSqr = cutoff * cutoff;

        //pre-calculated coefficients and forces for Buckingham potential
        public static float[, ,] preBuckinghamAcceleration;
        public static float[, ,] PreBuckinghamPotential;

        public static float[,] coeff_A = new float[3, 3];
        public static float[,] coeff_B = new float[3, 3];
        public static float[,] coeff_C = new float[3, 3];
        public static float[,] coeff_D = new float[3, 3];

        public static void preBuckingham()
        {
            //precompute sigma and acceleration coefficient for the Buckingham potential
            for (int i = 0; i < Atom.templateAtoms.Count; i++)
            {
                Atom firstAtom = Atom.templateAtoms[i];
                for (int j = 0; j < Atom.templateAtoms.Count; j++)
                {
                    Atom secondAtom = Atom.templateAtoms[j];

                    float currentA = (float)Math.Sqrt(firstAtom.buck_A * secondAtom.buck_A);
                    coeff_A[firstAtom.atomID, secondAtom.atomID] = currentA;

                    float currentB = (float)Math.Sqrt(firstAtom.buck_B * secondAtom.buck_B);
                    coeff_B[firstAtom.atomID, secondAtom.atomID] = currentB;

                    float currentC = (float)Math.Sqrt(firstAtom.buck_C * secondAtom.buck_C);
                    coeff_C[firstAtom.atomID, secondAtom.atomID] = currentC;

                    float currentD = (float)Math.Sqrt(firstAtom.buck_D * secondAtom.buck_D);
                    coeff_D[firstAtom.atomID, secondAtom.atomID] = currentD;
                }
            }
        }
    }
}
