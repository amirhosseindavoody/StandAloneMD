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

        //The mesh size for pre-calculating Lennard Jones force.
        private static float dR = 0.0001f;

        //pre-calculated coefficients and forces for Buckingham potential
        private static float[, ,] preBuckinghamAcceleration;
        private static float[, ,] PreBuckinghamPotential;

        private static float[,] coeff_A = new float[3, 3];
        private static float[,] coeff_B = new float[3, 3];
        private static float[,] coeff_C = new float[3, 3];
        private static float[,] coeff_D = new float[3, 3];

        public static void preBuckingham()
        {
            // precalculate the LennardJones potential and store it in preLennarJones array.
            int nR = (int)(cutoff / dR) + 1;
            preBuckinghamAcceleration = new float[3,3,nR];
            PreBuckinghamPotential = new float[3,3,nR];

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

                    for (int iR = 0; i < nR; i++)
                    {
                        float distance = (float)i * dR;
                        preBuckinghamAcceleration[firstAtom.atomID,secondAtom.atomID,iR] = calcForce(distance,firstAtom,secondAtom);
                        PreBuckinghamPotential[firstAtom.atomID, secondAtom.atomID, iR] = calcPotential(distance, firstAtom, secondAtom);
                    }
                }
            }
        }

        //the function returns the LennarJones force on the atom given the list of the atoms that are within range of it
        private static float calcForce(float distance,Atom firsAtom, Atom secondAtom)
        {
            float invDistance2 = 1.0f / distance / distance;
            float invDistance6 = invDistance2 * invDistance2 * invDistance2;
            float invCutoff2 = 1.0f / cutoff / cutoff;
            float invCutoff6 = invCutoff2 * invCutoff2 * invCutoff2;
            float r_min = rMinMultiplier;

            float forceMagnitude = 0.0f;

            if (distance > r_min)
            {
                forceMagnitude = invDistance2 * ((2.0f * invDistance6 * invDistance6 - invDistance6) - (invCutoff2 / invDistance2) * (2.0f * invCutoff6 * invCutoff6 - invCutoff6));
            }
            // Smooth the potential to go to a constant not infinity at r=0
            else
            {
                float invr_min = 1 / r_min;
                float invr_min2 = invr_min * invr_min;
                float invr_min6 = invr_min2 * invr_min2 * invr_min2;
                float magnitude_Vmin = invr_min2 * ((2.0f * invr_min6 * invr_min6 - invr_min6) - (invCutoff2 / invr_min2) * (2.0f * invCutoff6 * invCutoff6 - invCutoff6));

                float r_Vmax = r_min / 1.5f;
                float invr_Vmax2 = 1 / r_Vmax / r_Vmax;
                float invr_Vmax6 = invr_Vmax2 * invr_Vmax2 * invr_Vmax2;
                float magnitude_Vmax = invr_Vmax2 * ((2.0f * invr_Vmax6 * invr_Vmax6 - invr_Vmax6) - (invCutoff2 / invr_Vmax2) * (2.0f * invCutoff6 * invCutoff6 - invCutoff6));

                float part1 = (distance / r_min) * ((float)Math.Exp(distance - r_min));
                float part2 = magnitude_Vmax - magnitude_Vmin;
                forceMagnitude = magnitude_Vmax - (part1 * part2);
            }

            return forceMagnitude;
        }

        //the function returns the LennarJones force on the atom given the list of the atoms that are within range of it
        private static float calcPotential(float distance, Atom firsAtom, Atom secondAtom)
        {
            float invDistance2 = 1.0f / distance / distance;
            float invDistance6 = invDistance2 * invDistance2 * invDistance2;
            float invCutoff2 = 1.0f / cutoff / cutoff;
            float invCutoff6 = invCutoff2 * invCutoff2 * invCutoff2;
            float r_min = rMinMultiplier;

            float potential = 0.0f;

            if (distance > 0.0f)
            {
                potential = 4.0f * ((invDistance6 * invDistance6 - invDistance6) + (6.0f * invCutoff6 * invCutoff6 - 3.0f * invCutoff6) * (invCutoff2 / invDistance2) - 7.0f * invCutoff6 * invCutoff6 + 4.0f * invCutoff6);
            }

            return potential;
        }
    }
}
