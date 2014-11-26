using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace StandAloneMD
{
    class Buckingham : Potential
    {
        //Cutoff distance for calculating LennarJones force. This quantity is unit less and normalized to sigmaValue for atom pair
        public float cutoff = 10.0f; //[Angstroms]
        public float cutoffSqr;

        //The mesh size for pre-calculating Lennard Jones force.
        private float dR = 0.00001f;

        //pre-calculated coefficients and forces for Buckingham potential
        private float[, ,] preBuckinghamAcceleration;
        private float[, ,] PreBuckinghamPotential;

        private float[,] coeff_A = new float[3, 3];
        private float[,] coeff_B = new float[3, 3];
        private float[,] coeff_C = new float[3, 3];
        private float[,] coeff_D = new float[3, 3];

        public Buckingham()
        {
            cutoffSqr = cutoff * cutoff;
        }

        public override void preCompute()
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

                    for (int iR = 0; iR < nR; iR++)
                    {
                        float distance = (float)iR * dR;
                        if (distance < 0.5f)
                            distance = 0.5f;
                        preBuckinghamAcceleration[firstAtom.atomID,secondAtom.atomID,iR] = calcAcceleration(distance,firstAtom,secondAtom);
                        PreBuckinghamPotential[firstAtom.atomID, secondAtom.atomID, iR] = calcPotential(distance, firstAtom, secondAtom);
                    }
                }
            }
        }

        //the function returns the LennarJones force on the atom given the list of the atoms that are within range of it
        private float calcAcceleration(float distance,Atom firstAtom, Atom secondAtom)
        {
            float invDistance2 = 1.0f / distance / distance;
            float invDistance6 = invDistance2 * invDistance2 * invDistance2;
            float invDistance7 = invDistance6 / distance;
            float invDistance8 = invDistance2 * invDistance2 * invDistance2 * invDistance2;
            float invDistance9 = invDistance8 / distance;
            float invCutoff2 = 1.0f / cutoff / cutoff;
            float invCutoff6 = invCutoff2 * invCutoff2 * invCutoff2;
            float invCutoff7 = invCutoff6 / cutoff;
            float invCutoff8 = invCutoff2 * invCutoff2 * invCutoff2 * invCutoff2;
            float invCutoff9 = invCutoff8 / cutoff;
            
            float A = coeff_A [firstAtom.atomID,secondAtom.atomID];
            float B = coeff_B [firstAtom.atomID,secondAtom.atomID];
            float C = coeff_C [firstAtom.atomID,secondAtom.atomID];
            float D = coeff_D [firstAtom.atomID,secondAtom.atomID];

            float uPrime_r = -A * B * (float)Math.Exp(-B * distance) / StaticVariables.angstromsToMeters + 6.0f * C * invDistance7 / StaticVariables.angstromsToMeters + 8.0f * D * invDistance9 / StaticVariables.angstromsToMeters - firstAtom.Q_eff * secondAtom.Q_eff / (4.0f * StaticVariables.epsilon0 * (float)Math.PI * StaticVariables.angstromsToMeters * StaticVariables.angstromsToMeters) * invDistance2;
            float uPrime_rc = -A * B * (float)Math.Exp(-B * cutoff) / StaticVariables.angstromsToMeters + 6.0f * C * invCutoff7 / StaticVariables.angstromsToMeters + 8.0f * D * invCutoff9 / StaticVariables.angstromsToMeters - firstAtom.Q_eff * secondAtom.Q_eff / (4.0f * StaticVariables.epsilon0 * (float)Math.PI * StaticVariables.angstromsToMeters * StaticVariables.angstromsToMeters) * invCutoff2;

            //float forceMagnitude = -1.0f * uPrime_r / distance + uPrime_rc / cutoff;
            float forceMagnitude = -1.0f * uPrime_r / distance;
            float acceleration = forceMagnitude / (firstAtom.massamu * StaticVariables.amuToKg * StaticVariables.angstromsToMeters); //Units of [1 / second^2] when multiplied by deltaR gets units of [Angstrom / second^2]
            return acceleration;
        }

        //the function returns the LennarJones force on the atom given the list of the atoms that are within range of it
        private float calcPotential(float distance, Atom firstAtom, Atom secondAtom)
        {
            float invDistance2 = 1.0f / distance / distance;
            float invDistance6 = invDistance2 * invDistance2 * invDistance2;
            float invDistance7 = invDistance6 / distance;
            float invDistance8 = invDistance2 * invDistance2 * invDistance2 * invDistance2;
            float invDistance9 = invDistance8 / distance;
            float invCutoff2 = 1.0f / cutoff / cutoff;
            float invCutoff6 = invCutoff2 * invCutoff2 * invCutoff2;
            float invCutoff7 = invCutoff6 / cutoff;
            float invCutoff8 = invCutoff2 * invCutoff2 * invCutoff2 * invCutoff2;
            float invCutoff9 = invCutoff8 / cutoff;

            float A = coeff_A[firstAtom.atomID, secondAtom.atomID];
            float B = coeff_B[firstAtom.atomID, secondAtom.atomID];
            float C = coeff_C[firstAtom.atomID, secondAtom.atomID];
            float D = coeff_D[firstAtom.atomID, secondAtom.atomID];

            float uPrime_r = -A * B * (float)Math.Exp(-B * distance) / StaticVariables.angstromsToMeters + 6.0f * C * invDistance7 / StaticVariables.angstromsToMeters + 8.0f * D * invDistance9 / StaticVariables.angstromsToMeters - firstAtom.Q_eff * secondAtom.Q_eff / (4.0f * StaticVariables.epsilon0 * (float)Math.PI * StaticVariables.angstromsToMeters * StaticVariables.angstromsToMeters) * invDistance2;
            float uPrime_rc = -A * B * (float)Math.Exp(-B * cutoff) / StaticVariables.angstromsToMeters + 6.0f * C * invCutoff7 / StaticVariables.angstromsToMeters + 8.0f * D * invCutoff9 / StaticVariables.angstromsToMeters - firstAtom.Q_eff * secondAtom.Q_eff / (4.0f * StaticVariables.epsilon0 * (float)Math.PI * StaticVariables.angstromsToMeters * StaticVariables.angstromsToMeters) * invCutoff2;

            float u_r = A * (float)Math.Exp(-B * distance) - C * invDistance6 - D * invDistance8 + firstAtom.Q_eff * secondAtom.Q_eff / (4.0f * StaticVariables.epsilon0 * (float)Math.PI * StaticVariables.angstromsToMeters) / distance;
            float u_rc = A * (float)Math.Exp(-B * cutoff) - C * invCutoff6 - D * invCutoff8 + firstAtom.Q_eff * secondAtom.Q_eff / (4.0f * StaticVariables.epsilon0 * (float)Math.PI * StaticVariables.angstromsToMeters) / cutoff;

            //float potential = u_r - (uPrime_rc * cutoff * StaticVariables.angstromsToMeters / 2.0f) * (distance * distance / cutoff / cutoff) - u_rc + (uPrime_rc * cutoff * StaticVariables.angstromsToMeters / 2.0f) ; //Units of Joules
            float potential = u_r; //Units of Joules
            return potential;
        }

        //the function returns the Lennard-Jones force on the atom given the list of all the atoms in the simulation
        public override void getForce(Atom firstAtom, Atom secondAtom)
        {
            float[] firstAtomAcceleration = new float[3];
            float[] secondAtomAcceleration = new float[3];

            float[] deltaR = new float[3];
            for (int idx = 0; idx < 3; idx++)
            {
                deltaR[idx] = firstAtom.position[idx] - secondAtom.position[idx];
            }
            float distanceSqr = deltaR[0] * deltaR[0] + deltaR[1] * deltaR[1] + deltaR[2] * deltaR[2];

            //only get the forces of the atoms that are within the cutoff range
            if (distanceSqr <= cutoffSqr)
            {
                int iR = (int)((float)Math.Sqrt(distanceSqr) / (dR));
                for (int idx = 0; idx < 3; idx++)
                {
                    firstAtom.accelerationNew[idx] = firstAtom.accelerationNew[idx] + preBuckinghamAcceleration[firstAtom.atomID, secondAtom.atomID,iR] * deltaR[idx];
                    secondAtom.accelerationNew[idx] = secondAtom.accelerationNew[idx] - preBuckinghamAcceleration[secondAtom.atomID, firstAtom.atomID,iR] * deltaR[idx];
                }
            }
        }

        //the function returns the Lennard-Jones force on the atom given the list of all the atoms in the simulation
        public override float getPotential(Atom firstAtom, Atom secondAtom)
        {
            float potential = 0.0f;
            float[] deltaR = new float[3];
            for (int idx = 0; idx < 3; idx++)
            {
                deltaR[idx] = firstAtom.position[idx] - secondAtom.position[idx];
            }
            float distanceSqr = deltaR[0] * deltaR[0] + deltaR[1] * deltaR[1] + deltaR[2] * deltaR[2];

            //only get the forces of the atoms that are within the cutoff range
            if (distanceSqr <= cutoffSqr)
            {
                int iR = (int)((float)Math.Sqrt(distanceSqr) / (dR));
                potential = (PreBuckinghamPotential[firstAtom.atomID, secondAtom.atomID,iR] + PreBuckinghamPotential[firstAtom.atomID, secondAtom.atomID,iR]) / 2.0f ;
            }
            return potential;
        }

        public override void calculateVerletRadius()
        {
            for (int i = 0; i < Atom.AllAtoms.Count - 1; i++)
            {
                Atom currAtom = Atom.AllAtoms[i];
                currAtom.verletRadius = cutoff + 1.0f;
            }
        }

        //This function creates a list of all neighbor list for each atom
        public override void calculateNeighborList()
        {
            //clear the old neighborList
            for (int i = 0; i < Atom.AllAtoms.Count - 1; i++)
            {
                Atom currAtom = Atom.AllAtoms[i];
                currAtom.neighborList.Clear();
            }

            //create the new neighborList
            for (int i = 0; i < Atom.AllAtoms.Count - 1; i++)
            {
                Atom firstAtom = Atom.AllAtoms[i];
                for (int j = i + 1; j < Atom.AllAtoms.Count; j++)
                {
                    Atom secondAtom = Atom.AllAtoms[j];
                    float[] deltaR = new float[3];
                    for (int idx = 0; idx < 3; idx++)
                    {
                        deltaR[idx] = firstAtom.position[idx] - secondAtom.position[idx];
                    }
                    float distanceSqr = deltaR[0] * deltaR[0] + deltaR[1] * deltaR[1] + deltaR[2] * deltaR[2];
                    if (distanceSqr < (firstAtom.verletRadius * firstAtom.verletRadius))
                    {
                        firstAtom.neighborList.Add(secondAtom);
                    }
                }
            }
        }
    }
}
