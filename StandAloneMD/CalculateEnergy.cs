using System;
using System.Collections.Generic;
using System.IO;

namespace StandAloneMD
{
    public static class CalculateEnergy
    {
        public static void CalculateKineticEnergy()
        {
            StaticVariables.KineticEnergy = 0.0f;
            for (int i = 0; i < Atom.AllAtoms.Count; i++)
            {
                Atom currAtom = Atom.AllAtoms[i];
                float velocitySqr = currAtom.velocity[0] * currAtom.velocity[0] + currAtom.velocity[1] * currAtom.velocity[1] + currAtom.velocity[2] * currAtom.velocity[2];
                
                StaticVariables.KineticEnergy += 0.5f * velocitySqr / currAtom.massamu / StaticVariables.amuToKg;
            }
        }

        public static void CalculatePotentialEnergy()
        {
            StaticVariables.PotentialEnergy = 0.0f;

            for (int i = 0; i < Atom.AllAtoms.Count-1; i++)
            {
                Atom firstAtom = Atom.AllAtoms[i];
                for (int j = i+1; j < Atom.AllAtoms.Count; j++)
                {
                    Atom secondAtom = Atom.AllAtoms[j];
                    StaticVariables.PotentialEnergy +=GetLennarJonesPotential(firstAtom, secondAtom);
                }
            }

        }

        private static float GetLennarJonesPotential(Atom firstAtom, Atom secondAtom)
        {
            float LJPotential = 0.0f;
            
            float[] deltaR = new float[3];
            for (int idx = 0; idx < 3; idx++)
            {
                deltaR[idx] = firstAtom.position[idx]-secondAtom.position[idx];
            }

            float distanceSqr = deltaR[0] * deltaR[0] + deltaR[1] * deltaR[1] + deltaR[2] * deltaR[2];
            float normDistanceSqr = distanceSqr / StaticVariables.sigmaValues[firstAtom.atomID, secondAtom.atomID] / StaticVariables.sigmaValues[firstAtom.atomID, secondAtom.atomID];

            if (normDistanceSqr < StaticVariables.cutoffSqr)
            {
                LJPotential = 4.0f * firstAtom.epsilon * (1/(normDistanceSqr*normDistanceSqr*normDistanceSqr*normDistanceSqr*normDistanceSqr*normDistanceSqr) - 1/(normDistanceSqr*normDistanceSqr*normDistanceSqr));
            }

            return LJPotential;
        }
    }
}
