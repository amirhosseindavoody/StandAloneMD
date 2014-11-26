using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace StandAloneMD
{
    class PairDistributionFunction
    {
        //The mesh size for pair distribution function.
        private static float dR = 0.1f; //[Angstrom]
        private static float maxR = 60.0f; //[Angstrom]
        private static float[] pairDistribution = new float[(int)(maxR / dR)];
        private static float[] pairDistributionAverage = new float[(int)(maxR / dR)];
        private static float numberOfCalculations = 0.0f;
        private static float normCoefficient = 0.0f;

        public static float[] PairDistributionAverage
        {
            get
            {
                return pairDistributionAverage;
            }
        }

        public static void calculatePairDistribution()
        {
            if ((StaticVariables.iTime % StaticVariables.nVerlet == 0) && (StaticVariables.iTime > 20000))
            {
                normCoefficient = CreateEnvironment.myEnvironment.volume * CreateEnvironment.myEnvironment.volume / (Atom.AllAtoms.Count * 4.0f * (float)Math.PI * dR * dR * dR);
                updatePairDistribution();
                for (int iR = 0; iR < pairDistribution.Length; iR++)
                {
                    pairDistributionAverage[iR] = (pairDistributionAverage[iR] * numberOfCalculations + pairDistribution[iR] * (normCoefficient / (float)iR / (float)iR));
                    numberOfCalculations++;
                    pairDistributionAverage[iR] = pairDistributionAverage[iR] / numberOfCalculations;
                }
            }
        }

        private static void updatePairDistribution()
        {
            pairDistribution.Initialize();

            //create the new neighborList
            for (int i = 0; i < Atom.AllAtoms.Count-1; i++)
            {
                Atom firstAtom = Atom.AllAtoms[i];
                for (int j = i+1; j < Atom.AllAtoms.Count; j++)
                {

                    Atom secondAtom = Atom.AllAtoms[j];
                    float[] deltaR = new float[3];
                    for (int idx = 0; idx < 3; idx++)
                    {
                        deltaR[idx] = firstAtom.position[idx] - secondAtom.position[idx];
                    }
                    float distance = (float) Math.Sqrt(deltaR[0] * deltaR[0] + deltaR[1] * deltaR[1] + deltaR[2] * deltaR[2]);
                    int iR = (int)(distance / dR);
                    pairDistribution[iR]++;
                }
            }
        }
    }
}
