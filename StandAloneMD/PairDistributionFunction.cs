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

        public static void calculateAveragePairDistribution()
        {
            normCoefficient = CreateEnvironment.myEnvironment.volume / ((float)Atom.AllAtoms.Count * (float)Atom.AllAtoms.Count * 4.0f * (float)Math.PI * dR * dR * dR);
            for (int iR = 0; iR < pairDistribution.Length; iR++)
            {
                pairDistributionAverage[iR] = (pairDistributionAverage[iR] * numberOfCalculations + pairDistribution[iR] * normCoefficient / (float)iR / (float)iR);
                numberOfCalculations++;
                pairDistributionAverage[iR] = pairDistributionAverage[iR] / numberOfCalculations;
            }
        }

        public static void updatePairDistribution(float distance)
        {
            pairDistribution.Initialize();
            int iR = (int)Math.Floor(distance / dR);
            pairDistribution[iR] += 2.0f;
        }
    }
}
