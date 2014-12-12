using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace StandAloneMD
{
    class WriteData
    {
        int writeFlag;
        StreamWriter positionFile;
        StreamWriter energyFile;
        StreamWriter temperatureFile;

        public WriteData()
        {
            positionFile = new StreamWriter("position.txt");
            energyFile = new StreamWriter("energy.txt");
            temperatureFile = new StreamWriter("temperature.txt");

            writeFlag = 01;
        }

        public void WritePosition()
        {
            writeFlag--;
            if (writeFlag == 0)
            {
                Console.WriteLine("iTime = " + StaticVariables.iTime);
                for (int i = 0; i < Atom.AllAtoms.Count; i++)
                {
                    positionFile.WriteLine(Atom.AllAtoms[i].position[0] + "    " + Atom.AllAtoms[i].position[1] + "    " + Atom.AllAtoms[i].position[2]);
                    energyFile.WriteLine(StaticVariables.potentialEnergy + "    " + StaticVariables.kineticEnergy);
                    temperatureFile.WriteLine(StaticVariables.currentTemperature);
                }
                writeFlag = 20;
            }
        }

        public static void WritePairDistribution()
        {
            StreamWriter pairDistributionFile;
            pairDistributionFile = new StreamWriter("pairDistribution.txt");
            for (int i = 0; i < PairDistributionFunction.PairDistributionAverage.Length; i++)
            {
                pairDistributionFile.WriteLine(PairDistributionFunction.PairDistributionAverage[i]);
            }
        }

        public static void WritePotential(float[,,] myPotential)
        {
            StreamWriter potentialFile;
            potentialFile = new StreamWriter("potential.txt");
            int numAtomTypes = myPotential.GetLength(0);
            int nR = myPotential.GetLength(2);

            for (int iR = 0; iR < nR; iR++)
            {
                potentialFile.WriteLine(myPotential[0,0,iR] + "     " + myPotential[0,1,iR] + "     " + myPotential[1,0,iR] + "     " + myPotential[1,1,iR]);
            }
        }

        public static void WriteForce(float[, ,] myForce)
        {
            StreamWriter forceFile;
            forceFile = new StreamWriter("force.txt");
            int numAtomTypes = myForce.GetLength(0);
            int nR = myForce.GetLength(2);

            for (int iR = 0; iR < nR; iR++)
            {
                forceFile.WriteLine(myForce[0, 0, iR] + "     " + myForce[0, 1, iR] + "     " + myForce[1, 0, iR] + "     " + myForce[1, 1, iR]);
            }
        }
    }
}
