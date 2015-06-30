using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace StandAloneMD
{
    class InputOutput
    {
        private static int writeFlag;
        private static StreamWriter positionFile;
        private static StreamWriter energyFile;
        private static StreamWriter temperatureFile;
        private static StreamWriter pairDistributionFile;

        public static void OpenFiles()
        {
            positionFile = new StreamWriter("position.txt");
            energyFile = new StreamWriter("energy.txt");
            temperatureFile = new StreamWriter("temperature.txt");
            pairDistributionFile = new StreamWriter("pairDistribution.txt");

            writeFlag = 01;
        }

        public static void WritePosition()
        {
            writeFlag--;
            if (writeFlag == 0)
            {
                Console.WriteLine("iTime = " + StaticVariables.iTime);
                energyFile.WriteLine(StaticVariables.potentialEnergy + "    " + StaticVariables.kineticEnergy);
				temperatureFile.WriteLine(StaticVariables.currentTemperature);
				//for (int i = 0; i < Atom.AllAtoms.Count; i++)
				//{
				//	positionFile.WriteLine(Atom.AllAtoms[i].position[0] + "    " + Atom.AllAtoms[i].position[1] + "    " + Atom.AllAtoms[i].position[2]);
				//}
                writeFlag = 100;
            }
        }

        public static void WritePairDistribution()
        {   
            for (int i = 0; i < PairDistributionFunction.PairDistributionAverage.Length; i++)
            {
                pairDistributionFile.WriteLine(PairDistributionFunction.PairDistributionAverage[i]);
            }
        }

        public static void CloseFiles()
        {
            positionFile.Close();
            energyFile.Close();
            temperatureFile.Close();
            pairDistributionFile.Close();
        }

        public static void WritePotential(float[, ,] myPotential)
        {
            StreamWriter potentialFile;
            potentialFile = new StreamWriter("potential.txt");
            int numAtomTypes = myPotential.GetLength(0);
            int nR = myPotential.GetLength(2);

            for (int iR = 0; iR < nR; iR++)
            {
                for (int iAtom1 = 0; iAtom1 < numAtomTypes; iAtom1++)
                {
                    for (int iAtom2 = 0; iAtom2 < numAtomTypes; iAtom2++)
                    {
                        potentialFile.WriteLine(myPotential[iAtom1, iAtom2, iR].ToString("E6"));
                    }
                }
            }
            potentialFile.Close();
        }

        public static void WriteForce(float[, ,] myForce)
        {
            StreamWriter forceFile;
            forceFile = new StreamWriter("force.txt");
            int numAtomTypes = myForce.GetLength(0);
            int nR = myForce.GetLength(2);

            for (int iR = 0; iR < nR; iR++)
            {
                for (int iAtom1 = 0; iAtom1 < numAtomTypes; iAtom1++)
                {
                    for (int iAtom2 = 0; iAtom2 < numAtomTypes; iAtom2++)
                    {
                        forceFile.WriteLine(myForce[iAtom1, iAtom2, iR].ToString("E6"));
                    }
                }
            }
            forceFile.Close();
        }

        public static void ReadPotential(float[, ,] myPotential)
        {
            string line;
            int numAtomTypes = myPotential.GetLength(0);
            int nR = myPotential.GetLength(2);
            StreamReader potentialFile;
            potentialFile = new StreamReader("potential.txt");

            
            for (int iR = 0; iR < nR; iR++)
            {
                for (int iAtom1 = 0; iAtom1 < numAtomTypes; iAtom1++)
                {
                    for (int iAtom2 = 0; iAtom2 < numAtomTypes; iAtom2++)
                    {
                        line = potentialFile.ReadLine();
                        if (line == null)
                        {
                            Console.WriteLine("Input file does not match!");
                            Console.ReadLine();
                        }
                        else
                        {
                            myPotential[iAtom1, iAtom2, iR] = float.Parse(line);
                        }
                    }
                }
            }
            line = potentialFile.ReadLine();
            if (line != null)
            {
                Console.WriteLine("Input file does not match!");
                Console.ReadLine();
            }
            else
            {
                Console.WriteLine("Potential read successfully!!!");
            }
        }

        public static void ReadForce(float[, ,] myForce)
        {
            string line;
            int numAtomTypes = myForce.GetLength(0);
            int nR = myForce.GetLength(2);
            StreamReader potentialFile;
            potentialFile = new StreamReader("force.txt");


            for (int iR = 0; iR < nR; iR++)
            {
                for (int iAtom1 = 0; iAtom1 < numAtomTypes; iAtom1++)
                {
                    for (int iAtom2 = 0; iAtom2 < numAtomTypes; iAtom2++)
                    {
                        line = potentialFile.ReadLine();
                        if (line == null)
                        {
                            Console.WriteLine("Input file does not match!");
                            Console.ReadLine();
                        }
                        else
                        {
                            myForce[iAtom1, iAtom2, iR] = float.Parse(line);
                        }
                    }
                }
            }
            line = potentialFile.ReadLine();
            if (line != null)
            {
                Console.WriteLine("Input file does not match!");
                Console.ReadLine();
            }
            else
            {
                Console.WriteLine("Potential read successfully!!!");
            }
        }

    }
}
