using System;
using System.IO;
using System.Diagnostics;

namespace StandAloneMD
{
	class MainClass
	{
		public static void Main (string[] args)
		{
            WriteData myData = new WriteData();

            CreateEnvironment.myEnvironment = new CreateEnvironment();
            CreateEnvironment.myEnvironment.PreCompute();
            //CreateEnvironment.myEnvironment.InitAtoms();
            CreateEnvironment.myEnvironment.InitAtomsDebug();
            Potential.myPotential.calculateVerletRadius();

            Console.WriteLine("Number of atoms = " + CreateEnvironment.myEnvironment.numAtoms);

            float totalTime = 40000.0f * StaticVariables.MDTimestep;
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();
            while (StaticVariables.currentTime < totalTime)
            {            
                PhysicsEngine.VelocityVerlet();
                PhysicsEngine.ReflectFromWalls();
                PhysicsEngine.CalculateEnergy();
                //PairDistributionFunction.calculatePairDistribution();

                StaticVariables.currentTime = StaticVariables.currentTime + StaticVariables.MDTimestep;
                StaticVariables.iTime++;
                myData.WritePosition();
            }
            
            stopwatch.Stop();
            //myData.WritePairDistribution();
            Console.WriteLine("iTime = " + StaticVariables.iTime + "            Current Time = " + StaticVariables.currentTime);
            Console.WriteLine("Time elapsed: {0}", stopwatch.Elapsed);
            Console.ReadLine();
		}
	}
}
