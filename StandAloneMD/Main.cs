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

            StaticVariables.myEnvironment = new CreateEnvironment();
            StaticVariables.myEnvironment.PreCompute ();
            StaticVariables.myEnvironment.InitAtoms ();

            if (StaticVariables.currentPotential == StaticVariables.Potential.LennardJones)
                LennardJones.calculateVerletRadius();
            if (StaticVariables.currentPotential == StaticVariables.Potential.Buckingham)
                Buckingham.calculateVerletRadius();

            Console.WriteLine("Number of atoms = " + StaticVariables.myEnvironment.numAtoms);

            float totalTime = 40000.0f * StaticVariables.MDTimestep;
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();
            while (StaticVariables.currentTime < totalTime)
            {            
                PhysicsEngine.VelocityVerlet();
                PhysicsEngine.ReflectFromWalls();
                PhysicsEngine.CalculateEnergy();
                PairDistributionFunction.calculatePairDistribution();

                StaticVariables.currentTime = StaticVariables.currentTime + StaticVariables.MDTimestep;
                StaticVariables.iTime++;
                myData.WritePosition();
                //Console.WriteLine("iTime = " + StaticVariables.iTime);
            }
            
            stopwatch.Stop();
            myData.WritePairDistribution();
            Console.WriteLine("iTime = " + StaticVariables.iTime + "            Current Time = " + StaticVariables.currentTime);
            Console.WriteLine("Time elapsed: {0}", stopwatch.Elapsed);
            Console.ReadLine();
		}
	}
}
