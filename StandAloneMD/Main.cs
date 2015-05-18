using System;
using System.IO;
using System.Diagnostics;
using System.Collections.Generic;

namespace StandAloneMD
{
	class MainClass
	{
		public static void Main (string[] args)
		{
            Console.WriteLine("Testing REBO potential!!!");
            Potential.myPotential = new REBO();

            Console.WriteLine("Reached end of REBO potential!!!");
            Console.WriteLine("********************STOP THE PROGRAM NOW*********************");
            Console.ReadLine();

            Console.ReadLine();
            InputOutput.OpenFiles();

            CreateEnvironment.myEnvironment = new CreateEnvironment();
            CreateEnvironment.myEnvironment.PreCompute();
            CreateEnvironment.myEnvironment.InitAtoms();
            //CreateEnvironment.myEnvironment.InitAtomsDebug();
            Potential.myPotential.calculateVerletRadius();

            Console.WriteLine("Number of atoms = " + Atom.AllAtoms.Count);
            Console.ReadLine();

            float totalTime = 20000.0f * StaticVariables.MDTimestep;
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();
            while (StaticVariables.currentTime < totalTime)
            {            
                PhysicsEngine.VelocityVerlet();
                BoundaryCondition.myBoundary.Apply();
                PhysicsEngine.CalculateEnergy();
                PairDistributionFunction.calculateAveragePairDistribution();

                StaticVariables.currentTime = StaticVariables.currentTime + StaticVariables.MDTimestep;
                StaticVariables.iTime++;
                InputOutput.WritePosition();
            }
            
            stopwatch.Stop();
            InputOutput.WritePairDistribution();
            InputOutput.CloseFiles();
            Console.WriteLine("iTime = " + StaticVariables.iTime + "            Current Time = " + StaticVariables.currentTime);
            Console.WriteLine("Time elapsed: {0}", stopwatch.Elapsed);
            Console.ReadLine();
		}
	}
}
