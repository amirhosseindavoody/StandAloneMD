using System;
using System.IO;
using System.Diagnostics;

namespace StandAloneMD
{
	class MainClass
	{
		public static void Main (string[] args)
		{
            StaticVariables.myEnvironment = new CreateEnvironment();
            StaticVariables.myEnvironment.PreCompute ();
            StaticVariables.myEnvironment.InitAtoms ();
            PhysicsEngine.calculateVerletRadius();

            Console.WriteLine("Number of atoms = " + StaticVariables.myEnvironment.numAtoms);
            //Console.ReadLine();

            float totalTime = 10.0f * (float)Math.Pow(10, -12);
            int iTime = 0;
            //int writeFlag = 01;

            //using (StreamWriter file = new StreamWriter("position.txt"))
            //{
                Stopwatch stopwatch = new Stopwatch();
                stopwatch.Start();
                while (StaticVariables.currentTime < totalTime)
                {
                    if (iTime % StaticVariables.nVerlet == 0)
                        PhysicsEngine.calculateNeighborList();                    
                    PhysicsEngine.VelocityVerlet();
                    PhysicsEngine.ReflectFromWalls();
                    PhysicsEngine.CalculateEnergy();

                    StaticVariables.currentTime = StaticVariables.currentTime + StaticVariables.MDTimestep;
                    iTime++;

                    /*
                    writeFlag--;
                    if (writeFlag == 0)
                    {
                        for (int i = 0; i < Atom.AllAtoms.Count; i++)
                        {
                            file.WriteLine(Atom.AllAtoms[i].position[0] + "    " + Atom.AllAtoms[i].position[1] + "    " + Atom.AllAtoms[i].position[2] + "    " + StaticVariables.potentialEnergy + "    " + StaticVariables.kineticEnergy + "    " + StaticVariables.currentTemperature);
                        }
                        writeFlag = 20;
                    }
                    */


                    //Console.WriteLine("iTime = " + iTime + "            Current Time = " + StaticVariables.currentTime);

                }
                stopwatch.Stop();
                Console.WriteLine("Time elapsed: {0}", stopwatch.Elapsed);
                Console.WriteLine("iTime = " + iTime + "            Current Time = " + StaticVariables.currentTime);
                Console.ReadLine();
            //}
		}
	}
}
