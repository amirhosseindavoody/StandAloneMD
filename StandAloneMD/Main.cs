using System;
using System.IO;

namespace StandAloneMD
{
	class MainClass
	{
		public static void Main (string[] args)
		{
            CreateEnvironment environment = new CreateEnvironment();
            Console.WriteLine("Number of atoms = " + environment.numAtoms);

            environment.PreCompute ();
            environment.InitAtoms ();

            float totalTime = 100.0f * (float) Math.Pow (10,-12);
            int iTime = 0;
            int writeFlag = 01;

            using (StreamWriter file = new StreamWriter("force.txt"))
            {
                for (int i = 0; i < StaticVariables.preLennardJones.Length; i++)
                {
                    file.WriteLine(StaticVariables.preLennardJones[i]);
                }
            }

            using (StreamWriter file = new StreamWriter("position.txt"))
            {
                while (StaticVariables.currentTime < totalTime)
                {
                    StaticVariables.currentTime = StaticVariables.currentTime + StaticVariables.MDTimestep;
                    PhysicsEngine.VelocityVerlet();
                    iTime++;
                    writeFlag--;

                    /*
                    if (writeFlag == 0)
                    {
                        writeFlag = 100;
                        for (int i = 0; i < Atom.AllAtoms.Count; i++)
                        {
                            file.WriteLine(Atom.AllAtoms[i].position[0] + "    " + Atom.AllAtoms[i].position[1] + "    " + Atom.AllAtoms[i].position[2]);
                        }
                    }
                    */

                    Console.WriteLine("iTime = " + iTime + "            Current Time = " + StaticVariables.currentTime);

                }
            }

            Console.ReadLine();
		}
	}
}
