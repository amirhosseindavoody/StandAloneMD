using System;

namespace StandAloneMD
{
	class MainClass
	{
		public static void Main (string[] args)
		{
            CreateEnvironment environment = new CreateEnvironment();
            Console.WriteLine("Number of atoms = " + environment.numAtoms);

            environment.PreCompute();
            environment.InitAtoms();

		}
	}
}
