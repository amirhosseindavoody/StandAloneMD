/**
 * Class: Atom.cs
 * Created by: Amirhossein Davoody
 **/

//L-J potentials from Zhen and Davies, Phys. Stat. Sol. a, 78, 595 (1983)
//Symbol, epsilon/k_Boltzmann (K) n-m version, 12-6 version, sigma (Angstroms),
//     mass in amu, mass in (20 amu) for Unity 
//     FCC lattice parameter in Angstroms, expected NN bond (Angs)
//Au: 4683.0, 5152.9, 2.6367, 196.967, 9.848, 4.080, 2.88
//Cu: 3401.1, 4733.5, 2.3374,  63.546, 3.177, 3.610, 2.55
//Pt: 7184.2, 7908.7, 2.5394, 165.084, 8.254, 3.920, 2.77


using System;
using System.Collections;
using System.Collections.Generic;

namespace StandAloneMD
{
	public abstract class Atom
	{

    	//this is a list of all atoms
		protected static List<Atom> m_AllAtoms = new List<Atom> ();

		//variables that must be implemented because they are declared as abstract in the base class
		public abstract float epsilon{ get; } // J
		public abstract float sigma { get; }
		public abstract float massamu{ get; } //amu
		public abstract String atomName { get; }
		public abstract int atomID { get;}

		public abstract float buck_A { get; } // Buckingham potential coefficient
		public abstract float buck_B { get; } // Buckingham potential coefficient
		public abstract float buck_C { get; } // Buckingham potential coefficient
		public abstract float buck_D { get; } // Buckingham potential coefficient
		public abstract float Q_eff { get; } // Ion effective charge for use in Buckingham potential

        public float verletRadius = 0.0f;
        public List<Atom> neighborList = new List<Atom>();

		public float[] velocity = new float[3]{0.0f, 0.0f, 0.0f};
		public float[] position = new float[3]{0.0f, 0.0f, 0.0f};
		public float[] accelerationNew = new float[3] {0.0f, 0.0f, 0.0f};
        public float[] accelerationOld = new float[3] { 0.0f, 0.0f, 0.0f };

		public Atom()
		{
			m_AllAtoms.Add(this);
		}

		// method to extract the list of allMolecules
		public static List<Atom> AllAtoms { 
			get
			{
				return m_AllAtoms;
			}
		}
	}
}

