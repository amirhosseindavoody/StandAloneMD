/**
 * Class: Hydrogen.cs
 * Created by: Amirhossein Davoody
 * Description: This class defines anything that is copper specific, and NOT related to all of
 * the atoms. This class is derived from the base class of Atom.cs, and takes on all of its behavior.
 * It must override all of the abstract variables and functions that are defined in Atom.cs, such as
 * atomName, epsilon, sigma, massamu, SetSelected(), and SetTransparent().
 * 
 * 
 **/ 

using System.Collections;
using System;

namespace StandAloneMD
{
	public class Hydrogen : Atom
	{
	
		public override String atomName 
		{ 
			get{ return "Hydrogen"; } 
		}

		public override int atomID
		{
			get{ return 5;}
		}
	
		public override float epsilon
		{
			get { return 0.0f; } // J
		}
		
		public override float sigma
		{
			get { return 0.0f; }
		}
	
	    public override float massamu
		{
            get { return 1.0f; } //amu for Hydrogen
		}

		public override float buck_A 
		{
			get { return 0.0f; } //units of [J]
		}

		public override float buck_B 
		{
			get { return 0.0f; } //units of [1/Angstrom]
		}

		public override float buck_C
		{
			get { return 0.0f; } //units of [J.Anstrom^6]
		}

		public override float buck_D
		{
			get { return 0.0f; } //units of [J.Angstrom^8]
		}

		public override float Q_eff
		{
			get { return 0.0f; } //units of Coulomb
		}
	}
}
