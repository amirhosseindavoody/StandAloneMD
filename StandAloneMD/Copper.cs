/**
 * Class: Copper.cs
 * Created by: Justin Moeller
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
	public class Copper : Atom
	{
	
		public override String atomName 
		{ 
			get{ return "Copper"; } 
		}

		public override int atomID
		{
			get{ return 0;}
		}
	
		public override float epsilon
		{
			get { return ((float)(6.537 * Math.Pow(10, -20))); } // J
		}
		
		public override float sigma
		{
			get { return 2.3374f; }
		}
	
	    public override float massamu
		{
			//get { return 63.546f; } //amu
            get { return 22.98977f; } //amu for Sodium
		}

		// We assume copper to play the role of sodium
		public override float buck_A 
		{
			get { return 487.0f*1.6f* (float)Math.Pow(10,-19); } //units of [J]
		}

		public override float buck_B 
		{
			get { return 4.207408f; } //units of [1/Angstrom]
		}

		public override float buck_C
		{
			get { return 1.048f*1.6f* (float)Math.Pow(10,-19); } //units of [J.Anstrom^6]
		}

		public override float buck_D
		{
			get { return 0.499f*1.6f* (float)Math.Pow(10,-19); } //units of [J.Angstrom^8]
		}

		public override float Q_eff
		{
			get { return 1.0f*1.6f* (float)Math.Pow(10,-19); } //units of Coulomb
		}
	}
}
