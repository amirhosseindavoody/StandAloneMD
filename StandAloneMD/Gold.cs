/**
 * Class: Gold.cs
 * Created by: Justin Moeller
 * Description: This class defines anything that is gold specific, and NOT related to all of
 * the atoms. This class is derived from the base class of Atom.cs, and takes on all of its behavior.
 * It must override all of the abstract variables and functions that are defined in Atom.cs, such as
 * atomName, epsilon, sigma, massamu, SetSelected(), and SetTransparent().
 * 
 * 
 **/

using System;
using System.Collections;


namespace StandAloneMD
{
	public class Gold : Atom
	{
		public override String atomName
		{ 
			get{ return "Gold"; } 
		}

		public override int atomID
		{
			get{ return 1;}
		}

		public override int atomicNumber
		{
			get { return 79; }
		}
	
		public override float epsilon 
		{
			get { return 5152.9f * 1.381f * (float)Math.Pow (10, -23); } // J
		}
	
		public override float sigma 
		{
			get { return 2.6367f; }
		}
	
		public override float massamu 
		{
			get { return 196.967f; } //amu for Gold
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

