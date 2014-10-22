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
			get { return 196.967f; } //amu
		}

		// We assume gold to play the role of chloride
		public override float buck_A 
		{
			get { return 405774.0f*1.6f* (float)Math.Pow(10,-19); } //units of [J]
		}
	
		public override float buck_B 
		{
			get { return 4.207408f; } //units of [1/Angstrom]
		}
	
		public override float buck_C 
		{
			get { return 72.4f*1.6f* (float)Math.Pow(10,-19); } //units of [J.Anstrom^6]
		}
	
		public override float buck_D 
		{
			get { return 145.425f*1.6f* (float)Math.Pow(10,-19); } //units of [J.Angstrom^8]
		}
	
		public override float Q_eff 
		{
			get { return -1.0f*1.6f* (float)Math.Pow(10,-19); } //units of Coulomb
		}
	}
}

