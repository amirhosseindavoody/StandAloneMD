/**
 * Class: Platinum.cs
 * Created by: Justin Moeller
 * Description: This class defines anything that is platinum specific, and NOT related to all of
 * the atoms. This class is derived from the base class of Atom.cs, and takes on all of its behavior.
 * It must override all of the abstract variables and functions that are defined in Atom.cs, such as
 * atomName, epsilon, sigma, massamu, SetSelected(), and SetTransparent().
 * 
 * 
 **/

using System.Collections;
using System;

public class Platinum : Atom {
	
	public override String atomName { 
		get{ return "Platinum"; } 
	}

	public override int atomID {
		get{ return 2;}
	}
	
	public override float epsilon {
		get { return ((float)(1.0922 * Math.Pow(10, -19))); } // J
	}
		
	public override float sigma {
		get { return 2.5394f; }
	}
	
	protected override float massamu {
		get { return 195.084f; } //amu
	}

	// These are just dummy variables for platinum
	public override float buck_A {
		get { return 0.0f; } //units of [J]
	}
	
	public override float buck_B {
		get { return 0.0f; } //units of [1/Angstrom]
	}
	
	public override float buck_C {
		get { return 0.0f; } //units of [J.Anstrom^6]
	}
	
	public override float buck_D {
		get { return 0.0f; } //units of [J.Angstrom^8]
	}
	
	public override float Q_eff {
		get { return 0.0f; } //units of Coulomb
	}
}
