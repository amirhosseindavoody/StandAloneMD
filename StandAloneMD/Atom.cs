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


using System.Collections;
using System.Collections.Generic;
using System;

namespace StandAloneMD
{
	public abstract class Atom
	{

    	//this is a list of all atoms
		protected static List<Atom> m_AllAtoms = new List<Atom> ();

		//variables that must be implemented because they are declared as abstract in the base class
		public abstract float epsilon{ get; } // J
		public abstract float sigma { get; }
		protected abstract float massamu{ get; } //amu
		public abstract String atomName { get; }
		public abstract int atomID { get;}

		public abstract float buck_A { get; } // Buckingham potential coefficient
		public abstract float buck_B { get; } // Buckingham potential coefficient
		public abstract float buck_C { get; } // Buckingham potential coefficient
		public abstract float buck_D { get; } // Buckingham potential coefficient
		public abstract float Q_eff { get; } // Ion effective charge for use in Buckingham potential

		public float[] velocity = new float[3]{0.0f, 0.0f, 0.0f};
		public float[] position = new float[3]{0.0f, 0.0f, 0.0f};
		public float[] acceleration = new float[3] {0.0f, 0.0f, 0.0f};

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

    /*
    
	//the function returns the Lennard-Jones force on the atom given the list of all the atoms in the simulation
	Vector3 GetLennardJonesForce(List<Atom> objectsInRange){
		Vector3 finalForce = Vector3.zero;

		for (int i = 0; i < objectsInRange.Count; i++) {
			Vector3 deltaR = transform.position - objectsInRange [i].transform.position;
			float distanceSqr = deltaR.sqrMagnitude;

			//only get the forces of the atoms that are within the cutoff range
			if (objectsInRange[i].gameObject != gameObject && (distanceSqr < StaticVariables.cutoffSqr)) {

				float finalSigma = StaticVariables.sigmaValues[atomID * objectsInRange[i].atomID];

				int iR = (int) ((Mathf.Sqrt(distanceSqr)/finalSigma)/(StaticVariables.deltaR/StaticVariables.sigmaValueMax))+2;
				float magnitude = StaticVariables.preLennardJones[iR];
				magnitude = magnitude * 48.0f * epsilon / StaticVariables.angstromsToMeters/ finalSigma / finalSigma;
				finalForce += deltaR * magnitude;
			}
		}
		
		Vector3 adjustedForce = finalForce / StaticVariables.mass100amuToKg;
		adjustedForce = adjustedForce / StaticVariables.angstromsToMeters;
		adjustedForce = adjustedForce * StaticVariables.fixedUpdateIntervalToRealTime * StaticVariables.fixedUpdateIntervalToRealTime;
		return adjustedForce;
	}

	//the function returns the Buckingham force on the atom given the list of all the atoms in the simulation
	Vector3 GetBuckinghamForce(List<Atom> objectsInRange){
		Vector3 finalForce = Vector3.zero;
		
		for (int i = 0; i < objectsInRange.Count; i++) {
			Vector3 deltaR = transform.position - objectsInRange [i].transform.position;
			float distanceSqr = deltaR.sqrMagnitude;
			float magnitude = 0.0f;
			
			//only get the forces of the atoms that are within the cutoff range
			if ((objectsInRange[i].gameObject != gameObject) && (distanceSqr < StaticVariables.cutoffSqr)) {
				float distance = Mathf.Sqrt(distanceSqr);

				float final_A = StaticVariables.coeff_A [atomName + objectsInRange[i].atomName];
				float final_B = StaticVariables.coeff_B [atomName + objectsInRange[i].atomName];
				float final_C = StaticVariables.coeff_C [atomName + objectsInRange[i].atomName];
				float final_D = StaticVariables.coeff_D [atomName + objectsInRange[i].atomName];

				magnitude = final_A * final_B * Mathf.Exp (-final_B * distance) / distance;
				magnitude = magnitude - 6.0f * final_C / Mathf.Pow (distanceSqr, 4);
				magnitude = magnitude - 8.0f * final_D / Mathf.Pow (distanceSqr, 5);
				magnitude = magnitude / StaticVariables.angstromsToMeters;
				magnitude = magnitude + Q_eff * objectsInRange[i].Q_eff / (4.0f * Mathf.PI * StaticVariables.epsilon0 * distanceSqr * distance * StaticVariables.angstromsToMeters * StaticVariables.angstromsToMeters);

				finalForce += deltaR * magnitude;
			}
								
		}
		
		Vector3 adjustedForce = finalForce / StaticVariables.mass100amuToKg;
		adjustedForce = adjustedForce / StaticVariables.angstromsToMeters;
		adjustedForce = adjustedForce * StaticVariables.fixedUpdateIntervalToRealTime * StaticVariables.fixedUpdateIntervalToRealTime;
		return adjustedForce;
	}

    */
	}
}

