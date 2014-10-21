/**
 * Class: CreateEnvironment.cs
 * Created by: Justin Moeller
 * Description: This class handles anything that has to do with the creation of the environment.
 * (i.e the atoms, the box, the lines, and the width/height/depth text). Most of the functionality 
 * of this class happens when the game begins, however, the box, the lines, and the text all must be 
 * scaled when the user changes the volume of the box. Additionally, to cut down on some computation,
 * all of the sigma values are pre-computed and then stored in a dictionary. The key to the dictionary
 * is a String made up of the two atoms for which you want the sigma value for. (i.e "CopperCopper" or
 * "CopperGold" or "GoldPlatinum") All of the atoms are also named (given a number 0-numAtoms) for easier 
 * access and debugging. 
 * 
 * 
 **/ 

using System.Collections;
using System.Collections.Generic;
using System;

public class CreateEnvironment {

	public int numAtoms = 10;
    private int numAtomTypes = 3;
	public float width = 100.0f;
	public float height = 100.0f;
	public float depth = 100.0f;
	public float volume;

	public void PreCompute () {

		//the min sigma value and max sigma value are used for precalculating LJ forces.
		Atom[] atomSample = new Atom[numAtomTypes];
        atomSample[0] = new Copper();
        atomSample[1] = new Gold();
        atomSample[2] = new Platinum();

		//precompute all of the coefficients for potentials in the system so it doesnt have to be done dynamically
		for (int i = 0; i < numAtomTypes; i++) {
			Atom firstAtom = atomSample[i];
			for(int j = 0; j < numAtomTypes; j++){
				Atom secondAtom = atomSample[j];

				float currentSigma = (float)Math.Sqrt(firstAtom.sigma*secondAtom.sigma);
				StaticVariables.sigmaValues[firstAtom.atomID,secondAtom.atomID] = currentSigma;

				float currentA = (float)Math.Sqrt(firstAtom.buck_A*secondAtom.buck_A);
				StaticVariables.coeff_A[firstAtom.atomID,secondAtom.atomID] = currentA;

                float currentB = (float)Math.Sqrt(firstAtom.buck_B * secondAtom.buck_B);
                StaticVariables.coeff_A[firstAtom.atomID, secondAtom.atomID] = currentB;

                float currentC = (float)Math.Sqrt(firstAtom.buck_C * secondAtom.buck_C);
                StaticVariables.coeff_A[firstAtom.atomID, secondAtom.atomID] = currentC;

                float currentD = (float)Math.Sqrt(firstAtom.buck_D * secondAtom.buck_D);
                StaticVariables.coeff_A[firstAtom.atomID, secondAtom.atomID] = currentD;
			}
		}

		// precalculate the LennardJones potential and store it in preLennarJones list.

		int nR = (int)(StaticVariables.cutoff/StaticVariables.deltaR)+1;

		StaticVariables.preLennardJones = new float[nR];

		for (int i = 0; i < nR; i++) {
			float distance = (float)i * StaticVariables.deltaR;
			float magnitude = CalcLennardJonesForce (distance);
			StaticVariables.preLennardJones[i] = magnitude;
		}

	}

    //the function returns the LennarJones force on the atom given the list of the atoms that are within range of it
    private float CalcLennardJonesForce(float distance)
    {
        float invDistance2 = 1.0f / distance / distance;
        float invDistance6 = invDistance2 * invDistance2 * invDistance2;
        float r_min = StaticVariables.rMinMultiplier;

        float forceMagnitude = 0.0f;

        if (distance > r_min)
        {
            forceMagnitude = invDistance2 * invDistance6 * (invDistance6 - 0.5f);
        }
        // Smooth the potential to go to a constant not infinity at r=0
        else
        {
            float invr_min = 1 / r_min;
            float invr_min2 = invr_min * invr_min;
            float invr_min6 = invr_min2 * invr_min2 * invr_min2;

            float magnitude_Vmin = invr_min2 * invr_min6 * (invr_min6 - 0.5f);

            float r_Vmax = r_min / 1.5f;
            float invr_Vmax2 = 1 / r_Vmax / r_Vmax;
            float invr_Vmax6 = invr_Vmax2 * invr_Vmax2 * invr_Vmax2;

            float magnitude_Vmax = invr_Vmax2 * invr_Vmax6 * (invr_Vmax6 - 0.5f);

            float part1 = (distance / r_min) * ((float)Math.Exp(distance - r_min));
            float part2 = magnitude_Vmax - magnitude_Vmin;
            forceMagnitude = magnitude_Vmax - (part1 * part2);
        }
        return forceMagnitude;
    }


	
	//initialize the atoms to a random position and to the original number of atoms
	public void InitAtoms(){
		
		for (int i = 0; i < numAtoms; i++) {
			Atom currAtom = new Copper();
            Atom.AllAtoms.Add(currAtom);
            currAtom.velocity = new float[] { randomFloat(-1.0f, +1.0f), randomFloat(-1.0f, +1.0f), randomFloat(-1.0f, +1.0f)};
            currAtom.position = new float[] { randomFloat(0.0f, depth), randomFloat(0.0f, width), randomFloat(0.0f, height) };
		}
	}

    private float randomFloat(float minValue, float maxValue)
    {
        Random rnd = new Random();
        float rndFloat = (maxValue - minValue) * (float)rnd.NextDouble() + minValue;
        return rndFloat;
    }

    public void printAtomList()
    {
        Console.WriteLine("Total number of created atoms = " + Atom.AllAtoms.Count);
        for (int i = 0; i < Atom.AllAtoms.Count; i++)
        {
            Atom currAtom = Atom.AllAtoms[i];
            Console.WriteLine(currAtom.atomName + "     " + currAtom.sigma);
        }
        Console.ReadLine();
    }

}
