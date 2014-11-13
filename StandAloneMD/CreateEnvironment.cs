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

using System;
using System.Collections;
using System.Collections.Generic;

namespace StandAloneMD
{
	public class CreateEnvironment
	{
		public int numAtoms = 100;
    	private int numAtomTypes = 3;
		public float width = 30.0f;
		public float height = 30.0f;
		public float depth = 30.0f;
		public float volume;
		private Random rnd = new Random();

		public void PreCompute ()
		{
			//the min sigma value and max sigma value are used for precalculating LJ forces.
			Atom[] atomTemplate = new Atom[numAtomTypes];
        	atomTemplate[0] = new Copper();
			atomTemplate[1] = new Gold();
			atomTemplate[2] = new Platinum();

			//precompute all of the coefficients for potentials in the system so it doesnt have to be done dynamically
			for (int i = 0; i < numAtomTypes; i++)
			{
				Atom firstAtom = atomTemplate[i];
				for(int j = 0; j < numAtomTypes; j++)
				{
					Atom secondAtom = atomTemplate[j];

					float currentSigma = (float)Math.Sqrt(firstAtom.sigma*secondAtom.sigma);
					StaticVariables.sigmaValues[firstAtom.atomID,secondAtom.atomID] = currentSigma;

                    // when the pre-calculated normalized Lennard Jones force is multiplied by this coefficient the acceleration units is [Angstrom/second^2]
                    float currentAccelCoeff = 24.0f * firstAtom.epsilon / (currentSigma * currentSigma * StaticVariables.angstromsToMeters * StaticVariables.angstromsToMeters * firstAtom.massamu * StaticVariables.amuToKg);
                    StaticVariables.accelCoefficient[firstAtom.atomID, secondAtom.atomID] = currentAccelCoeff;

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

			// delete the template atoms from the list of m_AllAtoms in Atom.cs
			for (int i=0; i<numAtomTypes; i++)
			{
				Atom.AllAtoms.Remove(atomTemplate[i]);
			}

			// precalculate the LennardJones potential and store it in preLennarJones array.
			int nR = (int)(StaticVariables.cutoff/StaticVariables.deltaR)+1;
			StaticVariables.preLennardJonesForce = new float[nR];
            StaticVariables.preLennardJonesPotential = new float[nR];

			for (int i = 0; i < nR; i++)
			{
				float distance = (float)i * StaticVariables.deltaR;
                StaticVariables.preLennardJonesForce[i] = calcLennardJonesForce(distance);
                StaticVariables.preLennardJonesPotential[i] = calcLennardJonesPotential(distance);
			}

		}

    	//the function returns the LennarJones force on the atom given the list of the atoms that are within range of it
    	private float calcLennardJonesForce(float distance)
    	{
        	float invDistance2 = 1.0f / distance / distance;
        	float invDistance6 = invDistance2 * invDistance2 * invDistance2;
            float invCutoff2 = 1.0f / StaticVariables.cutoff / StaticVariables.cutoff;
            float invCutoff6 = invCutoff2 * invCutoff2 * invCutoff2;
        	float r_min = StaticVariables.rMinMultiplier;

        	float forceMagnitude = 0.0f;

        	if (distance > r_min)
        	{
                forceMagnitude = invDistance2 * ((2.0f * invDistance6 * invDistance6 - invDistance6) - (invCutoff2 / invDistance2) * (2.0f * invCutoff6 * invCutoff6 - invCutoff6 ));
        	}
        	// Smooth the potential to go to a constant not infinity at r=0
        	else
        	{
            	float invr_min = 1 / r_min;
            	float invr_min2 = invr_min * invr_min;
            	float invr_min6 = invr_min2 * invr_min2 * invr_min2;
                float magnitude_Vmin = invr_min2 * ((2.0f * invr_min6 * invr_min6 - invr_min6) - (invCutoff2 / invr_min2) * (2.0f * invCutoff6 * invCutoff6 - invCutoff6));

            	float r_Vmax = r_min / 1.5f;
            	float invr_Vmax2 = 1 / r_Vmax / r_Vmax;
            	float invr_Vmax6 = invr_Vmax2 * invr_Vmax2 * invr_Vmax2;
                float magnitude_Vmax = invr_Vmax2 * ((2.0f * invr_Vmax6 * invr_Vmax6 - invr_Vmax6) - (invCutoff2 / invr_Vmax2) * (2.0f * invCutoff6 * invCutoff6 - invCutoff6));

            	float part1 = (distance / r_min) * ((float)Math.Exp(distance - r_min));
            	float part2 = magnitude_Vmax - magnitude_Vmin;
            	forceMagnitude = magnitude_Vmax - (part1 * part2);
        	}
            
        	return forceMagnitude;
    	}

        //the function returns the LennarJones force on the atom given the list of the atoms that are within range of it
        private float calcLennardJonesPotential(float distance)
        {
            float invDistance2 = 1.0f / distance / distance;
            float invDistance6 = invDistance2 * invDistance2 * invDistance2;
            float invCutoff2 = 1.0f / StaticVariables.cutoff / StaticVariables.cutoff;
            float invCutoff6 = invCutoff2 * invCutoff2 * invCutoff2;
            float r_min = StaticVariables.rMinMultiplier;

            float potential = 0.0f;

            if (distance > 0.0f)
            {
                potential = 4.0f * ((invDistance6 * invDistance6 - invDistance6) + (6.0f * invCutoff6 * invCutoff6 - 3.0f * invCutoff6) * (invCutoff2 / invDistance2) - 7.0f * invCutoff6 * invCutoff6 + 4.0f * invCutoff6);
            }
            
            return potential;
        }


	
		//initialize the atoms to a random position and to the original number of atoms
		public void InitAtoms()
		{
			for (int i = 0; i < numAtoms; i++)
			{
				Atom currAtom = new Copper();
                bool proximityFlag = false;
 
                while (proximityFlag == false)
                {
                    float maxInitialVelocity = 3.0f * (float)Math.Pow(10,12);
                    currAtom.velocity = new float[] { randomFloat(-1.0f * maxInitialVelocity, +1.0f * maxInitialVelocity), randomFloat(-1.0f * maxInitialVelocity, +1.0f * maxInitialVelocity), randomFloat(-1.0f * maxInitialVelocity, +1.0f * maxInitialVelocity) };
                    //currAtom.position = new float[] { randomFloat(-depth / 2.0f, depth / 2.0f), randomFloat(-width / 2.0f, width / 2.0f), 0.0f };
                    currAtom.position = new float[] { randomFloat(-depth / 2.0f, depth / 2.0f), randomFloat(-width / 2.0f, width / 2.0f), randomFloat(-height / 2.0f, height / 2.0f) };
                    proximityFlag = checkProximity(currAtom);
                    
                }
			}	
		}

        //this method generates a random float number between the minValue and maxValue
        private float randomFloat(float minValue, float maxValue)
    	{
        	float rndFloat = (maxValue - minValue) * (float)rnd.NextDouble() + minValue;
        	return rndFloat;
    	}

        //check the distance between the atoms, if it is larger than the equilibrium position move accept the random number, otherwise pick another set of random positions.
        private bool checkProximity(Atom currAtom)
        {
            bool proximityFlag = true;
            for (int i = 0; i < Atom.AllAtoms.Count - 1; i++)
            {
                Atom otherAtom = Atom.AllAtoms[i];
                float[] deltaR = new float[3];
                for (int idx = 0; idx < 3; idx++)
                {
                   deltaR[idx] = currAtom.position[idx] - otherAtom.position[idx];
                }
                float distanceSqr = deltaR[0] * deltaR[0] + deltaR[1] * deltaR[1] + deltaR[2] * deltaR[2];
                float finalSigma = StaticVariables.sigmaValues[currAtom.atomID, otherAtom.atomID];
                float normDistanceSqr = distanceSqr / finalSigma / finalSigma; // this is normalized distanceSqr to the sigmaValue

                //only get the forces of the atoms that are within the cutoff range
                if (normDistanceSqr < 1.259921)
                {
                    proximityFlag = false;
                }
            }
            return proximityFlag;
        }

	}
}
