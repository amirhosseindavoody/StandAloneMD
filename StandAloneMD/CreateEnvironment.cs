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
		public int numAtoms = 10;
		public float width = 10.0f;
		public float height = 10.0f;
		public float depth = 10.0f;
		public float volume;
		private Random rnd = new Random();

		public void PreCompute ()
		{
            Atom.templateAtoms.Add(new Copper());
            Atom.templateAtoms.Add(new Gold());
            Atom.templateAtoms.Add(new Platinum());

            // delete the template atoms from the list of m_AllAtoms in Atom.cs
            for (int i = 0; i < Atom.templateAtoms.Count; i++)
            {
                Atom.AllAtoms.Clear();
            }

            LennardJones.preLennardJones();
            Buckingham.preBuckingham();
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
                float minDistance = 3.0f;
                Atom otherAtom = Atom.AllAtoms[i];
                float[] deltaR = new float[3];
                for (int idx = 0; idx < 3; idx++)
                {
                   deltaR[idx] = currAtom.position[idx] - otherAtom.position[idx];
                }
                float distanceSqr = deltaR[0] * deltaR[0] + deltaR[1] * deltaR[1] + deltaR[2] * deltaR[2];
                
                if (distanceSqr < (minDistance * minDistance))
                    proximityFlag = false;
            }
            return proximityFlag;
        }

	}
}
