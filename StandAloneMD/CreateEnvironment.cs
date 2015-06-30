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
		public float width;
		public float height;
		public float depth;
		public float volume;
		private Random rnd = new Random();

        public static CreateEnvironment myEnvironment;

		public void PreCompute ()
		{
            Atom.templateAtoms.Add(new Copper());
            Atom.templateAtoms.Add(new Gold());
            Atom.templateAtoms.Add(new Platinum());
            Atom.templateAtoms.Add(new Sodium());
            Atom.templateAtoms.Add(new Chlorine());
			Atom.templateAtoms.Add(new Hydrogen());
			Atom.templateAtoms.Add(new Carbon());

            // delete the template atoms from the list of m_AllAtoms in Atom.cs
            for (int i = 0; i < Atom.templateAtoms.Count; i++)
            {
                Atom.AllAtoms.Clear();
            }

            switch (Potential.currentPotential)
            {
                case Potential.potentialType.LennardJones:
                    Potential.myPotential = new LennardJones();
                    break;
                case Potential.potentialType.Buckingham:
                    Potential.myPotential = new Buckingham();
                    break;
				case Potential.potentialType.REBO:
					Potential.myPotential = new REBO();
					break;
            }

            Potential.myPotential.preCompute();
		}


	
		//initialize the atoms to a random position and to the original number of atoms
        public void InitAtoms()
		{
            //set the values for the initialization of atoms, this will later change to the real box size
            width = 5.0f;
            depth = width;
            height = width;
            volume = width * depth * height;

            if (Potential.currentPotential == Potential.potentialType.LennardJones)
            {
                for (int i = 0; i < numAtoms; i++)
                {
                    Atom currAtom = new Copper();
                    bool proximityFlag = false;

                    while (proximityFlag == false)
                    {
                        float maxInitialVelocity = 3.0f * (float)Math.Pow(10, 12);
                        currAtom.velocity = new float[] { randomFloat(-1.0f * maxInitialVelocity, +1.0f * maxInitialVelocity), randomFloat(-1.0f * maxInitialVelocity, +1.0f * maxInitialVelocity), randomFloat(-1.0f * maxInitialVelocity, +1.0f * maxInitialVelocity) };
                        currAtom.position = new float[] { randomFloat(-depth / 2.0f, depth / 2.0f), randomFloat(-width / 2.0f, width / 2.0f), randomFloat(-height / 2.0f, height / 2.0f) };
                        proximityFlag = checkProximity(currAtom);
                    }
                }
            }
            else if ((Potential.currentPotential == Potential.potentialType.Buckingham))
            {
                for (int i = 0; i < numAtoms/2; i++)
                {
                    Atom currAtom = new Sodium();
                    bool proximityFlag = false;

                    while (proximityFlag == false)
                    {
                        currAtom.velocity = new float[] { 0.0f, 0.0f, 0.0f };
                        currAtom.position = new float[] { randomFloat(-depth / 2.0f, depth / 2.0f), randomFloat(-width / 2.0f, width / 2.0f), randomFloat(-height / 2.0f, height / 2.0f) };
                        proximityFlag = checkProximity(currAtom);
                    }
                }
                for (int i = numAtoms/2; i < numAtoms; i++)
                {
                    Atom currAtom = new Chlorine();
                    bool proximityFlag = false;

                    while (proximityFlag == false)
                    {
                        currAtom.velocity = new float[] { 0.0f, 0.0f, 0.0f };
                        currAtom.position = new float[] { randomFloat(-depth / 2.0f, depth / 2.0f), randomFloat(-width / 2.0f, width / 2.0f), randomFloat(-height / 2.0f, height / 2.0f) };
                        proximityFlag = checkProximity(currAtom);
                    }
                }
            }
			else if ((Potential.currentPotential == Potential.potentialType.REBO))
			{
				//for (int i = 0; i < numAtoms / 2; i++)
				//{
				//	Atom currAtom = new Hydrogen();
				//	bool proximityFlag = false;

				//	while (proximityFlag == false)
				//	{
				//		currAtom.velocity = new float[] { 0.0f, 0.0f, 0.0f };
				//		currAtom.position = new float[] { randomFloat(-depth / 2.0f, depth / 2.0f), randomFloat(-width / 2.0f, width / 2.0f), randomFloat(-height / 2.0f, height / 2.0f) };
				//		proximityFlag = checkProximity(currAtom);
				//	}
				//}
				//for (int i = numAtoms / 2; i < numAtoms; i++)
				//{
				//	Atom currAtom = new Carbon();
				//	bool proximityFlag = false;

				//	while (proximityFlag == false)
				//	{
				//		currAtom.velocity = new float[] { 0.0f, 0.0f, 0.0f };
				//		currAtom.position = new float[] { randomFloat(-depth / 2.0f, depth / 2.0f), randomFloat(-width / 2.0f, width / 2.0f), randomFloat(-height / 2.0f, height / 2.0f) };
				//		proximityFlag = checkProximity(currAtom);
				//	}
				//}

				for (int i = 0; i < numAtoms; i++)
				{
					Atom currAtom = new Carbon();
					bool proximityFlag = false;

					while (proximityFlag == false)
					{
						currAtom.velocity = new float[] { 0.0f, 0.0f, 0.0f };
						currAtom.position = new float[] { randomFloat(-depth / 2.0f, depth / 2.0f), randomFloat(-width / 2.0f, width / 2.0f), randomFloat(-height / 2.0f, height / 2.0f) };
						proximityFlag = checkProximity(currAtom);
					}
				}
			}

            //now, set the values for the real box size which will be used to reflect atoms from walls.
            width = 10.0f;
            depth = width;
            height = width;
            volume = width * depth * height;

		}
        

        //initialize two atoms to fixed positions for debuging perposes.
        public void InitAtomsDebug()
        {
            Atom copperAtom = new Copper();
            copperAtom.velocity = new float[] { 0.0f, 0.0f, 0.0f };
            copperAtom.position = new float[] { -2.0f, 0.0f, 0.0f };

            Atom goldAtom = new Gold();
            goldAtom.velocity = new float[] { 0.0f, 0.0f, 0.0f };
            goldAtom.position = new float[] { +2.0f, 0.0f, 0.0f };

            width = 20.0f;
            depth = width;
            height = width;
            volume = width * depth * height;

            myEnvironment.numAtoms = Atom.AllAtoms.Count;
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
                //float minDistance = 2.35f;
				float minDistance = 1.0f;
                Atom otherAtom = Atom.AllAtoms[i];
                float[] deltaR = BoundaryCondition.myBoundary.deltaPosition(currAtom, otherAtom);
                float distanceSqr = deltaR[0] * deltaR[0] + deltaR[1] * deltaR[1] + deltaR[2] * deltaR[2];
                
                if (distanceSqr < (minDistance * minDistance))
                    proximityFlag = false;
            }
            return proximityFlag;
        }

	}
}
