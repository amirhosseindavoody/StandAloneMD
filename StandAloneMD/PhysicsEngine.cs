using System;
namespace StandAloneMD
{
	public class PhysicsEngine
	{

		public static void VelocityVerlet()
		{
			// update the position of all atoms then initialize the acceleration to be updated
			for (int i=0; i< Atom.AllAtoms.Count; i++)
            {
				Atom currAtom = Atom.AllAtoms[i];
                for (int idx = 0; idx < 3; idx++)
                {
                    currAtom.position[idx] = currAtom.position[idx] + currAtom.velocity [idx] * StaticVariables.MDTimestep + 0.5f * StaticVariables.MDTimestepSqr * currAtom.accelerationNew [idx];
                }
                currAtom.accelerationOld = currAtom.accelerationNew;
                currAtom.accelerationNew = new float[3] { 0.0f, 0.0f, 0.0f };
			}

            if (StaticVariables.currentPotential == StaticVariables.Potential.LennardJones)
            {
                if (StaticVariables.iTime % StaticVariables.nVerlet == 0)
                {
                    LennardJones.calculateNeighborList();
                }
                // update the acceleration of all atoms
			    for (int i=0; i< Atom.AllAtoms.Count-1; i++) {
				    Atom firstAtom = Atom.AllAtoms[i];
                    for (int j=0; j<firstAtom.neighborList.Count; j++) {
				    	Atom secondAtom = firstAtom.neighborList[j];
                        LennardJones.getForce(firstAtom, secondAtom);
				    }
			    }
            } 
            else if (StaticVariables.currentPotential == StaticVariables.Potential.Buckingham)
            {
                if (StaticVariables.iTime % StaticVariables.nVerlet == 0)
                {
                    Buckingham.calculateNeighborList();
                }
                // update the acceleration of all atoms
			    for (int i=0; i< Atom.AllAtoms.Count-1; i++) {
				    Atom firstAtom = Atom.AllAtoms[i];
                    for (int j=0; j<firstAtom.neighborList.Count; j++) {
				    	Atom secondAtom = firstAtom.neighborList[j];
                        Buckingham.getForce(firstAtom, secondAtom);
				    }
			    }
            }

            // update the velocity of all atoms
            for (int i = 0; i < Atom.AllAtoms.Count; i++)
            {
                Atom currAtom = Atom.AllAtoms[i];
                for (int idx = 0; idx < 3; idx++)
                {
                    currAtom.velocity[idx] = currAtom.velocity[idx] + 0.5f * (currAtom.accelerationOld[idx] + currAtom.accelerationNew[idx]) * StaticVariables.MDTimestep;
                    currAtom.velocity[idx] = currAtom.velocity[idx] * StaticVariables.sqrtAlpha;
                }
            }
		}

		

        //reflect the atoms from the walls
        public static void ReflectFromWalls()
        {
            float[] boxDimension = new float[3] {StaticVariables.myEnvironment.depth, StaticVariables.myEnvironment.width , StaticVariables.myEnvironment.height};
            
            for (int i = 0; i < Atom.AllAtoms.Count; i++)
            {
                Atom currAtom = Atom.AllAtoms[i];
                for (int idx = 0; idx < 3; idx++)
                {
                    float sign = Math.Sign(currAtom.position[idx]);
                    float firstRemainder = ((Math.Abs(currAtom.position[idx]) + boxDimension[idx]/2.0f) % (2.0f * boxDimension[idx]));
                    if (firstRemainder < boxDimension[idx])
                    {
                        currAtom.position[idx] = firstRemainder - boxDimension[idx] / 2.0f;
                        currAtom.velocity[idx] = currAtom.velocity[idx];
                    }
                    else
                    {
                        currAtom.position[idx] = 3.0f * boxDimension[idx] / 2.0f - firstRemainder;
                        currAtom.velocity[idx] = -1.0f * currAtom.velocity[idx];
                    }
                    currAtom.position[idx] = sign * currAtom.position[idx];
                }
                
            }
        }

        public static void CalculateEnergy()
        {
            StaticVariables.potentialEnergy = 0.0f;
            StaticVariables.kineticEnergy = 0.0f;
            StaticVariables.currentTemperature = 0.0f;

            for (int i = 0; i < Atom.AllAtoms.Count - 1; i++)
            {
                Atom firstAtom = Atom.AllAtoms[i];
                
                // calculate kinetic energy of each atom
                float velocitySqr = firstAtom.velocity[0] * firstAtom.velocity[0] + firstAtom.velocity[1] * firstAtom.velocity[1] + firstAtom.velocity[2] * firstAtom.velocity[2];
                StaticVariables.kineticEnergy += 0.5f * firstAtom.massamu * StaticVariables.amuToKg * velocitySqr * StaticVariables.angstromsToMeters * StaticVariables.angstromsToMeters;

                // calculate potential energy between each pair of atoms
                if (StaticVariables.currentPotential == StaticVariables.Potential.LennardJones)
                {
                    for (int j = 0; j < firstAtom.neighborList.Count; j++)
                    {
                        Atom secondAtom = firstAtom.neighborList[j];
                        StaticVariables.potentialEnergy += LennardJones.getPotential(firstAtom, secondAtom);
                    }
                }
                else if (StaticVariables.currentPotential == StaticVariables.Potential.Buckingham)
                {
                    for (int j = 0; j < firstAtom.neighborList.Count; j++)
                    {
                        Atom secondAtom = firstAtom.neighborList[j];
                        StaticVariables.potentialEnergy += Buckingham.getPotential(firstAtom, secondAtom);
                    }
                }
            }

            StaticVariables.currentTemperature = StaticVariables.kineticEnergy / 1.5f / (float)Atom.AllAtoms.Count / StaticVariables.kB;
            calculateSqrtAlpha();

        }

        //this function calculates the coefficient that scales the velocity of atoms to conserve the temperature
        private static void calculateSqrtAlpha()
        {
            float alpha = StaticVariables.desiredTemperature / StaticVariables.currentTemperature;
            float draggedAlpha = 0.0f;
            float draggedTemperature = 0.0f;

            if (StaticVariables.currentTemperature < 0.000000000001f)
            {
                draggedAlpha = 1.0f;
            }
            else if (StaticVariables.currentTemperature > 5000.0f)
            {
                draggedAlpha = alpha;
            }
            else if (alpha > 1)
            {
                draggedTemperature = (StaticVariables.desiredTemperature - StaticVariables.currentTemperature) * StaticVariables.alphaDrag + StaticVariables.currentTemperature;
                draggedAlpha = draggedTemperature / StaticVariables.currentTemperature;
            }
            else if (alpha < 1)
            {
                draggedTemperature = StaticVariables.currentTemperature - ((StaticVariables.currentTemperature - StaticVariables.desiredTemperature) * StaticVariables.alphaDrag);
                draggedAlpha = draggedTemperature / StaticVariables.currentTemperature;
            }
            else
            {
                draggedAlpha = 1.0f;
            }
            StaticVariables.sqrtAlpha = (float)Math.Pow(draggedAlpha, 0.5f);
        }
	}
}

