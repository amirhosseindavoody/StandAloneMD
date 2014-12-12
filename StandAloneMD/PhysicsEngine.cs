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

            /*
            if (StaticVariables.iTime % StaticVariables.nVerlet == 0)
            {
                Potential.myPotential.calculateNeighborList(distance, firstAtom, secondAtom);
            }
             */

            updateNeighborList();
            
            // update the acceleration of all atoms
            for (int i = 0; i < Atom.AllAtoms.Count - 1; i++)
            {
                Atom firstAtom = Atom.AllAtoms[i];
                for (int j = 0; j < firstAtom.neighborList.Count; j++)
                {
                    Atom secondAtom = firstAtom.neighborList[j];
                    Potential.myPotential.getForce(firstAtom, secondAtom);
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

        private static void updateNeighborList()
        {
            if (StaticVariables.iTime % StaticVariables.nVerlet == 0)
            {
                //clear the old neighborList
                for (int i = 0; i < Atom.AllAtoms.Count - 1; i++)
                {
                    Atom currAtom = Atom.AllAtoms[i];
                    currAtom.neighborList.Clear();
                }

                //create the new neighborList
                for (int i = 0; i < Atom.AllAtoms.Count - 1; i++)
                {
                    Atom firstAtom = Atom.AllAtoms[i];
                    for (int j = i + 1; j < Atom.AllAtoms.Count; j++)
                    {
                        Atom secondAtom = Atom.AllAtoms[j];
                        float[] deltaR = Boundary.deltaPosition(firstAtom, secondAtom);
                        float distance = (float)Math.Sqrt(deltaR[0] * deltaR[0] + deltaR[1] * deltaR[1] + deltaR[2] * deltaR[2]);
                        Potential.myPotential.calculateNeighborList(distance, firstAtom, secondAtom);
                        PairDistributionFunction.updatePairDistribution(distance);
                    }
                }
                PairDistributionFunction.calculateAveragePairDistribution();
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
                for (int j = 0; j < firstAtom.neighborList.Count; j++)
                {
                    Atom secondAtom = firstAtom.neighborList[j];
                    StaticVariables.potentialEnergy += Potential.myPotential.getPotential(firstAtom, secondAtom);
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

