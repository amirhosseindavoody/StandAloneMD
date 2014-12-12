using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace StandAloneMD
{
    class Boundary
    {
        //this varaible keeps track of the current potential that is being used. (Note: only Lennard-Jones is currently implemented)
        private static boundaryType currentBoundary = boundaryType.Box;

        //Types of potential in the simulation
        public enum boundaryType
        {
            Box,
            Periodic
        };
        //reflect the atoms from the walls
        public static void Box()
        {
            float[] boxDimension = new float[3] { CreateEnvironment.myEnvironment.depth, CreateEnvironment.myEnvironment.width, CreateEnvironment.myEnvironment.height };

            for (int i = 0; i < Atom.AllAtoms.Count; i++)
            {
                Atom currAtom = Atom.AllAtoms[i];
                for (int idx = 0; idx < 3; idx++)
                {
                    float sign = Math.Sign(currAtom.position[idx]);
                    float firstRemainder = ((Math.Abs(currAtom.position[idx]) + boxDimension[idx] / 2.0f) % (2.0f * boxDimension[idx]));
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

        //apply the boundary condition on atoms that are out of bound
        public static void Periodic()
        {
            float[] boxDimension = new float[3] { CreateEnvironment.myEnvironment.depth, CreateEnvironment.myEnvironment.width, CreateEnvironment.myEnvironment.height };

            for (int i = 0; i < Atom.AllAtoms.Count; i++)
            {
                Atom currAtom = Atom.AllAtoms[i];
                for (int idx = 0; idx < 3; idx++)
                {
                    float sign = (float)Math.Sign(currAtom.position[idx]);
                    currAtom.position[idx] = sign * (((Math.Abs(currAtom.position[idx]) + boxDimension[idx] / 2.0f) % (boxDimension[idx])) - boxDimension[idx] / 2.0f);
                }
            }
        }


        public static float[] deltaPosition(Atom firstAtom, Atom secondAtom)
        {
            float[] boxDimension = new float[3] { CreateEnvironment.myEnvironment.depth, CreateEnvironment.myEnvironment.width, CreateEnvironment.myEnvironment.height };

            float[] deltaR = new float[3];
            for (int idx = 0; idx < 3; idx++)
            {
                deltaR[idx] = firstAtom.position[idx] - secondAtom.position[idx];
                if ((Math.Abs(deltaR[idx]) > (boxDimension[idx] / 2.0f)) && (currentBoundary == boundaryType.Periodic))
                {
                    float sign = (float)Math.Sign(deltaR[idx]);
                    deltaR[idx] = deltaR[idx] - sign * boxDimension[idx];
                }
            }
            return deltaR;
        }
    }
}
