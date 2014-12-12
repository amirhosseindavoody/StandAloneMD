using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace StandAloneMD
{
    class ReflectingBoundary : BoundaryCondition
    {
        //reflect the atoms from the walls
        public override void Apply()
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

        public override float[] deltaPosition(Atom firstAtom, Atom secondAtom)
        {
            float[] deltaR = new float[3];
            for (int idx = 0; idx < 3; idx++)
            {
                deltaR[idx] = firstAtom.position[idx] - secondAtom.position[idx];
            }
            return deltaR;
        }

    }
}
