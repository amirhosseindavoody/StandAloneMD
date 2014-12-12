using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace StandAloneMD
{
    class PeriodicBoundary : BoundaryCondition
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
                    float sign = (float)Math.Sign(currAtom.position[idx]);
                    currAtom.position[idx] = sign * (((Math.Abs(currAtom.position[idx]) + boxDimension[idx] / 2.0f) % (boxDimension[idx])) - boxDimension[idx] / 2.0f);
                }
            }
        }

        public override float[] deltaPosition(Atom firstAtom, Atom secondAtom)
        {
            float[] boxDimension = new float[3] { CreateEnvironment.myEnvironment.depth, CreateEnvironment.myEnvironment.width, CreateEnvironment.myEnvironment.height };

            float[] deltaR = new float[3];
            for (int idx = 0; idx < 3; idx++)
            {
                deltaR[idx] = firstAtom.position[idx] - secondAtom.position[idx];
                if (Math.Abs(deltaR[idx]) > boxDimension[idx] / 2.0f)
                {
                    float sign = (float)Math.Sign(deltaR[idx]);
                    deltaR[idx] = deltaR[idx] - sign * boxDimension[idx];
                }
            }
            return deltaR;
        }
    }
}
