using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace StandAloneMD
{
    abstract class BoundaryCondition
    {
        //this varaible keeps track of the current potential that is being used. (Note: only Lennard-Jones is currently implemented)
        public static BoundaryCondition myBoundary = new ReflectingBoundary();

        
        //reflect the atoms from the walls
        public abstract void Apply();
        public abstract float[] deltaPosition(Atom firstAtom, Atom secondAtom);


        
    }
}
