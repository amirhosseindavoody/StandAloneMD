using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace StandAloneMD
{
    public abstract class Potential
    {
        //this varaible keeps track of the current potential that is being used. (Note: only Lennard-Jones is currently implemented)
        public static potentialType currentPotential = potentialType.Buckingham;

        //Types of potential in the simulation
        public enum potentialType
        {
            LennardJones,
            Buckingham,
			REBO
        };

        public static Potential myPotential;

        abstract public void preCompute();
        abstract public void getForce(Atom firstAtom, Atom secondAtom);
        abstract public float getPotential(Atom firstAtom, Atom secondAtom);
        abstract public void calculateVerletRadius();
        abstract public void calculateNeighborList();
    }
}
