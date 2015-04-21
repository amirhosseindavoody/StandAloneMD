using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace StandAloneMD
{
    class REBO : Potential
    {
        private float[,] REBO_firstneigh;
        private float[, ,] preBuckinghamAcceleration;
        
        public REBO()
        {

        }

        public override void preCompute()
        {
            
        }

        //the function returns the LennarJones force on the atom given the list of the atoms that are within range of it
        private float calcAcceleration(float distance, Atom firstAtom, Atom secondAtom)
        {
            
            return 0.0f;
        }

        //the function returns the LennarJones force on the atom given the list of the atoms that are within range of it
        private float calcPotential(float distance, Atom firstAtom, Atom secondAtom)
        {
            
            return 0.0f;
        }

        //the function returns the Lennard-Jones force on the atom given the list of all the atoms in the simulation
        public override void getForce(Atom firstAtom, Atom secondAtom)
        {
            
        }

        //the function returns the Lennard-Jones force on the atom given the list of all the atoms in the simulation
        public override float getPotential(Atom firstAtom, Atom secondAtom)
        {
            return 0.0f;
        }

        public override void calculateVerletRadius()
        {
            
        }

        //This function creates a list of all neighbor list for each atom
        public override void calculateNeighborList()
        {
            
        }


        private void REBO_neigh()
        {


        }



        /* ----------------------------------------------------------------------
        cutoff function Sprime
        return cutoff and dX = derivative
        no side effects
        ------------------------------------------------------------------------- */
        private float Sp(float Xij, float Xmin, float Xmax, ref float dX) 
        {
            float cutoff;

            float t = (Xij-Xmin) / (Xmax-Xmin);
            if (t <= 0.0f) 
            {
                cutoff = 1.0f;
                dX = 0.0f;
            } else if (t >= 1.0f) {
            cutoff = 0.0f;
            dX = 0.0f;
            } else {
                cutoff = 0.5f * (1.0f+(float)(Math.Cos(t*MathConst.MY_PI)));
                dX = (-0.5f*MathConst.MY_PI*(float)(Math.Sin(t*MathConst.MY_PI))) / (Xmax-Xmin);
            }
            return cutoff;
        }
    }
}
