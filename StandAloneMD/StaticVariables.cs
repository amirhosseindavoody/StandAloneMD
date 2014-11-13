/**
 * Class: StaticVariables.cs
 * Created By: Justin Moeller
 * Description: This class is simply a list of static variables and static functions
 * that can be called from any class. The static variables in this list are either
 * constants or variables that can be controlled from across the entire system of
 * atom (i.e currentPotential). There are two functions in this class, DrawLine and
 * DrawQuad. DrawLine draw a line in 3D space and DrawQuad draws a quad in 3D space.
 * To use these functions in 2D space, the coordinates given to the function must be
 * translated such that they are rotated based on the rotation of the camera. 
 * 
 * 
 **/ 


using System.Collections;
using System;
using System.Collections.Generic;

namespace StandAloneMD
{
	public class StaticVariables 
	{
		//this is the molecular dynamic simulation timestep
    	public static float MDTimestep = 0.5f * (float) Math.Pow (10, -15);
		public static float MDTimestepSqr = MDTimestep * MDTimestep;
    	public static float fixedDeltaTime = MDTimestep;

    	//this variable keeps track of the amount of simulation time that has passed
    	public static float currentTime = 0.0f;

		//do not scale temperature all at once
		//public static float alphaDrag = 1.0f * fixedDeltaTime;
        public static float alphaDrag = 1.0f * 0.001f;

		//Boltzmann constant in J/K
		public static float kB = 1.381f * (float) Math.Pow(10,-23);

		//Permittivity of free space
		public static float epsilon0 = 8.85f * (float) Math.Pow (10, -12);

    	//Convert units of 1 amu to kg
    	public static float amuToKg = 1.6605f * (float)Math.Pow(10, -27); 

		//Convert units of 100 amu to kg
		public static float mass100amuToKg = 100f * amuToKg; 

		//Convert units of Angstroms to meters
		public static float angstromsToMeters = (float) Math.Pow (10,-10);
		
		//Cutoff distance for calculating LennarJones force. This quantity is unit less and normalized to sigmaValue for atom pair
		public static float cutoff = 2.5f; //[unit less]
		public static float cutoffSqr = cutoff * cutoff;

        //Number of MD timesteps to update verlet list
        public static int nVerlet = 100;

		//The mesh size for pre-calculating Lennard Jones force.
		public static float deltaR = 0.000001f; 

		//When r_ij is small, the Lennard-Jones potential is extremely large.
		//At a certain r_min, we will substitute the L-J potential with a function that
		//curves to a constant as r_ij goes to zero.

		//Multiplier for transition between actual L-J potential and curve to constant
		//    This number will be multiplied by sigma to find the transition distance
		public static float rMinMultiplier = 0.75f;
		
		//Temperature slider bounds in K
		public static float tempRangeLow = 0.01f;
		public static float tempRangeHigh = 5000.0f;

        public static float desiredTemperature = 300.0f;

		//pre-calculated coefficients and forces for Lennard-Jones potential
		public static float[,] sigmaValues = new float[3,3];
        public static float[,] accelCoefficient = new float[3, 3]; // this is the coefficient that is multiplied by the preLennardJones vector to get the acceleration of each atom for each combinations
        public static float[] preLennardJonesForce; //This is the pre-calculated value of LennardJones force for some mesh points.
        public static float[] preLennardJonesPotential; //This is the pre-calculated value of LennardJones potential for some mesh points.

        //pre-calculated coefficients and forces for Buckingham potential
		public static float[,] coeff_A = new float[3, 3];
    	public static float[,] coeff_B = new float[3, 3];
    	public static float[,] coeff_C = new float[3, 3];
    	public static float[,] coeff_D = new float[3, 3];

		//this varaible keeps track of the current potential that is being used. (Note: only Lennard-Jones is currently implemented)
		public static Potential currentPotential = Potential.LennardJones;
		
		//There are three potentials, but currently Lennard-Jones is the only one that is implemented so changing
		//between these potentials doesnt do anything
		public enum Potential
		{
			LennardJones,
			Brenner,
			Buckingham
		};

        public static CreateEnvironment myEnvironment;
        public static float kineticEnergy = 0.0f;  // units in Joules
        public static float potentialEnergy = 0.0f;  // units in Joules
        public static float currentTemperature = 0.0f;  // units in Kelvin

        public static float sqrtAlpha = 1.0f;
	}
}
