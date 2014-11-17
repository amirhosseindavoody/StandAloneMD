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
        public static int iTime = 0;

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

        //Number of MD timesteps to update verlet list
        public static int nVerlet = 100;
		
		//Temperature slider bounds in K
		public static float tempRangeLow = 0.01f;
		public static float tempRangeHigh = 5000.0f;
        public static float desiredTemperature = 300.0f;

		//this varaible keeps track of the current potential that is being used. (Note: only Lennard-Jones is currently implemented)
		public static Potential currentPotential = Potential.Buckingham;
		
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
