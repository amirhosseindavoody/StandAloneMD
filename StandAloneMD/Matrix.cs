using System;
using System.Collections;
using System.Collections.Generic;

namespace StandAloneMD
{
    public static class Matirx
    {
        // method to add two vectors with each other
        public static float[] Add (float[] vector1, float[] vector2)
        {
            float[] result = new float[vector1.Length];
            for (int i = 0; i < vector1.Length; i++)
            {
                result[i] = vector1[i] + vector2[i];
            }
            return result;
        }

        // method to subtract two vectors from each other
        public static float[] Subtract(float[] vector1, float[] vector2)
        {
            float[] result = new float[vector1.Length];
            for (int i = 0; i < vector1.Length; i++)
            {
                result[i] = vector1[i] - vector2[i];
            }
            return result;
        }

        // method to multiply a number to a vector
        public static float[] scalarMultiply(float scalar, float[] vector)
        {
            float[] result = new float[vector.Length];
            for (int i = 0; i < vector.Length; i++)
            {
                result[i] = scalar * vector[i];
            }
            return result;
        }

    }
}
