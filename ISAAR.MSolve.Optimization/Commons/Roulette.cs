using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Troschuetz.Random;

namespace ISAAR.MSolve.Optimization.Commons
{
    class Roulette
    {
        /// <summary>
        /// Contains the cummulative sizes of each pocket. These are normalized to belong to the interval [0, 1].
        /// No need to store the last one.
        /// </summary>
        private readonly double[] wheel;
        private readonly IGenerator rng;

        /// <summary>
        /// Creates a new instance of <see cref="Roulette"/>
        /// </summary>
        /// <param name="pocketSizes">The sizes of each pocket of the roulette wheel. 
        ///                           They must be normalized, namely positive with a sum of 1.0</param>
        /// <param name="randomNumberGenerator">The random number generator that this object will use</param>
        private Roulette(double[] pocketSizes, IGenerator randomNumberGenerator)
        {
            wheel = new double[pocketSizes.Length - 1];
            wheel[0] = pocketSizes[0];
            for (int i = 1; i < pocketSizes.Length - 1; ++i)
            {
                wheel[i] = wheel[i] + pocketSizes[i-1];
            }

            this.rng = randomNumberGenerator;
        }

        /// <summary>
        /// Creates a new instance of <see cref="Roulette"/> with the raw values passed in as pocket sizes. 
        /// Only use this factory method if you are sure the pocket sizes are positive and their sum = 1.0
        /// </summary>
        /// <param name="positiveNormalizedPocketSizes">The raw values of the pocket sizes</param>
        /// <param name="randomNumberGenerator">The random number generator that this object will use</param>
        /// <returns>A new instance of <see cref="Roulette"/></returns>
        public static Roulette CreateWithNormalized(double[] positiveNormalizedPocketSizes, IGenerator randomNumberGenerator)
        {
            return new Roulette(positiveNormalizedPocketSizes, randomNumberGenerator);
        }

        /// <summary>
        /// Creates a new <see cref="Roulette"/> after normalizing the raw values passed in as pocket sizes. 
        /// Only use this factory method if you are sure the pocket sizes are positive.
        /// </summary>
        /// <param name="positivePocketSizes">The raw values of the pocket sizes</param>
        /// <param name="randomNumberGenerator">The random number generator that this object will use</param>
        /// <returns>A new instance of <see cref="Roulette"/></returns>
        public static Roulette CreateFromPositive(double[] positivePocketSizes, IGenerator randomNumberGenerator)
        {
            Normalize(positivePocketSizes);
            return new Roulette(positivePocketSizes, randomNumberGenerator);
        }

        /// <summary>
        /// Creates a new <see cref="Roulette"/> after checking and normalizing the raw values passed in as pocket sizes.
        /// These values need to be positive, otherwise an <see cref="ArgumentException"/> will be thrown.
        /// </summary>
        /// <param name="pocketSizes">The raw values of the pocket sizes. They must be positive.</param>
        /// <param name="randomNumberGenerator">The random number generator that this object will use</param>
        /// <returns></returns>
        public static Roulette CreateSafely(double[] pocketSizes, IGenerator randomNumberGenerator)
        {
            // Check if pocket sizes are positive
            for (int i = 0; i < pocketSizes.Length; ++i)
            {
                if (pocketSizes[i] <= 0)
                {
                    throw new ArgumentException("Pocket siezes must be positive, but at index = " + i + 
                                                " pocket size = " + pocketSizes[i] );
                }
            }

            // Normalize them so that their sum = 1.0
            double[] normalizedPocketSizes = new double[pocketSizes.Length];
            Normalize(pocketSizes);
            return new Roulette(normalizedPocketSizes, randomNumberGenerator);
        }

        /// <summary>
        /// Simulates spinning the roulette wheel while rolling a ball.
        /// Returns a random index from 0 to pocketsCount-1. 
        /// The probability of each index equals the size of its corresponding pocket.
        /// </summary>
        /// <returns>A random index from 0 to pocketsCount-1</returns>
        public int SpinWheelWithBall()
        {
            double ball = rng.NextDouble();
            return LookUp(ball);
        }

        /// <summary>
        /// Simulates spinning the roulette wheel with many pointers, placed at equal distances, and then shuffles the results. 
        /// The probability that a pointer lands in a pocket is equal to the pocket's size.
        /// Useful for universal stochastic sampling.
        /// </summary>
        /// <param name="pointersCount">The number of pointers. 
        ///                             It must be >=2, otherwise an <see cref="ArgumentException"/> will be thrown</param>
        /// <returns>An array of indexes from 0 to pocketsCount-1</returns>
        public int[] SpinWheelWithPointers(int pointersCount)
        {
            if (pointersCount < 2) throw new ArgumentException("There must be at least 2 pointers, but there were "
                                                                + pointersCount);
            double pointerDistance = 1.0 / pointersCount;
            double firstPointer = rng.NextDouble(pointerDistance);

            int[] indexes = new int[pointersCount];
            for (int i = 0; i < pointersCount; ++i)
            {
                double pointer = firstPointer + i * pointerDistance;
                indexes[i] = LookUp(pointer);
            }
            RandomNumberGenerationUtilities.Shuffle<int>(indexes);
            return indexes;
        }

        // TODO: employ (iterative) binary search to get to O(log(pocketsSize))
        private int LookUp(double value)
        {
            int pocket;
            for (pocket = 0; pocket < wheel.Length; ++pocket)
            {
                if (value < wheel[pocket]) return pocket;
            }
            return pocket; // pocket == wheel.Length at this point
        }

        private static void Normalize(double[] pocketSizes)
        {
            double sum = 0.0;
            for (int i = 0; i < pocketSizes.Length; ++i) sum += pocketSizes[i];
            for (int i = 0; i < pocketSizes.Length; ++i) pocketSizes[i] /= sum;
        }
    }
}
