using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: When the general classes achieve the efficiency of the 2D uniform hardcoded classes (such as this), then the latter
//      must be removed.
namespace ISAAR.MSolve.Optimization.Structural.Topology.SIMP.Filtering
{
    public class MeshIndependentSensitivityFilter2DUniform: IDensityFilter
    {
        private const double minDensityInDenominator = 1E-3; //TODO: Should the user be able to change this? 

        private readonly int numElementsX, numElementsY;
        private readonly double filterAreaRadius;

        public MeshIndependentSensitivityFilter2DUniform(int numElementsX, int numElementsY, double filterAreaRadius)
        {
            this.numElementsX = numElementsX;
            this.numElementsY = numElementsY;
            this.filterAreaRadius = filterAreaRadius;
        }

        public void FilterSensitivities(Vector densities, ref Vector sensitivities)
        {
            int dimension = densities.Length;
            var filtered = Vector.CreateZero(dimension);
            int radiusRounded = (int)Math.Round(filterAreaRadius);
            for (int i = 1; i <= numElementsX; ++i)
            {
                int kStart = Math.Max(i - radiusRounded, 1);
                int kEnd = Math.Min(i + radiusRounded, numElementsX);
                for (int j = 1; j <= numElementsY; ++j)
                {
                    int mStart = Math.Max(j - radiusRounded, 1);
                    int mEnd = Math.Min(j + radiusRounded, numElementsY);
                    int elementIdx = Find1DIndex(i - 1, j - 1);
                    double sum = 0.0;
                    for (int k = kStart; k <= kEnd; ++k)
                    {
                        for (int m = mStart; m <= mEnd; ++m)
                        {
                            int otherElementIdx = Find1DIndex(k - 1, m - 1);
                            double fac = filterAreaRadius - Math.Sqrt((i - k) * (i - k) + (j - m) * (j - m));
                            if (fac > 0)
                            {
                                sum += fac;
                                filtered[elementIdx] += fac * densities[otherElementIdx] * sensitivities[otherElementIdx];
                            }
                        }
                    }
                    filtered[elementIdx] /= Math.Max(minDensityInDenominator, densities[elementIdx]) * sum;
                }
            }
            sensitivities = filtered;
        }

        /// <summary>
        /// Returns the (0-based) index into a 1D array/vector of the element with the provided 2D indices (also 0-based).
        /// </summary>
        /// <param name="xIdx">0-based</param>
        /// <param name="yIdx">0-based</param>
        private int Find1DIndex(int xIdx, int yIdx) => xIdx + yIdx * numElementsX;
    }
}
