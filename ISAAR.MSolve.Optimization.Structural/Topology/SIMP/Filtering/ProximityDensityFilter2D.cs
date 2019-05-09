using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: concolution operations can be parallelized very efficiently
namespace ISAAR.MSolve.Optimization.Structural.Topology.SIMP.Filtering
{
    /// <summary>
    /// Filters the density of each finite element by taking into account the density and compliance derivative of other 
    /// elements near it. See "An efficient 3D topology optimization code written in Matlab - K. Liu, A. Tovar, 2014".
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class ProximityDensityFilter2D : IDensityFilter
    {
        private const double minDensityInDenominator = 1E-3; //TODO: Should the user be able to change this? 

        private readonly double filterAreaRadius;
        private readonly ContinuumElement2D[] elements;
        private readonly ElementNeighborhood[] elementNeighborhoods;
        private readonly int numElements;

        /// <summary>
        /// </summary>
        /// <param name="model"></param>
        /// <param name="filterAreaRadius">
        /// Typically 10% of the smaller dimension of the design domain. If 0, the unfiltered sensitivities are retained. 
        /// </param>
        public ProximityDensityFilter2D(Model model, double filterAreaRadius)
        {
            this.filterAreaRadius = filterAreaRadius;

            // Cast the continuum elements. //TODO: this process is already performed in LinearFemAnalysis2DGeneral.
            IList<Element> modelElements = model.Elements;
            numElements = modelElements.Count;
            elements = new ContinuumElement2D[numElements];
            for (int e = 0; e < numElements; ++e)
            {
                if (modelElements[e].ElementType is ContinuumElement2D continuum) elements[e] = continuum;
                else throw new ArgumentException("2D topology optimization only works with 2D continuum elements,"
                    + $" but the element with ID = {modelElements[e].ID} was not.");
            }

            // Find the neighbors of each element and precalculate the quantities that are constant
            elementNeighborhoods = new ElementNeighborhood[numElements];
            for (int e = 0; e < numElements; ++e) elementNeighborhoods[e] = FindNeighbors(filterAreaRadius, elements, e);
        }

        public void FilterSensitivities(Vector densities, ref Vector sensitivities)
        {
            var filteredSensitivities = Vector.CreateZero(numElements);
            for (int e = 0; e < numElements; ++e)
            {
                double convolution = 0.0;
                int[] neighbors = elementNeighborhoods[e].NeighborElementIndices;
                double[] weightFactors = elementNeighborhoods[e].WeightFactors;
                for (int f = 0; f < neighbors.Length; ++f)
                {
                    convolution += weightFactors[f] * densities[neighbors[f]] * sensitivities[neighbors[f]];
                }
                double denominator = Math.Max(minDensityInDenominator, densities[e]) * elementNeighborhoods[e].WeightFactorSum;
                filteredSensitivities[e] = convolution / denominator;
            }
            sensitivities = filteredSensitivities;
        }

        private static ElementNeighborhood FindNeighbors(double radius, ContinuumElement2D[] allElements, 
            int elementIdx)
        {
            //TODO: Using appropriate mesh classes, we do not have to process all elements. Instead we can focus on neighbouring
            //      once until the radius is exceeded.
            var neighborIndices = new List<int>();
            var weightFactors = new List<double>();
            CartesianPoint thisCentroid = allElements[elementIdx].FindCentroid();
            for (int f = 0; f < allElements.Length; ++f)
            {
                CartesianPoint otherCentroid = allElements[f].FindCentroid();
                double distance = otherCentroid.CalculateDistanceFrom(thisCentroid);
                if (distance <= radius)
                {
                    neighborIndices.Add(f);
                    weightFactors.Add(radius - distance);
                }
            }
            return new ElementNeighborhood(neighborIndices.ToArray(), weightFactors.ToArray());
        }

        private class ElementNeighborhood
        {
            internal ElementNeighborhood(int[] neighborElementIndices, double[] weightFactors)
            {
                this.NeighborElementIndices = neighborElementIndices;
                this.WeightFactors = weightFactors;
                this.WeightFactorSum = weightFactors.Sum();
            }

            internal int[] NeighborElementIndices { get; }
            internal double[] WeightFactors { get; }
            internal double WeightFactorSum { get; }
        }
    }
}
