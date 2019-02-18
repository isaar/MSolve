using System;
using System.Collections.Generic;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Numerical.LinearAlgebra;

namespace ISAAR.MSolve.FEM.Interpolation
{
    /// <summary>
    /// Authors: Serafeim Bakalakos, Dimitris Tsapetis
    /// </summary>
    public static class InterpolationUtilities
    {
        public static CartesianPoint2D TransformPointNaturalToGlobalCartesian(IReadOnlyList<Node_v2> nodes, 
            Vector shapeFunctionsAtNaturalPoint)
        {
            int numFuncs = shapeFunctionsAtNaturalPoint.Length;
            if (nodes.Count != numFuncs) throw new ArgumentException(
                $"There are {numFuncs} evaluated shape functions stored, but {nodes.Count} were passed in.");
            double x = 0, y = 0;
            for (int i = 0; i < numFuncs; ++i)
            {
                Node_v2 node = nodes[i];
                double val = shapeFunctionsAtNaturalPoint[i];
                x += val * node.X;
                y += val * node.Y;
            }
            return new CartesianPoint2D(x, y);
        }

        public static CartesianPoint3D TransformPointToGlobalCartesian(IReadOnlyList<Node_v2> nodes,
            Vector shapeFunctionsAtNaturalPoint)
        {
            int numFuncs = shapeFunctionsAtNaturalPoint.Length;
            if (nodes.Count != numFuncs) throw new ArgumentException(
                $"There are {numFuncs} evaluated shape functions stored, but {nodes.Count} were passed in.");
            double x = 0, y = 0, z = 0;
            for (int i = 0; i < numFuncs; ++i)
            {
                Node_v2 node = nodes[i];
                double val = shapeFunctionsAtNaturalPoint[i];
                x += val * node.X;
                y += val * node.Y;
                z += val * node.Z;
            }

            return new CartesianPoint3D(x, y, z);
        }
    }
}
