using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Numerical.LinearAlgebra;

//TODO: add tensor product for LinearAlgebra.Vectors.IVectorView. It is useful in mass matrix calculations.
namespace ISAAR.MSolve.FEM.Interpolation
{
    /// <summary>
    /// Stores the shape functions of an interpolation, evaluated at a certain natural point of a finite element.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class EvalShapeFunctions2D
    {
        private readonly double[] shapeFunctions;

        public EvalShapeFunctions2D(double[] shapeFunctions)
        {
            this.shapeFunctions = shapeFunctions;
        }

        /// <summary>
        /// The value of the stored shape function that corresponds to the node with local index <paramref name="nodeIdx"/>.
        /// </summary>
        /// <param name="nodeIdx">The local index of the node, namely its order among the nodes of the finite element.</param>
        /// <returns></returns>
        public double this[int nodeIdx] => shapeFunctions[nodeIdx];

        /// <summary>
        /// The shape function matrix is 2-by-2n, where n = is the number of shape functions. Row 0 corresponds to dof X, while
        /// row 1 to dof Y.
        /// </summary>
        /// <returns></returns>
        public Matrix2D BuildShapeFunctionMatrix()
        {
            var array2D = new double[2, 2 * shapeFunctions.Length];
            for (int i = 0; i < shapeFunctions.Length; ++i)
            {
                array2D[0, 2 * i] = shapeFunctions[i];
                array2D[1, 2 * i + 1] = shapeFunctions[i];
            }
            return new Matrix2D(array2D);
        }

        public CartesianPoint2D TransformPointNaturalToGlobalCartesian(IReadOnlyList<Node2D> nodes)
        {
            if (nodes.Count != shapeFunctions.Length) throw new ArgumentException(
                $"There are {shapeFunctions.Length} evaluated shape functions stored, but {nodes.Count} were passed in.");
            double x = 0, y = 0;
            for (int i = 0; i < shapeFunctions.Length; ++i)
            {
                Node2D node = nodes[i];
                double val = shapeFunctions[i];
                x += val * node.X;
                y += val * node.Y;
            }
            return new CartesianPoint2D(x, y);
        }
    }
}
