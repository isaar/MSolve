using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Numerical.LinearAlgebra;

//TODO: In XFEM I used dictionaries with nodes as keys, but it is less efficient and offers little extra safety, since the
//      shape functions and derivatives will be used by classes that have direct access to the nodes.
namespace ISAAR.MSolve.FEM.Interpolation
{
    /// <summary>
    /// Stores the shape functions, 1st order derivatives with respect to the global cartesian coordinates and the Jacobian
    /// of an interpolation, evaluated at a certain natural point of a finite element. These quantities are needed in many 
    /// places, thus passing an instance of this class is less verbose and error prone.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class EvalInterpolation2D
    {
        private readonly double[] shapeFunctions;
        private readonly double[][] shapeGradientsCartesian;

        public EvalInterpolation2D(double[] shapeFunctions, double[,] shapeGradientsNatural, Jacobian2D jacobian)
        {
            int numNodes = shapeFunctions.Length;
            if (shapeGradientsNatural.GetLength(0) != numNodes) throw new ArgumentException($"There are {shapeFunctions.Length}"
               + $" evaluated shape functions, but {shapeGradientsNatural.GetLength(0)} evaluated natural shape derivatives.");
            this.shapeFunctions = shapeFunctions;
            this.shapeGradientsCartesian = new double[numNodes][];
            for (int i = 0; i < numNodes; ++i)
            {
                this.shapeGradientsCartesian[i] = jacobian.TransformNaturalDerivativesToCartesian(
                    shapeGradientsNatural[i, 0], shapeGradientsNatural[i, 1]);
            }
            this.Jacobian = jacobian;
        }

        /// <summary>
        /// The inverse Jacobian matrix of the interpolation and its determinant.
        /// </summary>
        public Jacobian2D Jacobian { get; }

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

        /// <summary>
        /// The value of the stored shape function that corresponds to the node with local index <paramref name="nodeIdx"/>.
        /// </summary>
        /// <param name="nodeIdx">The local index of the node, namely its order among the nodes of the finite element.</param>
        /// <returns></returns>
        public double GetShapeFunction(int nodeIdx) => shapeFunctions[nodeIdx];

        /// <summary>
        /// The values of the stored shape function derivatives, with respect to the global cartesian coordinates, that 
        /// correspond to the node with local index <paramref name="nodeIdx"/>.
        /// </summary>
        /// <param name="nodeIdx">The local index of the node, namely its order among the nodes of the finite element.</param>
        /// <returns></returns>
        public IReadOnlyList<double> GetShapeGradientCartesian(int nodeIdx) => shapeGradientsCartesian[nodeIdx];

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
