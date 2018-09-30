using System;
using System.Collections.Generic;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation.Jacobians;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Numerical.LinearAlgebra;

namespace ISAAR.MSolve.FEM.Interpolation
{
    /// <summary>
    /// Stores the shape functions, 1st order derivatives with respect to the global cartesian coordinates and the Jacobian
    /// of an interpolation, evaluated at a certain natural point of a finite element. These quantities are needed in many 
    /// places, thus passing an instance of this class is less verbose and error prone.
    /// Authors: Dimitris Tsapetis
    /// </summary>
    public class EvalInterpolation3D
    {
        public EvalInterpolation3D(Vector shapeFunctions, Matrix2D shapeGradientsNatural, IsoparametricJacobian3D jacobian)
        {
            int numberOfNodes = shapeFunctions.Length;
            if (shapeGradientsNatural.Rows != numberOfNodes) throw new ArgumentException($"There are {shapeFunctions.Length}"
                + $" evaluated shape functions, but {shapeGradientsNatural.Rows} evaluated natural shape derivatives.");
            this.ShapeFunctions = shapeFunctions;
            this.ShapeGradientsNatural = shapeGradientsNatural;
            this.Jacobian = jacobian;
            this.ShapeGradientsCartesian = jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural);
        }

        /// <summary>
        /// The inverse Jacobian matrix of the interpolation and its determinant.
        /// </summary>
        public IsoparametricJacobian3D Jacobian { get; }

        /// <summary>
        /// A vector that contains the shape functions in the same order as the nodes of the interpolation.
        /// </summary>
        public Vector ShapeFunctions { get; }

        /// <summary>
        /// A matrix that contains the 1st order shape function derivatives with respect to the global cartesian coordinate 
        /// system at the integration points defined by a given quadrature. Each row corresponds to the gradient of a single 
        /// shape function. Each column corresponds to the derivatives of all shape functions with respect to a single 
        /// coordinate.
        /// </summary>
        public Matrix2D ShapeGradientsCartesian { get; }

        /// <summary>
        /// A matrix that contains the 1st order shape function derivatives with respect to the natural coordinate 
        /// system at the integration points defined by a given quadrature. Each row corresponds to the gradient of a single 
        /// shape function. Each column corresponds to the derivatives of all shape functions with respect to a single 
        /// coordinate.
        /// </summary>
        public Matrix2D ShapeGradientsNatural { get; }

        public CartesianPoint3D TransformPointNaturalToGlobalCartesian(IReadOnlyList<Node3D> nodes, NaturalPoint3D naturalCoordinates)
        {
            if (nodes.Count != ShapeFunctions.Length) throw new ArgumentException(
                $"There are {ShapeFunctions.Length} evaluated shape functions stored, but {nodes.Count} were passed in.");
            double x = 0, y = 0, z = 0;
            for (int i = 0; i < ShapeFunctions.Length; i++)
            {
                Node3D node = nodes[i];
                x += ShapeFunctions[i] * node.X;
                y += ShapeFunctions[i] * node.Y;
                z += ShapeFunctions[i] * node.Z;
            }
            return new CartesianPoint3D(x,y,z);
        }
    }
}
