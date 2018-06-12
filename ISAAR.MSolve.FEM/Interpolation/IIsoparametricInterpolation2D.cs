using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation.Inverse;
using ISAAR.MSolve.Geometry.Coordinates;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Interpolation
{
    /// <summary>
    /// Interpolation (or shape or basis) functions used by isoparametric finite elements. Instances are able to compute the 
    /// values of the shape functions, their derivatives, the Jacobian of the matrix, etc.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IIsoparametricInterpolation2D
    {
        /// <summary>
        /// The coordinates of the finite element's nodes in the natural (element local) coordinate system. The order of these
        /// nodes matches the order of the shape functions and is always the same for each element.
        /// </summary>
        IReadOnlyList<NaturalPoint2D> NodalNaturalCoordinates { get; }

        /// <summary>
        /// The number of shape functions that define this interpolation.
        /// </summary>
        int NumFunctions { get; }

        /// <summary>
        /// The inverse mapping of this interpolation, namely from global cartesian to natural (element local) coordinate system.
        /// </summary>
        /// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
        /// <returns></returns>
        IInverseInterpolation2D CreateInverseMappingFor(IReadOnlyList<Node2D> nodes);

        /// <summary>
        /// Evaluate the shape functions and 1st order derivatives with respect to the global cartesian coordinate system 
        /// at a given natural point.
        /// </summary>
        /// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
        /// <param name="naturalPoint">The coordinates of the point in the natural (element local) coordinate system of the 
        ///     element.</param>
        /// <returns></returns>
        EvalInterpolation2D EvaluateAllAt(IReadOnlyList<Node2D> nodes, NaturalPoint2D naturalPoint);

        /// <summary>
        /// Evaluate the shape functions at a given natural point.
        /// </summary>
        /// <param name="naturalPoint">The coordinates of the point in the natural (element local) coordinate system of the 
        ///     element.</param>
        /// <returns></returns>
        EvalShapeFunctions2D EvaluateFunctionsAt(NaturalPoint2D naturalPoint);

        /// <summary>
        /// Evaluate the 1st order shape function derivatives with respect to the global cartesian coordinate system 
        /// at a given natural point.
        /// </summary>
        /// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
        /// <param name="naturalPoint">The coordinates of the point in the natural (element local) coordinate system of the 
        ///     element.</param>
        /// <returns></returns>
        EvalShapeGradients2D EvaluateGradientsAt(IReadOnlyList<Node2D> nodes, NaturalPoint2D naturalPoint);

        CartesianPoint2D TransformNaturalToCartesian(IReadOnlyList<Node2D> nodes, NaturalPoint2D naturalPoint);
    }
}
