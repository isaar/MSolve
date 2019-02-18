using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Integration.Points;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation.Inverse;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.Numerical.LinearAlgebra;

namespace ISAAR.MSolve.FEM.Interpolation
{
    /// <summary>
    /// Interpolation (or shape or basis) functions used by isoparametric finite elements. Instances are able to compute the 
    /// values of the shape functions, their derivatives, the Jacobian of the matrix, etc.
    /// Authors: Dimitris Tsapetis, Serafeim Bakalakos
    /// </summary>
    public interface IIsoparametricInterpolation3D_OLD
    {
        /// <summary>
        /// The shape of a cell. Useful for interacting with other modules and software.
        /// </summary>
        CellType CellType { get; }

        /// <summary>
        /// The coordinates of the finite element's nodes in the natural (element local) coordinate system. The order of these
        /// nodes matches the order of the shape functions and is always the same for each element.
        /// </summary>
        IReadOnlyList<NaturalPoint3D> NodalNaturalCoordinates { get; }

        /// <summary>
        /// The number of shape functions that define this interpolation.
        /// </summary>
        int NumFunctions { get; }

        /// <summary>
        /// The inverse mapping of this interpolation, namely from global cartesian to natural (element local) coordinate system.
        /// </summary>
        /// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
        IInverseInterpolation3D CreateInverseMappingFor(IReadOnlyList<Node3D> nodes);

        /// <summary>
        /// Evaluate the shape functions, the 1st order derivatives with respect to the global cartesian coordinate system 
        /// and the jacobian at a given natural point.
        /// </summary>
        /// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
        /// <param name="naturalPoint">The coordinates of the point in the natural (element local) coordinate system of the 
        ///     element.</param>
        EvalInterpolation3D_OLD EvaluateAllAt(IReadOnlyList<Node3D> nodes, NaturalPoint3D naturalPoint);

        /// <summary>
        /// Evaluate the shape functions, the 1st order derivatives with respect to the global cartesian coordinate system and
        /// the jacobian at the integration points defined by a given quadrature. 
        /// This method caches all possible quantities from previous calls. Use it instead of when the integration points of an 
        /// <see cref="EvaluateAllAt(IReadOnlyList{Node3D}, NaturalPoint3D)"/>
        /// element are the same across multiple elements or multiple iterations of a non linear or dynamic analysis (in the 
        /// latter cases we could also cache the <see cref="EvalInterpolation3D"/> at the integration points of each element).
        /// The <see cref="EvalInterpolation3D"/>s per integration point are returned in the same order as the integration 
        /// points in <paramref name="quadrature"/>.<see cref="IQuadrature3D.IntegrationPoints"/>.
        /// </summary>
        /// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
        /// <param name="quadrature">The integration rule that defines integration points where shape functions and derivatives  
        ///     are calculated. The integration points of this instance of <see cref="IQuadrature3D"/> are always the 
        ///     same.</param>
        IReadOnlyList<EvalInterpolation3D_OLD> EvaluateAllAtGaussPoints(IReadOnlyList<Node3D> nodes, IQuadrature3D quadrature);

        /// <summary>
        /// Evaluate the shape functions at a given natural point and returns them in a vector in the same order as the nodes
        /// of the interpolation.
        /// </summary>
        /// <param name="naturalPoint">The coordinates of the point in the natural (element local) coordinate system of the 
        ///     element.</param>
        Vector EvaluateFunctionsAt(NaturalPoint3D naturalPoint);

        /// <summary>
        /// Evaluate the shape functions at the integration points defined by a given quadrature. 
        /// This method caches all possible quantities from previous calls. Use it instead of 
        /// <see cref="EvaluateFunctionsAt(NaturalPoint3D)"/> when the integration points of an element are the same across 
        /// multiple elements or multiple iterations of a non linear or dynamic analysis.
        /// The shape functions vectors per integration point are returned in the same order as the integration 
        /// points in <paramref name="quadrature"/>.<see cref="IQuadrature2D.IntegrationPoints"/>. Each vector contains the
        /// shape functions in the same order as the nodes of the interpolation.
        /// </summary>
        /// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
        /// <param name="quadrature">The integration rule that defines integration points where shape functions are calculated. 
        ///     The integration points of this instance of <see cref="IQuadrature2D"/> are always the same.</param>
        IReadOnlyList<Vector> EvaluateFunctionsAtGaussPoints(IQuadrature3D quadrature);

        /// <summary>
        /// Evaluate the 1st order shape function derivatives with respect to the natural coordinate system 
        /// at a given natural point.
        /// Each row of the returned matrix corresponds to the gradient of a single shape function. Each column corresponds to 
        /// the derivatives of all shape functions with respect to a single coordinate.
        /// </summary>
        /// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
        /// <param name="naturalPoint">The coordinates of the point in the natural (element local) coordinate system of the 
        ///     element.</param>
        Matrix2D EvaluateNaturalGradientsAt(NaturalPoint3D naturalPoint);

        /// <summary>
        /// Evaluate the 1st order shape function derivatives with respect to the natural coordinate system 
        /// at the integration points defined by a given quadrature. 
        /// This method caches all possible quantities from previous calls. Use it instead of  
        /// <see cref="EvaluateNaturalGradientsAt(NaturalPoint3D)"/> when the integration points of an element are the same
        /// across multiple elements or multiple iterations of a non linear or dynamic analysis (in the latter cases we could
        /// also cache the cartesian gradients at the integration points of each element).
        /// The shape gradients matrices per integration point are returned in the same order as the integration 
        /// points in <paramref name="quadrature"/>.<see cref="IQuadrature2D.IntegrationPoints"/>. Each row of a matrix  
        /// corresponds to the gradient of a single shape function. Each column corresponds to the derivatives of all shape
        /// functions with respect to a single coordinate.
        /// </summary>
        /// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
        /// <param name="quadrature">The integration rule that defines integration points where shape function derivatives are 
        ///     calculated. The integration points of this instance of <see cref="IQuadrature3D"/> are always the same.</param>
        IReadOnlyList<Matrix2D> EvaluateNaturalGradientsAtGaussPoints(IQuadrature3D quadrature);

        /// <summary>
        /// Transforms the coordinates from the natural (element local) coordinate system to the the global
        /// coordinate system of a point that is internal to the finite element.
        /// </summary>
        /// <param name="nodes">The coordinates of the finite element's nodes in the global cartesian system.</param>
        /// <param name="naturalPoint">The coordinates in the natural system of a point that is internal to the finite 
        ///     element.</param>
        CartesianPoint3D TransformNaturalToCartesian(IReadOnlyList<Node3D> nodes, NaturalPoint3D naturalPoint);
    }
}
