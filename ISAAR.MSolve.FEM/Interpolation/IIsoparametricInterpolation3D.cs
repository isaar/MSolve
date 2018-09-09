using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Integration.Points;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation.Inverse;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.FEM.Interpolation
{
	/// <summary>
	/// Interpolation (or shape or basis) functions used by isoparametric finite elements. Instances are able to compute the 
	/// values of the shape functions, their derivatives, the Jacobian of the matrix, etc.
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public interface IIsoparametricInterpolation3D
    {
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
		/// <returns></returns>
		IInverseInterpolation3D CreateInverseMappingFor(IReadOnlyList<Node3D> nodes);

		/// <summary>
		/// Evaluate the shape functions and 1st order derivatives with respect to the global cartesian coordinate system 
		/// at a given natural point.
		/// </summary>
		/// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
		/// <param name="naturalPoint">The coordinates of the point in the natural (element local) coordinate system of the 
		///     element.</param>
		/// <returns></returns>
		EvalInterpolation3D EvaluateAllAt(IReadOnlyList<Node3D> nodes, NaturalPoint3D naturalPoint);

		/// <summary>
		/// Evaluate the shape functions and 1st order derivatives with respect to the global cartesian coordinate system 
		/// at the integration points defined by a given quadrature. This method caches all possible quantities from previous
		/// calls. Use it instead of <see cref="EvaluateAllAt(IReadOnlyList{Node3D}, NaturalPoint3D)"/>
		/// when the integration points of an element are the same across multiple elements or multiple iterations of a
		/// non linear or dynamic analysis (in the latter cases we could also cache the <see cref="EvalInterpolation3D"/> at 
		/// the integration points of each element).
		/// </summary>
		/// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
		/// <param name="quadrature">The integration rule that defines integration points where shape functions and derivatives  
		///     are calculated. The integration points of this instance of <see cref="IQuadrature3D"/> are always the 
		///     same.</param>
		/// <returns></returns>
		Dictionary<GaussPoint3D, EvalInterpolation3D> EvaluateAllAtGaussPoints(IReadOnlyList<Node3D> nodes,
		    IQuadrature3D quadrature);

		/// <summary>
		/// Evaluate the shape functions at a given natural point.
		/// </summary>
		/// <param name="naturalPoint">The coordinates of the point in the natural (element local) coordinate system of the 
		///     element.</param>
		/// <returns></returns>
		EvalShapeFunctions3D EvaluateFunctionsAt(NaturalPoint3D naturalPoint);

		/// <summary>
		/// Evaluate the shape functions at the integration points defined by a given quadrature. This method caches all possible
		/// quantities from previous calls. Use it instead of <see cref="EvaluateAllAt(IReadOnlyList{Node3D}, NaturalPoint3D)"/>
		/// when the integration points of an element are the same across multiple elements or multiple iterations of a
		/// non linear or dynamic analysis.
		/// </summary>
		/// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
		/// <param name="quadrature">The integration rule that defines integration points where shape functions are calculated. 
		///     The integration points of this instance of <see cref="IQuadrature2D"/> are always the same.</param>
		/// <returns></returns>
		Dictionary<GaussPoint3D, EvalShapeFunctions3D> EvaluateFunctionsAtGaussPoints(IReadOnlyList<Node3D> nodes,
		    IQuadrature3D quadrature);

		/// <summary>
		/// Evaluate the 1st order shape function derivatives with respect to the global cartesian coordinate system 
		/// at a given natural point.
		/// </summary>
		/// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
		/// <param name="naturalPoint">The coordinates of the point in the natural (element local) coordinate system of the 
		///     element.</param>
		/// <returns></returns>
		EvalShapeGradients3D EvaluateGradientsAt(IReadOnlyList<Node3D> nodes, NaturalPoint3D naturalPoint);

		/// <summary>
		/// Evaluate the 1st order shape function derivatives with respect to the global cartesian coordinate system 
		/// at the integration points defined by a given quadrature. This method caches all possible quantities from previous
		/// calls. Use it instead of <see cref="EvaluateGradientsAt(IReadOnlyList{Node3D}, NaturalPoint3D)"/> 
		/// when the integration points of an element are the same across multiple elements or multiple iterations of a
		/// non linear or dynamic analysis (in the latter cases we could also cache the <see cref="EvalShapeGradients3D"/> at 
		/// the integration points of each element).
		/// </summary>
		/// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
		/// <param name="quadrature">The integration rule that defines integration points where shape function derivatives are 
		///     calculated. The integration points of this instance of <see cref="IQuadrature3D"/> are always the same.</param>
		/// <returns></returns>
		Dictionary<GaussPoint3D, EvalShapeGradients3D> EvaluateGradientsAtGaussPoints(IReadOnlyList<Node3D> nodes,
		    IQuadrature3D quadrature);

		/// <summary>
		/// Transforms the coordinates from the natural (element local) coordinate system to the the global
		/// coordinate system of a point that is internal to the finite element.
		/// </summary>
		/// <param name="nodes">The coordinates of the finite element's nodes in the global cartesian system.</param>
		/// <param name="naturalPoint">The coordinates in the natural system of a point that is internal to the finite 
		///     element.</param>
		/// <returns></returns>
		CartesianPoint3D TransformNaturalToCartesian(IReadOnlyList<Node3D> nodes, NaturalPoint3D naturalPoint);
    }
}
