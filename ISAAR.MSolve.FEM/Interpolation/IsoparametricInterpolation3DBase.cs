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
	/// Basic implementation of <see cref="IIsoparametricInterpolation3D"/>
	/// Authors: Dimitris Tsapetis
	/// </summary>
    public abstract class IsoparametricInterpolation3DBase : IIsoparametricInterpolation3D
	{
		private readonly Dictionary<IQuadrature3D, Dictionary<GaussPoint3D, double[,]>> cachedNaturalShapederivativesAtGPs;
		private readonly Dictionary<IQuadrature3D, Dictionary<GaussPoint3D, double[]>> cachedRawShapeFunctionsAtGPs;
		private readonly Dictionary<IQuadrature3D, Dictionary<GaussPoint3D, EvalShapeFunctions3D>> cachedShapeFunctionsAtGPs;

		public IsoparametricInterpolation3DBase(int numFunctions)
		{
			this.NumFunctions = numFunctions;
			this.cachedNaturalShapederivativesAtGPs = new Dictionary<IQuadrature3D, Dictionary<GaussPoint3D, double[,]>>();
			this.cachedRawShapeFunctionsAtGPs = new Dictionary<IQuadrature3D, Dictionary<GaussPoint3D, double[]>>();
			this.cachedShapeFunctionsAtGPs = new Dictionary<IQuadrature3D, Dictionary<GaussPoint3D, EvalShapeFunctions3D>>();
		}

		/// <summary>
		/// The number of shape functions that define this interpolation.
		/// </summary>
		public int NumFunctions { get; }

		/// <summary>
		/// The coordinates of the finite element's nodes in the natural (element local) coordinate system. The order of these
		/// nodes matches the order of the shape functions and is always the same for each element.
		/// </summary>
		public abstract IReadOnlyList<NaturalPoint3D> NodalNaturalCoordinates { get; }

		/// <summary>
		/// The inverse mapping of this interpolation, namely from global cartesian to natural (element local) coordinate system.
		/// </summary>
		/// <param name="node">The nodes of the finite element in the global cartesian coordinate system.</param>
		/// <returns></returns>
		public abstract IInverseInterpolation3D CreateInverseMappingFor(IReadOnlyList<Node3D> node);

		/// <summary>
		/// Evaluate the shape functions and 1st order derivatives with respect to the global cartesian coordinate system 
		/// at a given natural point.
		/// </summary>
		/// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
		/// <param name="naturalPoint">The coordinates of the point in the natural (element local) coordinate system of the 
		///     element.</param>
		/// <returns></returns>
		public EvalInterpolation3D EvaluateAllAt(IReadOnlyList<Node3D> nodes, NaturalPoint3D naturalPoint)
		{
			double xi = naturalPoint.Xi;
			double eta = naturalPoint.Eta;
			double zeta = naturalPoint.Zeta;
			double[,] naturalShapeDerivatives = EvaluateGradientsAt(xi, eta, zeta);
			return new EvalInterpolation3D(EvaluateAt(xi, eta, zeta), naturalShapeDerivatives,
				new Jacobian3D(nodes, naturalShapeDerivatives));
		}

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
		///     are calculated. The integration points of this instance of <see cref="IQuadrature2D"/> are always the 
		///     same.</param>
		/// <returns></returns>
		public Dictionary<GaussPoint3D, EvalInterpolation3D> EvaluateAllAtGaussPoints(IReadOnlyList<Node3D> nodes,
			IQuadrature3D quadrature)
		{
			Dictionary<GaussPoint3D, double[]> shapeFunctionsAtGPs = EvaluateShapeFunctionsAtGPs(quadrature);
			Dictionary<GaussPoint3D, double[,]> naturalShapeDerivativesAtGPs = EvaluateNaturalShapeDerivativesAtGPs(quadrature);

			var interpolationAtGPs = new Dictionary<GaussPoint3D, EvalInterpolation3D>();
			foreach (var gaussPoint in quadrature.IntegrationPoints)
			{
				interpolationAtGPs[gaussPoint] = new EvalInterpolation3D(shapeFunctionsAtGPs[gaussPoint],
					naturalShapeDerivativesAtGPs[gaussPoint], new Jacobian3D(nodes, naturalShapeDerivativesAtGPs[gaussPoint]));
			}

			return interpolationAtGPs;
		}

		/// <summary>
		/// Evaluate the shape functions at a given natural point.
		/// </summary>
		/// <param name="naturalPoint">The coordinates of the point in the natural (element local) coordinate system of the 
		///     element.</param>
		/// <returns></returns>
		public EvalShapeFunctions3D EvaluateFunctionsAt(NaturalPoint3D naturalPoint)
		{
			return new EvalShapeFunctions3D(EvaluateAt(naturalPoint.Xi, naturalPoint.Eta, naturalPoint.Zeta));
		}

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
		public Dictionary<GaussPoint3D, EvalShapeFunctions3D> EvaluateFunctionsAtGaussPoints(IReadOnlyList<Node3D> nodes,
			IQuadrature3D quadrature)
		{
			bool isCached = cachedShapeFunctionsAtGPs.TryGetValue(quadrature,
				out Dictionary<GaussPoint3D, EvalShapeFunctions3D> shapeFunctionsAtGPs);
			if (!isCached)
			{
				shapeFunctionsAtGPs= new Dictionary<GaussPoint3D, EvalShapeFunctions3D>();
				foreach (var gaussPoint in quadrature.IntegrationPoints)
				{
					shapeFunctionsAtGPs[gaussPoint]=new EvalShapeFunctions3D(EvaluateAt(gaussPoint.Xi, gaussPoint.Eta, gaussPoint.Zeta));
				}
				cachedShapeFunctionsAtGPs.Add(quadrature,shapeFunctionsAtGPs);
			}

			return shapeFunctionsAtGPs;
		}

		/// <summary>
		/// Evaluate the 1st order shape function derivatives with respect to the global cartesian coordinate system 
		/// at a given natural point.
		/// </summary>
		/// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
		/// <param name="naturalPoint">The coordinates of the point in the natural (element local) coordinate system of the 
		///     element.</param>
		/// <returns></returns>
		public EvalShapeGradients3D EvaluateGradientsAt(IReadOnlyList<Node3D> nodes, NaturalPoint3D naturalPoint)
		{
			double[,] naturalShapeDerivatives = EvaluateGradientsAt(naturalPoint.Xi, naturalPoint.Eta, naturalPoint.Zeta);
			return new EvalShapeGradients3D(naturalShapeDerivatives, new Jacobian3D(nodes, naturalShapeDerivatives));
		}

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
		public Dictionary<GaussPoint3D, EvalShapeGradients3D> EvaluateGradientsAtGaussPoints(IReadOnlyList<Node3D> nodes,
			IQuadrature3D quadrature)
		{
			Dictionary<GaussPoint3D, double[,]> naturalShapeDerivativesAtGPs = EvaluateNaturalShapeDerivativesAtGPs(quadrature);

			var shapeGradientsAtGPs = new Dictionary<GaussPoint3D, EvalShapeGradients3D>();
			foreach (var gpNaturalDerivativesPair in naturalShapeDerivativesAtGPs)
			{
				GaussPoint3D gaussPoint = gpNaturalDerivativesPair.Key;
				double[,] naturalShapeDerivatives = gpNaturalDerivativesPair.Value;
				shapeGradientsAtGPs[gaussPoint] =
					new EvalShapeGradients3D(naturalShapeDerivatives, new Jacobian3D(nodes, naturalShapeDerivatives));
			}

			return shapeGradientsAtGPs;
		}

		/// <summary>
		/// Transforms the coordinates from the natural (element local) coordinate system to the the global
		/// coordinate system of a point that is internal to the finite element.
		/// </summary>
		/// <param name="nodes">The coordinates of the finite element's nodes in the global cartesian system.</param>
		/// <param name="naturalPoint">The coordinates in the natural system of a point that is internal to the finite 
		///     element.</param>
		/// <returns></returns>
		public CartesianPoint3D TransformNaturalToCartesian(IReadOnlyList<Node3D> nodes, NaturalPoint3D naturalPoint)
		{
			double[] shapeFunctionValues = EvaluateAt(naturalPoint.Xi, naturalPoint.Eta, naturalPoint.Zeta);
			double x = 0, y = 0, z = 0;
			for (int i = 0; i < nodes.Count; i++)
			{
				x += shapeFunctionValues[i] * nodes[i].X;
				y += shapeFunctionValues[i] * nodes[i].Y;
				z += shapeFunctionValues[i] * nodes[i].Z;
			}

			return new CartesianPoint3D(x, y, z);
		}

		/// <summary>
		/// Evaluate shape function at a given point expressed in the natural coordinate system. Each entry corresponds to a
		/// different shape function.
		/// </summary>
		/// <param name="xi">The coordinate of the point along local axis Xi.</param>
		/// <param name="eta">The coordinate of the point along local axis Eta.</param>
		/// <param name="zeta">The coordinate of the point along local axis Zeta.</param>
		/// <returns></returns>
		protected abstract double[] EvaluateAt(double xi, double eta, double zeta);

		/// <summary>
		/// Evaluate derivatives of shape functions with repsect to natural coordinates at a given point expressed in the 
		/// natural coordinate system. Each row corresponds to a different shape function, column 0 corresponds to derivatives
		/// with respect to Xi coordinate and column 1 corresponds to derivatives with respect to Eta coordinate.
		/// </summary>
		/// <param name="xi">The coordinate of the point along local axis Xi.</param>
		/// <param name="eta">The coordinate of the point along local axis Eta.</param>
		/// <param name="zeta">The coordinate of the point along local axis Zeta.</param>
		/// <returns></returns>
		protected abstract double[,] EvaluateGradientsAt(double xi, double eta, double zeta);

		private Dictionary<GaussPoint3D, double[,]> EvaluateNaturalShapeDerivativesAtGPs(IQuadrature3D quadrature)
		{
			bool isCached = cachedNaturalShapederivativesAtGPs.TryGetValue(quadrature,
				out Dictionary<GaussPoint3D, double[,]> naturalShapeDerivativesAtGPs);
			if (!isCached)
			{
				naturalShapeDerivativesAtGPs= new Dictionary<GaussPoint3D, double[,]>();
				foreach (var gaussPoint in quadrature.IntegrationPoints)
				{
					naturalShapeDerivativesAtGPs[gaussPoint] = EvaluateGradientsAt(gaussPoint.Xi, gaussPoint.Eta, gaussPoint.Zeta);
				}
				cachedNaturalShapederivativesAtGPs.Add(quadrature, naturalShapeDerivativesAtGPs);
			}

			return naturalShapeDerivativesAtGPs;
		}

		private Dictionary<GaussPoint3D, double[]> EvaluateShapeFunctionsAtGPs(IQuadrature3D quadrature)
		{
			bool isCached =
				cachedRawShapeFunctionsAtGPs.TryGetValue(quadrature, out Dictionary<GaussPoint3D, double[]> shapeFuctionsAtGPs);
			if (!isCached)
			{
				shapeFuctionsAtGPs= new Dictionary<GaussPoint3D, double[]>();
				foreach (var gaussPoint in quadrature.IntegrationPoints)
				{
					shapeFuctionsAtGPs[gaussPoint] = EvaluateAt(gaussPoint.Xi, gaussPoint.Eta, gaussPoint.Zeta);
				}
				cachedRawShapeFunctionsAtGPs.Add(quadrature,shapeFuctionsAtGPs);
			}

			return shapeFuctionsAtGPs;
		}
	}
}
