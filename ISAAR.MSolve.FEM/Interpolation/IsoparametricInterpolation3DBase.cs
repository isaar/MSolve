using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Integration.Points;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation.Inverse;
using ISAAR.MSolve.FEM.Interpolation.Jacobians;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.Numerical.LinearAlgebra;

namespace ISAAR.MSolve.FEM.Interpolation
{
    /// <summary>
    /// Basic implementation of <see cref="IIsoparametricInterpolation3D"/>
    /// Authors: Dimitris Tsapetis, Serafeim Bakalakos
    /// </summary>
    public abstract class IsoparametricInterpolation3DBase : IIsoparametricInterpolation3D
    {
        private readonly Dictionary<IQuadrature3D, Dictionary<GaussPoint3D, Vector>> cachedFunctionsAtGPs;
        private readonly Dictionary<IQuadrature3D, Dictionary<GaussPoint3D, Matrix2D>> cachedNaturalGradientsAtGPs;

        public IsoparametricInterpolation3DBase(int numFunctions)
        {
            this.NumFunctions = numFunctions;
            this.cachedFunctionsAtGPs = new Dictionary<IQuadrature3D, Dictionary<GaussPoint3D, Vector>>();
            this.cachedNaturalGradientsAtGPs = new Dictionary<IQuadrature3D, Dictionary<GaussPoint3D, Matrix2D>>();
        }

        /// <summary>
        /// See <see cref="IIsoparametricInterpolation3D.CellType"/>.
        /// </summary>
        public CellType3D CellType { get; }

        /// <summary>
        /// See <see cref="IIsoparametricInterpolation3D.NumFunctions"/>.
        /// </summary>
        public int NumFunctions { get; }

        /// <summary>
        /// See <see cref="IIsoparametricInterpolation3D.NodalNaturalCoordinates"/>.
        /// </summary>
        public abstract IReadOnlyList<NaturalPoint3D> NodalNaturalCoordinates { get; }

        /// <summary>
        /// See <see cref="IIsoparametricInterpolation3D.CreateInverseMappingFor(IReadOnlyList{Node3D})"/>.
        /// </summary>
        public abstract IInverseInterpolation3D CreateInverseMappingFor(IReadOnlyList<Node3D> nodes);

        /// <summary>
        /// See <see cref="IIsoparametricInterpolation3D.EvaluateAllAt(IReadOnlyList{Node3D}, NaturalPoint3D)"/>.
        /// </summary>
        public EvalInterpolation3D EvaluateAllAt(IReadOnlyList<Node3D> nodes, NaturalPoint3D naturalPoint)
        {
            double xi = naturalPoint.Xi;
            double eta = naturalPoint.Eta;
            double zeta = naturalPoint.Zeta;
            var shapeFunctions = new Vector(EvaluateAt(xi, eta, zeta));
            var naturalShapeDerivatives = new Matrix2D(EvaluateGradientsAt(xi, eta, zeta));
            return new EvalInterpolation3D(shapeFunctions, naturalShapeDerivatives,
                new IsoparametricJacobian3D(nodes, naturalShapeDerivatives));
        }

        /// <summary>
        /// See <see cref="IIsoparametricInterpolation3D.EvaluateAllAtGaussPoints(IReadOnlyList{Node3D}, IQuadrature3D)"/>
        /// </summary>
        public Dictionary<GaussPoint3D, EvalInterpolation3D> EvaluateAllAtGaussPoints(IReadOnlyList<Node3D> nodes,
            IQuadrature3D quadrature)
        {
            // The shape functions and natural derivatives at each Gauss point are probably cached from previous calls
            Dictionary<GaussPoint3D, Vector> shapeFunctionsAtGPs = EvaluateFunctionsAtGaussPoints(quadrature);
            Dictionary<GaussPoint3D, Matrix2D> naturalShapeDerivativesAtGPs = EvaluateNaturalGradientsAtGaussPoints(quadrature);

            // Calculate the Jacobians and shape derivatives w.r.t. global cartesian coordinates at each Gauss point
            var interpolationsAtGPs = new Dictionary<GaussPoint3D, EvalInterpolation3D>();
            foreach (var gaussPoint in quadrature.IntegrationPoints)
            {
                interpolationsAtGPs[gaussPoint] = new EvalInterpolation3D(shapeFunctionsAtGPs[gaussPoint],
                    naturalShapeDerivativesAtGPs[gaussPoint], new IsoparametricJacobian3D(nodes, naturalShapeDerivativesAtGPs[gaussPoint]));
            }
            return interpolationsAtGPs;
        }

        /// <summary>
        /// See <see cref="IIsoparametricInterpolation3D.EvaluateFunctionsAt(NaturalPoint3D)"/>.
        /// </summary>
        public Vector EvaluateFunctionsAt(NaturalPoint3D naturalPoint)
            => new Vector(EvaluateAt(naturalPoint.Xi, naturalPoint.Eta, naturalPoint.Zeta));

        /// <summary>
        /// See <see cref="IIsoparametricInterpolation3D.EvaluateFunctionsAtGaussPoints(IQuadrature3D)"/>.
        /// </summary>
        public Dictionary<GaussPoint3D, Vector> EvaluateFunctionsAtGaussPoints(IQuadrature3D quadrature)
        {
            bool isCached = cachedFunctionsAtGPs.TryGetValue(quadrature,
                out Dictionary<GaussPoint3D, Vector> shapeFunctionsAtGPs);
            if (!isCached)
            {
                shapeFunctionsAtGPs = new Dictionary<GaussPoint3D, Vector>();
                foreach (var gaussPoint in quadrature.IntegrationPoints)
                {
                    shapeFunctionsAtGPs[gaussPoint] = new Vector(EvaluateAt(gaussPoint.Xi, gaussPoint.Eta, gaussPoint.Zeta));
                }
                cachedFunctionsAtGPs.Add(quadrature, shapeFunctionsAtGPs);
            }
            return shapeFunctionsAtGPs;
        }

        /// <summary>
        /// See <see cref="IIsoparametricInterpolation3D.EvaluateNaturalGradientsAt(NaturalPoint3D)".
        /// </summary>
        public Matrix2D EvaluateNaturalGradientsAt(NaturalPoint3D naturalPoint)
            => new Matrix2D(EvaluateGradientsAt(naturalPoint.Xi, naturalPoint.Eta, naturalPoint.Zeta));

        /// <summary>
        /// See <see cref="IIsoparametricInterpolation3D.EvaluateNaturalGradientsAtGaussPoints(IQuadrature3D)"/>.
        /// </summary>
        /// <param name="quadrature"></param>
        public Dictionary<GaussPoint3D, Matrix2D> EvaluateNaturalGradientsAtGaussPoints(IQuadrature3D quadrature)
        {
            bool isCached = cachedNaturalGradientsAtGPs.TryGetValue(quadrature,
                out Dictionary<GaussPoint3D, Matrix2D> naturalGradientsAtGPs);
            if (!isCached)
            {
                naturalGradientsAtGPs = new Dictionary<GaussPoint3D, Matrix2D>();
                foreach (var gaussPoint in quadrature.IntegrationPoints)
                {
                    naturalGradientsAtGPs[gaussPoint] = new Matrix2D(
                        EvaluateGradientsAt(gaussPoint.Xi, gaussPoint.Eta, gaussPoint.Zeta));
                }
                cachedNaturalGradientsAtGPs.Add(quadrature, naturalGradientsAtGPs);
            }
            return naturalGradientsAtGPs;
        }

        /// <summary>
        /// See <see cref="IIsoparametricInterpolation3D.TransformNaturalToCartesian(IReadOnlyList{Node3D}, NaturalPoint3D)"/>.
        /// </summary>
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
        /// Evaluate derivatives of shape functions with respect to natural coordinates at a given point expressed in the 
        /// natural coordinate system. Each row corresponds to a different shape function, column 0 corresponds to derivatives
        /// with respect to Xi coordinate, column 1 corresponds to derivatives with respect to Eta coordinate, etc.
        /// </summary>
        /// <param name="xi">The coordinate of the point along local axis Xi.</param>
        /// <param name="eta">The coordinate of the point along local axis Eta.</param>
        /// <param name="zeta">The coordinate of the point along local axis Zeta.</param>
        /// <returns></returns>
        protected abstract double[,] EvaluateGradientsAt(double xi, double eta, double zeta);
    }
}
