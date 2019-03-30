using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Integration.Points;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation.Inverse;
using ISAAR.MSolve.FEM.Interpolation.Jacobians;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.FEM.Interpolation
{
    /// <summary>
    /// Basic implementation of <see cref="IIsoparametricInterpolation2D"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public abstract class IsoparametricInterpolation2DBase : IIsoparametricInterpolation2D
    {
        private readonly Dictionary<IQuadrature2D, IReadOnlyList<double[]>> cachedFunctionsAtGPs;
        private readonly Dictionary<IQuadrature2D, IReadOnlyList<Matrix>> cachedNaturalGradientsAtGPs;

        protected IsoparametricInterpolation2DBase(CellType cellType, int numFunctions)
        {
            this.CellType = cellType;
            this.NumFunctions = numFunctions;
            this.cachedFunctionsAtGPs = new Dictionary<IQuadrature2D, IReadOnlyList<double[]>>();
            this.cachedNaturalGradientsAtGPs = new Dictionary<IQuadrature2D, IReadOnlyList<Matrix>>();
        }

        /// <summary>
        /// See <see cref="IIsoparametricInterpolation2D.CellType"/>.
        /// </summary>
        public CellType CellType { get; }

        /// <summary>
        /// See <see cref="IIsoparametricInterpolation2D.NumFunctions"/>.
        /// </summary>
        public int NumFunctions { get; }

        /// <summary>
        /// See <see cref="IIsoparametricInterpolation2D.NodalNaturalCoordinates"/>.
        /// </summary>
        public abstract IReadOnlyList<NaturalPoint2D> NodalNaturalCoordinates { get; }

        /// <summary>
        /// See <see cref="IIsoparametricInterpolation2D.CreateInverseMappingFor(IReadOnlyList{Node_v2})"/>.
        /// </summary>
        public abstract IInverseInterpolation2D CreateInverseMappingFor(IReadOnlyList<Node_v2> nodes);

        /// <summary>
        /// See <see cref="IIsoparametricInterpolation2D.EvaluateAllAt(IReadOnlyList{Node_v2}, NaturalPoint2D)"/>.
        /// </summary>
        public EvalInterpolation2D EvaluateAllAt(IReadOnlyList<Node_v2> nodes, NaturalPoint2D naturalPoint)
        {
            double xi = naturalPoint.Xi;
            double eta = naturalPoint.Eta;
            var shapeFunctions = EvaluateAt(xi, eta);
            Matrix naturalShapeDerivatives = EvaluateGradientsAt(xi, eta);
            return new EvalInterpolation2D(shapeFunctions, naturalShapeDerivatives,
                new IsoparametricJacobian2D(nodes, naturalShapeDerivatives));
        }

        /// <summary>
        /// See <see cref="IIsoparametricInterpolation2D.EvaluateAllAtGaussPoints(IReadOnlyList{Node_v2}, IQuadrature2D)"/>.
        /// </summary>
        public IReadOnlyList<EvalInterpolation2D> EvaluateAllAtGaussPoints(IReadOnlyList<Node_v2> nodes, IQuadrature2D quadrature)
        {
            // The shape functions and natural derivatives at each Gauss point are probably cached from previous calls
            IReadOnlyList<double[]> shapeFunctionsAtGPs = EvaluateFunctionsAtGaussPoints(quadrature);
            IReadOnlyList<Matrix> naturalShapeDerivativesAtGPs = EvaluateNaturalGradientsAtGaussPoints(quadrature);

            // Calculate the Jacobians and shape derivatives w.r.t. global cartesian coordinates at each Gauss point
            int numGPs = quadrature.IntegrationPoints.Count;
            var interpolationsAtGPs = new EvalInterpolation2D[numGPs];
            for (int gp = 0; gp < numGPs; ++gp)
            {
                interpolationsAtGPs[gp] = new EvalInterpolation2D(shapeFunctionsAtGPs[gp],
                    naturalShapeDerivativesAtGPs[gp], new IsoparametricJacobian2D(nodes, naturalShapeDerivativesAtGPs[gp]));
            }
            return interpolationsAtGPs;
        }

        /// <summary>
        /// See <see cref="IIsoparametricInterpolation2D.EvaluateFunctionsAt(NaturalPoint2D)"/>.
        /// </summary>
        public double[] EvaluateFunctionsAt(NaturalPoint2D naturalPoint)
            => EvaluateAt(naturalPoint.Xi, naturalPoint.Eta);

        /// <summary>
        /// See <see cref="IIsoparametricInterpolation2D.EvaluateFunctionsAtGaussPoints(IQuadrature2D)"/>.
        /// </summary>
        public IReadOnlyList<double[]> EvaluateFunctionsAtGaussPoints(IQuadrature2D quadrature)
        {
            bool isCached = cachedFunctionsAtGPs.TryGetValue(quadrature,
                out IReadOnlyList<double[]> shapeFunctionsAtGPs);
            if (isCached) return shapeFunctionsAtGPs;
            else
            {
                int numGPs = quadrature.IntegrationPoints.Count;
                var shapeFunctionsAtGPsArray = new double[numGPs][];
                for (int gp = 0; gp < numGPs; ++gp)
                {
                    GaussPoint2D gaussPoint = quadrature.IntegrationPoints[gp];
                    shapeFunctionsAtGPsArray[gp] = EvaluateAt(gaussPoint.Xi, gaussPoint.Eta);
                }
                cachedFunctionsAtGPs.Add(quadrature, shapeFunctionsAtGPsArray);
                return shapeFunctionsAtGPsArray;
            }
        }

        /// <summary>
        /// See <see cref="IIsoparametricInterpolation2D.EvaluateNaturalGradientsAt(NaturalPoint2D)".
        /// </summary>
        public Matrix EvaluateNaturalGradientsAt(NaturalPoint2D naturalPoint)
            => EvaluateGradientsAt(naturalPoint.Xi, naturalPoint.Eta);

        /// <summary>
        /// See <see cref="IIsoparametricInterpolation2D.EvaluateNaturalGradientsAtGaussPoints(IQuadrature2D)"/>.
        /// </summary>
        /// <param name="quadrature"></param>
        public IReadOnlyList<Matrix> EvaluateNaturalGradientsAtGaussPoints(IQuadrature2D quadrature)
        {
            bool isCached = cachedNaturalGradientsAtGPs.TryGetValue(quadrature,
                out IReadOnlyList<Matrix> naturalGradientsAtGPs);
            if (isCached) return naturalGradientsAtGPs;
            else
            {
                int numGPs = quadrature.IntegrationPoints.Count;
                var naturalGradientsAtGPsArray = new Matrix[numGPs];
                for (int gp = 0; gp < numGPs; ++gp)
                {
                    GaussPoint2D gaussPoint = quadrature.IntegrationPoints[gp];
                    naturalGradientsAtGPsArray[gp] = EvaluateGradientsAt(gaussPoint.Xi, gaussPoint.Eta);
                }
                cachedNaturalGradientsAtGPs.Add(quadrature, naturalGradientsAtGPsArray);
                return naturalGradientsAtGPsArray;
            }
        }

        /// <summary>
        /// See <see cref="IIsoparametricInterpolation2D.TransformNaturalToCartesian(IReadOnlyList{Node_v2}, NaturalPoint2D)"/>.
        /// </summary>
        public CartesianPoint2D TransformNaturalToCartesian(IReadOnlyList<Node_v2> nodes, NaturalPoint2D naturalPoint)
        {
            double[] shapeFunctionValues = EvaluateAt(naturalPoint.Xi, naturalPoint.Eta);
            double x = 0, y = 0;
            for (int i = 0; i < nodes.Count; ++i)
            {
                x += shapeFunctionValues[i] * nodes[i].X;
                y += shapeFunctionValues[i] * nodes[i].Y;
            }
            return new CartesianPoint2D(x, y);
        }

        /// <summary>
        /// Evaluate shape function at a given point expressed in the natural coordinate system. Each entry corresponds to a
        /// different shape function.
        /// </summary>
        /// <param name="xi">The coordinate of the point along local axis Xi.</param>
        /// <param name="eta">The coordinate of the point along local axis Eta.</param>
        protected abstract double[] EvaluateAt(double xi, double eta);

        /// <summary>
        /// Evaluate derivatives of shape functions with respect to natural coordinates at a given point expressed in the 
        /// natural coordinate system. Each row corresponds to a different shape function, column 0 corresponds to derivatives
        /// with respect to Xi coordinate and column 1 corresponds to derivatives with respect to Eta coordinate.
        /// </summary>
        /// <param name="xi">The coordinate of the point along local axis Xi.</param>
        /// <param name="eta">The coordinate of the point along local axis Eta.</param>
        protected abstract Matrix EvaluateGradientsAt(double xi, double eta);
    }
}
