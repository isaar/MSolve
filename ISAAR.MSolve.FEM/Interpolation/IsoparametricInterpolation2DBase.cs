using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Integration.Points;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation.Inverse;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;

namespace ISAAR.MSolve.FEM.Interpolation
{
    /// <summary>
    /// Basic implementation of <see cref="IIsoparametricInterpolation2D"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public abstract class IsoparametricInterpolation2DBase: IIsoparametricInterpolation2D
    {
        private readonly Dictionary<IQuadrature2D, Dictionary<GaussPoint2D, double[,]>> cachedNaturalShapeDerivativesAtGPs;
        private readonly Dictionary<IQuadrature2D, Dictionary<GaussPoint2D, double[]>> cachedRawShapeFunctionsAtGPs;
        private readonly Dictionary<IQuadrature2D, Dictionary<GaussPoint2D, EvalShapeFunctions2D>> cachedShapeFunctionsAtGPs;

        protected IsoparametricInterpolation2DBase(CellType2D cellType, int numFunctions)
        {
            this.CellType = cellType;
            this.NumFunctions = numFunctions;
            this.cachedNaturalShapeDerivativesAtGPs = new Dictionary<IQuadrature2D, Dictionary<GaussPoint2D, double[,]>>();
            this.cachedRawShapeFunctionsAtGPs = new Dictionary<IQuadrature2D, Dictionary<GaussPoint2D, double[]>>();
            this.cachedShapeFunctionsAtGPs = new Dictionary<IQuadrature2D, Dictionary<GaussPoint2D, EvalShapeFunctions2D>>();
        }

        /// <summary>
        /// The shape of a cell. Useful for interacting with other modules and software.
        /// </summary>
        public CellType2D CellType { get; }

        /// <summary>
        /// The number of shape functions that define this interpolation.
        /// </summary>
        public int NumFunctions { get; }

        /// <summary>
        /// The coordinates of the finite element's nodes in the natural (element local) coordinate system. The order of these
        /// nodes matches the order of the shape functions and is always the same for each element.
        /// </summary>
        public abstract IReadOnlyList<NaturalPoint2D> NodalNaturalCoordinates { get; }

        /// <summary>
        /// The inverse mapping of this interpolation, namely from global cartesian to natural (element local) coordinate system.
        /// </summary>
        /// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
        /// <returns></returns>
        public abstract IInverseInterpolation2D CreateInverseMappingFor(IReadOnlyList<Node2D> nodes);

        /// <summary>
        /// Evaluate the shape functions and 1st order derivatives with respect to the global cartesian coordinate system 
        /// at a given natural point.
        /// </summary>
        /// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
        /// <param name="naturalPoint">The coordinates of the point in the natural (element local) coordinate system of the 
        ///     element.</param>
        /// <returns></returns>
        public EvalInterpolation2D EvaluateAllAt(IReadOnlyList<Node2D> nodes, NaturalPoint2D naturalPoint)
        {
            double xi = naturalPoint.Xi;
            double eta = naturalPoint.Eta;
            double[,] naturalShapeDerivatives = EvaluateGradientsAt(xi, eta);
            return new EvalInterpolation2D(EvaluateAt(xi, eta), naturalShapeDerivatives, 
                new Jacobian2D(nodes, naturalShapeDerivatives));
        }

        /// <summary>
        /// Evaluate the shape functions and 1st order derivatives with respect to the global cartesian coordinate system 
        /// at the integration points defined by a given quadrature. This method caches all possible quantities from previous
        /// calls. Use it instead of <see cref="EvaluateAllAt(IReadOnlyList{Node2D}, NaturalPoint2D)"/>
        /// when the integration points of an element are the same across multiple elements or multiple iterations of a
        /// non linear or dynamic analysis (in the latter cases we could also cache the <see cref="EvalInterpolation2D"/> at 
        /// the integration points of each element).
        /// </summary>
        /// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
        /// <param name="quadrature">The integration rule that defines integration points where shape functions and derivatives  
        ///     are calculated. The integration points of this instance of <see cref="IQuadrature2D"/> are always the 
        ///     same.</param>
        /// <returns></returns>
        public Dictionary<GaussPoint2D, EvalInterpolation2D> EvaluateAllAtGaussPoints(IReadOnlyList<Node2D> nodes,
            IQuadrature2D quadrature)
        {
            // The shape functions and derivatives at each Gauss point are probably cached from previous calls of this method
            Dictionary<GaussPoint2D, double[]> shapeFunctionsAtGPs = EvaluateShapeFunctionsAtGPs(quadrature);
            Dictionary<GaussPoint2D, double[,]> naturalShapeDerivativesAtGPs = EvaluateNaturalShapeDerivativesAtGPs(quadrature);

            // Calculate the Jacobians and shape derivatives w.r.t. global cartesian coordinates at each Gauss point
            var interpolationsAtGPs = new Dictionary<GaussPoint2D, EvalInterpolation2D>();
            foreach (GaussPoint2D gaussPoint in quadrature.IntegrationPoints)
            {
                interpolationsAtGPs[gaussPoint] = new EvalInterpolation2D(shapeFunctionsAtGPs[gaussPoint],
                    naturalShapeDerivativesAtGPs[gaussPoint], new Jacobian2D(nodes, naturalShapeDerivativesAtGPs[gaussPoint]));
            }
            return interpolationsAtGPs;
        }

        /// <summary>
        /// Evaluate the shape functions at a given natural point.
        /// </summary>
        /// <param name="naturalPoint">The coordinates of the point in the natural (element local) coordinate system of the 
        ///     element.</param>
        /// <returns></returns>
        public EvalShapeFunctions2D EvaluateFunctionsAt(NaturalPoint2D naturalPoint)
        {
            return new EvalShapeFunctions2D(EvaluateAt(naturalPoint.Xi, naturalPoint.Eta));
        }

        /// <summary>
        /// Evaluate the shape functions at the integration points defined by a given quadrature. This method caches all possible
        /// quantities from previous calls. Use it instead of <see cref="EvaluateAllAt(IReadOnlyList{Node2D}, NaturalPoint2D)"/>
        /// when the integration points of an element are the same across multiple elements or multiple iterations of a
        /// non linear or dynamic analysis.
        /// </summary>
        /// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
        /// <param name="quadrature">The integration rule that defines integration points where shape functions are calculated. 
        ///     The integration points of this instance of <see cref="IQuadrature2D"/> are always the same.</param>
        /// <returns></returns>
        public Dictionary<GaussPoint2D, EvalShapeFunctions2D> EvaluateFunctionsAtGaussPoints(IReadOnlyList<Node2D> nodes,
            IQuadrature2D quadrature)
        {
            bool isCached = cachedShapeFunctionsAtGPs.TryGetValue(quadrature,
                out Dictionary<GaussPoint2D, EvalShapeFunctions2D> shapeFunctionsAtGPs);
            if (!isCached)
            {
                shapeFunctionsAtGPs = new Dictionary<GaussPoint2D, EvalShapeFunctions2D>();
                foreach (GaussPoint2D gaussPoint in quadrature.IntegrationPoints)
                {
                    shapeFunctionsAtGPs[gaussPoint] = new EvalShapeFunctions2D(EvaluateAt(gaussPoint.Xi, gaussPoint.Eta));
                }
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
        public EvalShapeGradients2D EvaluateGradientsAt(IReadOnlyList<Node2D> nodes, NaturalPoint2D naturalPoint)
        {
            double[,] naturalShapeDerivatives = EvaluateGradientsAt(naturalPoint.Xi, naturalPoint.Eta);
            return new EvalShapeGradients2D(naturalShapeDerivatives, new Jacobian2D(nodes, naturalShapeDerivatives));
        }

        /// <summary>
        /// Evaluate the 1st order shape function derivatives with respect to the global cartesian coordinate system 
        /// at the integration points defined by a given quadrature. This method caches all possible quantities from previous
        /// calls. Use it instead of <see cref="EvaluateGradientsAt(IReadOnlyList{Node2D}, NaturalPoint2D)"/> 
        /// when the integration points of an element are the same across multiple elements or multiple iterations of a
        /// non linear or dynamic analysis (in the latter cases we could also cache the <see cref="EvalShapeGradients2D"/> at 
        /// the integration points of each element).
        /// </summary>
        /// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
        /// <param name="quadrature">The integration rule that defines integration points where shape function derivatives are 
        ///     calculated. The integration points of this instance of <see cref="IQuadrature2D"/> are always the same.</param>
        /// <returns></returns>
        public Dictionary<GaussPoint2D, EvalShapeGradients2D> EvaluateGradientsAtGaussPoints(IReadOnlyList<Node2D> nodes, 
            IQuadrature2D quadrature)
        {
            // The shape function derivatives at each Gauss point are probably cached from previous calls of this method
            Dictionary<GaussPoint2D, double[,]> naturalShapeDerivativesAtGPs = EvaluateNaturalShapeDerivativesAtGPs(quadrature);

            // Calculate the Jacobians and shape derivatives w.r.t. global cartesian coordinates at each Gauss point
            var shapeGradientsAtGPs = new Dictionary<GaussPoint2D, EvalShapeGradients2D>();
            foreach (var gpNaturalDerivativesPair in naturalShapeDerivativesAtGPs)
            {
                GaussPoint2D gaussPoint = gpNaturalDerivativesPair.Key;
                double[,] naturalShapeDerivatives = gpNaturalDerivativesPair.Value;
                shapeGradientsAtGPs[gaussPoint] = new EvalShapeGradients2D(naturalShapeDerivatives,
                     new Jacobian2D(nodes, naturalShapeDerivatives));
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
        public CartesianPoint2D TransformNaturalToCartesian(IReadOnlyList<Node2D> nodes, NaturalPoint2D naturalPoint)
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
        /// <returns></returns>
        protected abstract double[] EvaluateAt(double xi, double eta);

        /// <summary>
        /// Evaluate derivatives of shape functions with repsect to natural coordinates at a given point expressed in the 
        /// natural coordinate system. Each row corresponds to a different shape function, column 0 corresponds to derivatives
        /// with respect to Xi coordinate and column 1 corresponds to derivatives with respect to Eta coordinate.
        /// </summary>
        /// <param name="xi">The coordinate of the point along local axis Xi.</param>
        /// <param name="eta">The coordinate of the point along local axis Eta.</param>
        /// <returns></returns>
        protected abstract double[,] EvaluateGradientsAt(double xi, double eta);

        private Dictionary<GaussPoint2D, double[,]> EvaluateNaturalShapeDerivativesAtGPs(IQuadrature2D quadrature)
        {
            bool isCached = cachedNaturalShapeDerivativesAtGPs.TryGetValue(quadrature,
                out Dictionary<GaussPoint2D, double[,]> naturalShapeDerivativesAtGPs);
            if (!isCached)
            {
                naturalShapeDerivativesAtGPs = new Dictionary<GaussPoint2D, double[,]>();
                foreach (GaussPoint2D gaussPoint in quadrature.IntegrationPoints)
                {
                    naturalShapeDerivativesAtGPs[gaussPoint] = EvaluateGradientsAt(gaussPoint.Xi, gaussPoint.Eta);
                }
            }
            return naturalShapeDerivativesAtGPs;
        }

        private Dictionary<GaussPoint2D, double[]> EvaluateShapeFunctionsAtGPs(IQuadrature2D quadrature)
        {
            bool isCached = cachedRawShapeFunctionsAtGPs.TryGetValue(quadrature,
                out Dictionary<GaussPoint2D, double[]> shapeFunctionsAtGPs);
            if (!isCached)
            {
                shapeFunctionsAtGPs = new Dictionary<GaussPoint2D, double[]>();
                foreach (GaussPoint2D gaussPoint in quadrature.IntegrationPoints)
                {
                    shapeFunctionsAtGPs[gaussPoint] = EvaluateAt(gaussPoint.Xi, gaussPoint.Eta);
                }
            }
            return shapeFunctionsAtGPs;
        }
    }
}
