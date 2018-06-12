using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation.Inverse;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.FEM.Interpolation
{
    /// <summary>
    /// Basic implementation of <see cref="IIsoparametricInterpolation2D"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public abstract class IsoparametricInterpolation2DBase: IIsoparametricInterpolation2D
    {
        protected IsoparametricInterpolation2DBase(int numFunctions)
        {
            this.NumFunctions = numFunctions;
        }

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
            double[,] naturalDerivatives = EvaluateGradientsAt(xi, eta);
            return new EvalInterpolation2D(EvaluateAt(xi, eta), naturalDerivatives, new Jacobian2D(nodes, naturalDerivatives));
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
        /// Evaluate the 1st order shape function derivatives with respect to the global cartesian coordinate system 
        /// at a given natural point.
        /// </summary>
        /// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
        /// <param name="naturalPoint">The coordinates of the point in the natural (element local) coordinate system of the 
        ///     element.</param>
        /// <returns></returns>
        public EvalShapeGradients2D EvaluateGradientsAt(IReadOnlyList<Node2D> nodes, NaturalPoint2D naturalPoint)
        {
            double[,] naturalDerivatives = EvaluateGradientsAt(naturalPoint.Xi, naturalPoint.Eta);
            return new EvalShapeGradients2D(naturalDerivatives, new Jacobian2D(nodes, naturalDerivatives));
        }

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
    }
}
