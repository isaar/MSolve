using ISAAR.MSolve.Geometry.Coordinates;
using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Integration.Quadratures;

namespace ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation
{
    /// Calculates extrapolations of scalar, vector and tensor fields from the integration points of symmetric Gauss quadrature
    /// for triangles with 1 Gauss point. This can be done at any point, but utility methods for directly outputting the 
    /// extrapolated fields at the nodes of finite elements are also provided. Note that since there is only 1 Gauss point,
    /// the scalar, vector and tensor fields are constant at all points and equal to their values at the Gauss point.
    /// Implements Singleton pattern.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class ExtrapolationGaussTriangular1Point : IGaussPointExtrapolation2D
    {
        private static readonly ExtrapolationGaussTriangular1Point uniqueInstance = new ExtrapolationGaussTriangular1Point();

        /// <summary>
        /// The integration rule which defines the integration points used for extrapolating values and defining an auxiliary 
        /// coordinate system.
        /// </summary>
        public IQuadrature2D Quadrature { get { return TriangleQuadratureSymmetricGaussian.Order1Point1; } }

        /// <summary>
        /// Get the unique <see cref="ExtrapolationGaussTriangular1Point"/> object for the whole program. Thread safe.
        /// </summary>
        public static ExtrapolationGaussTriangular1Point UniqueInstance { get { return uniqueInstance; } }

        /// <summary>
        /// Calculates a scalar quantity at a given point by extrapolating (or interpolating) its known values at 
        /// the integration points.
        /// </summary>
        /// <param name="scalarsAtGaussPoints">Their order must be the same as the order of integration points defined by 
        ///     <see cref="Quadrature"/>.</param>
        /// <param name="point">The point where the scalar will be computed. Its coordinates are expressed in the natural
        ///     (element local) system, instead of the coordinate system defined by the integration points.</param>
        /// <returns></returns>
        public double ExtrapolateScalarFromGaussPoints(IReadOnlyList<double> scalarsAtGaussPoints, NaturalPoint2D point)
            => scalarsAtGaussPoints[0];

        /// <summary>
        /// Calculates a scalar quantity at the nodes of a finite element by extrapolating its known values at the integration 
        /// points.
        /// </summary>
        /// <param name="scalarsAtGaussPoints">Their order must be the same as the order of integration points defined by 
        ///     <see cref="Quadrature"/>.</param>
        /// <param name="interpolation">Defines the natural coordinates of the finite element's nodes.</param>
        /// <returns></returns>
        public IReadOnlyList<double> ExtrapolateScalarFromGaussPointsToNodes(IReadOnlyList<double> scalarsAtGaussPoints, 
            IIsoparametricInterpolation2D interpolation)
        {
            var nodalScalars = new double[interpolation.NumFunctions];
            for (int i = 0; i < nodalScalars.Length; ++i) nodalScalars[i] = scalarsAtGaussPoints[0];
            return nodalScalars;
        }

        /// <summary>
        /// Calculates a tensor quantity at a given point by extrapolating (or interpolating) its known values at 
        /// the integration points.
        /// </summary>
        /// <param name="tensorsAtGaussPoints">Their order must be the same as the order of integration points defined by 
        ///     <see cref="Quadrature"/>.</param>
        /// <param name="point">The point where the tensor will be computed. Its coordinates are expressed in the natural
        ///     (element local) system, instead of the coordinate system defined by the integration points.</param>
        /// <returns></returns>
        public double[] ExtrapolateTensorFromGaussPoints(IReadOnlyList<double[]> tensorsAtGaussPoints, NaturalPoint2D point)
            => tensorsAtGaussPoints[0];

        /// <summary>
        /// Calculates a tensor quantity at the nodes of a finite element by extrapolating its known values at the integration 
        /// points.
        /// </summary>
        /// <param name="tensorsAtGaussPoints">Their order must be the same as the order of integration points defined by 
        ///     <see cref="Quadrature"/>.</param>
        /// <param name="interpolation">Defines the natural coordinates of the finite element's nodes.</param>
        /// <returns></returns>
        public IReadOnlyList<double[]> ExtrapolateTensorFromGaussPointsToNodes(IReadOnlyList<double[]> tensorsAtGaussPoints, 
            IIsoparametricInterpolation2D interpolation)
        {
            var nodalTensors = new double[interpolation.NumFunctions][];
            for (int i = 0; i < nodalTensors.Length; ++i) nodalTensors[i] = tensorsAtGaussPoints[0];
            return nodalTensors;
        }

        /// <summary>
        /// Calculates a vector quantity at a given point by extrapolating (or interpolating) its known values at 
        /// the integration points.
        /// </summary>
        /// <param name="vectorsAtGaussPoints">Their order must be the same as the order of integration points defined by 
        ///     <see cref="Quadrature"/>.</param>
        /// <param name="point">The point where the tensor will be computed. Its coordinates are expressed in the natural
        ///     (element local) system, instead of the coordinate system defined by the integration points.</param>
        /// <returns></returns>
        public double[] ExtrapolateVectorFromGaussPoints(IReadOnlyList<double[]> vectorsAtGaussPoints, NaturalPoint2D point)
            => vectorsAtGaussPoints[0];

        /// <summary>
        /// Calculates a vector quantity at the nodes of a finite element by extrapolating its known values at the integration 
        /// points.
        /// </summary>
        /// <param name="vectorsAtGaussPoints">Their order must be the same as the order of integration points defined by 
        ///     <see cref="Quadrature"/>.</param>
        /// <param name="interpolation">Defines the natural coordinates of the finite element's nodes.</param>
        /// <returns></returns>
        public IReadOnlyList<double[]> ExtrapolateVectorFromGaussPointsToNodes(IReadOnlyList<double[]> vectorsAtGaussPoints, 
            IIsoparametricInterpolation2D interpolation)
        {
            var nodalVectors = new double[interpolation.NumFunctions][];
            for (int i = 0; i < nodalVectors.Length; ++i) nodalVectors[i] = vectorsAtGaussPoints[0];
            return nodalVectors;
        }
    }
}
