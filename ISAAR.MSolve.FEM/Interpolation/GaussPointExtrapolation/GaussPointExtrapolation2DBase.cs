using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Geometry.Coordinates;

//TODO: Perhaps I should use dictionaries for matching input and output at nodes or at least gauss points.
//TODO: Perhaps I should use generic methods when the same thing is done for scalar, vector, tensor fields. However that would 
//      complicate the interface only to avoid a few lines of boilerplate code.
//TODO: prevent triangular quadratures from accessing quadrilateral nodes and vice-versa
//TODO: cache the shape functions at the nodes.
namespace ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation
{
    /// <summary>
    /// Basic implementation of <see cref="IGaussPointExtrapolation2D"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public abstract class GaussPointExtrapolation2DBase : IGaussPointExtrapolation2D
    {
        /// <summary>
        /// Each <see cref="IIsoparametricInterpolation2D"/> is mapped to a 2D array that contains the values of the 
        /// extrapolation functions calculated at the nodes of the finite element using that interpolation. Each row corresponds
        /// to a different node. Each columns corresponds to a different extrapolation function.
        /// </summary>
        private Dictionary<IIsoparametricInterpolation2D, double[][]> cachedExtrapolationFunctionsAtNodes;

        protected GaussPointExtrapolation2DBase(IQuadrature2D quadrature)
        {
            this.cachedExtrapolationFunctionsAtNodes = new Dictionary<IIsoparametricInterpolation2D, double[][]>();
            this.Quadrature = quadrature;
        }

        /// <summary>
        /// The integration rule which defines the integration points used for extrapolating values and defining an auxiliary 
        /// coordinate system.
        /// </summary>
        public IQuadrature2D Quadrature { get; }

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
        {
            double[] extrapolationFunctions = EvaluateExtrapolationFunctionsAt(point);
            double scalar = 0;
            for (int gp = 0; gp < Quadrature.IntegrationPoints.Count; ++gp)
            {
                scalar += extrapolationFunctions[gp] * scalarsAtGaussPoints[gp];
            }
            return scalar;
        }

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
            double[][] nodalExtrapolationFunctions = EvaluateExtrapolationFunctionsAtNodes(interpolation);
            IReadOnlyList<NaturalPoint2D> nodes = interpolation.NodalNaturalCoordinates;
            var nodalScalars = new double[nodes.Count];
            for (int i = 0; i < nodes.Count; ++i)
            {
                //nodalScalars[i] = ExtrapolateScalarFromGaussPoints(scalarsAtGaussPoints, nodes[i]); // for debugging
                double scalar = 0;
                for (int gp = 0; gp < Quadrature.IntegrationPoints.Count; ++gp)
                {
                    scalar += nodalExtrapolationFunctions[i][gp] * scalarsAtGaussPoints[gp];
                }
                nodalScalars[i] = scalar;
            }
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
        {
            double[] extrapolationFunctions = EvaluateExtrapolationFunctionsAt(point);
            var tensor = new double[3]; //In 2D problems, symmetric tensors have 3 entries. TODO: replace with Tensor2D class.
            for (int gp = 0; gp < Quadrature.IntegrationPoints.Count; ++gp)
            {
                tensor[0] += extrapolationFunctions[gp] * tensorsAtGaussPoints[gp][0];
                tensor[1] += extrapolationFunctions[gp] * tensorsAtGaussPoints[gp][1];
                tensor[2] += extrapolationFunctions[gp] * tensorsAtGaussPoints[gp][2];
            }
            return tensor;
        }

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
            double[][] nodalExtrapolationFunctions = EvaluateExtrapolationFunctionsAtNodes(interpolation);
            IReadOnlyList<NaturalPoint2D> nodes = interpolation.NodalNaturalCoordinates;
            var nodalTensors = new double[nodes.Count][];
            for (int i = 0; i < nodes.Count; ++i)
            {
                //nodalTensors[i] = ExtrapolateVectorFromGaussPoints(tensorsAtGaussPoints, nodes[i]); // for debugging
                var tensor = new double[3];
                for (int gp = 0; gp < Quadrature.IntegrationPoints.Count; ++gp)
                {
                    tensor[0] += nodalExtrapolationFunctions[i][gp] * tensorsAtGaussPoints[gp][0];
                    tensor[1] += nodalExtrapolationFunctions[i][gp] * tensorsAtGaussPoints[gp][1];
                    tensor[2] += nodalExtrapolationFunctions[i][gp] * tensorsAtGaussPoints[gp][2];
                }
                nodalTensors[i] = tensor;
            }
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
        {
            double[] extrapolationFunctions = EvaluateExtrapolationFunctionsAt(point);
            var vector = new double[2]; //In 2D problems, vector fields have 2 entries. TODO: replace with Vector2 class.
            for (int gp = 0; gp < Quadrature.IntegrationPoints.Count; ++gp)
            {
                vector[0] += extrapolationFunctions[gp] * vectorsAtGaussPoints[gp][0];
                vector[1] += extrapolationFunctions[gp] * vectorsAtGaussPoints[gp][1];
            }
            return vector;
        }

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
            double[][] nodalExtrapolationFunctions = EvaluateExtrapolationFunctionsAtNodes(interpolation);
            IReadOnlyList<NaturalPoint2D> nodes = interpolation.NodalNaturalCoordinates;
            var nodalVectors = new double[nodes.Count][];
            for (int i = 0; i < nodes.Count; ++i)
            {
                //nodalVectors[i] = ExtrapolateVectorFromGaussPoints(tensorsAtGaussPoints, nodes[i]); // for debugging
                var vector = new double[2];
                for (int gp = 0; gp < Quadrature.IntegrationPoints.Count; ++gp)
                {
                    vector[0] += nodalExtrapolationFunctions[i][gp] * vectorsAtGaussPoints[gp][0];
                    vector[1] += nodalExtrapolationFunctions[i][gp] * vectorsAtGaussPoints[gp][1];
                }
                nodalVectors[i] = vector;
            }
            return nodalVectors;
        }

        /// <summary>
        /// Calculates the functions used for extrapolating quantities from the integration points to a given point, at the 
        /// given point.
        /// </summary>
        /// <param name="point">The coordinates of the point where the extrapolation functions will be calculated, in the 
        ///     natural (element local) system.</param>
        /// <returns></returns>
        protected abstract double[] EvaluateExtrapolationFunctionsAt(NaturalPoint2D point);

        private double[][] EvaluateExtrapolationFunctionsAtNodes(IIsoparametricInterpolation2D interpolation)
        {
            bool isCached = cachedExtrapolationFunctionsAtNodes.TryGetValue(interpolation,
                out double[][] nodalExtrapolationFunctions);
            if (!isCached)
            {
                IReadOnlyList<NaturalPoint2D> nodes = interpolation.NodalNaturalCoordinates;
                nodalExtrapolationFunctions = new double[nodes.Count][];
                for (int i = 0; i < nodes.Count; ++i)
                {
                    nodalExtrapolationFunctions[i] = EvaluateExtrapolationFunctionsAt(nodes[i]);
                }
                cachedExtrapolationFunctionsAtNodes.Add(interpolation, nodalExtrapolationFunctions);
            }
            return nodalExtrapolationFunctions;
        }
    }
}
