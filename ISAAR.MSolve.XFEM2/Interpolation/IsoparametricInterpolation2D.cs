using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Interpolation.InverseMappings;

namespace ISAAR.MSolve.XFEM.Interpolation
{
    abstract class IsoparametricInterpolation2D
    {
        #region enum values
        public static readonly IsoparametricInterpolation2D Quad4 = new Quad4ShapeFunctions();
        public static readonly IsoparametricInterpolation2D Quad9 = new Quad9ShapeFunctions();
        public static readonly IsoparametricInterpolation2D Tri3 = new Tri3ShapeFunctions();
        #endregion

        #region implemented members
        // Prevents non-nested derived classes. What happens if the class gets too big?
        // i) Use an interface and many similar (and interchangeable) enum classes (e.g. 1 for quads, 1 for triangles).
        // ii) Relax the requirement that derived classes must be nested. 
        // iii) Keep the nested requirement but use partial class
        private IsoparametricInterpolation2D(int functionsCount) { this.FunctionsCount = functionsCount; }

        public int FunctionsCount { get; }

        public EvaluatedInterpolation2D EvaluateAt(IReadOnlyList<Node2D> nodes, INaturalPoint2D naturalPoint)
        {
            double xi = naturalPoint.Xi;
            double eta = naturalPoint.Eta;
            double[,] naturalDerivatives = EvaluateDerivativesAt(xi, eta);
            // TODO: perhaps check that the nodes match the values and derivatives
            return new EvaluatedInterpolation2D(nodes, EvaluateAt(xi, eta),
                naturalDerivatives, new Jacobian2D(nodes, naturalDerivatives));
        }

        public EvaluatedInterpolation2D EvaluateOnlyDerivativesAt(IReadOnlyList<Node2D> nodes, INaturalPoint2D naturalPoint)
        {
            double[,] naturalDerivatives = EvaluateDerivativesAt(naturalPoint.Xi, naturalPoint.Eta);
            // TODO: perhaps check that the nodes match the values and derivatives
            return new EvaluatedInterpolation2D(nodes, naturalDerivatives, new Jacobian2D(nodes, naturalDerivatives));
        }

        public ICartesianPoint2D TransformNaturalToCartesian(IReadOnlyList<Node2D> nodes, INaturalPoint2D naturalPoint)
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
        #endregion

        #region abstract members
        protected abstract double[] EvaluateAt(double xi, double eta);
        protected abstract double[,] EvaluateDerivativesAt(double xi, double eta);
        public abstract IInverseMapping2D CreateInverseMappingFor(IReadOnlyList<Node2D> nodes);
        #endregion

        #region concrete (private) classes
        private class Quad4ShapeFunctions : IsoparametricInterpolation2D
        {
            public Quad4ShapeFunctions() : base(4) { }

            protected override sealed double[] EvaluateAt(double xi, double eta)
            {
                double[] values = new double[4];
                values[0] = 0.25 * (1 - xi) * (1 - eta);
                values[1] = 0.25 * (1 + xi) * (1 - eta);
                values[2] = 0.25 * (1 + xi) * (1 + eta);
                values[3] = 0.25 * (1 - xi) * (1 + eta);
                return values;
            }

            protected override sealed double[,] EvaluateDerivativesAt(double xi, double eta)
            {
                double[,] derivatives = new double[4, 2];
                derivatives[0, 0] = -0.25 * (1 - eta);
                derivatives[0, 1] = -0.25 * (1 - xi);
                derivatives[1, 0] = 0.25 * (1 - eta);
                derivatives[1, 1] = -0.25 * (1 + xi);
                derivatives[2, 0] = 0.25 * (1 + eta);
                derivatives[2, 1] = 0.25 * (1 + xi);
                derivatives[3, 0] = -0.25 * (1 + eta);
                derivatives[3, 1] = 0.25 * (1 - xi);
                return derivatives;
            }

            public override IInverseMapping2D CreateInverseMappingFor(IReadOnlyList<Node2D> nodes)
            {
                return new InverseQuad4Mapping(nodes);
            }
        }

        private class Quad9ShapeFunctions : IsoparametricInterpolation2D
        {
            public Quad9ShapeFunctions() : base(9) { }

            protected override sealed double[] EvaluateAt(double xi, double eta)
            {
                double xiEtaOver4 = 0.25 * xi * eta;
                double xi_2 = xi * xi;
                double eta_2 = eta * eta;

                double[] values = new double[9];
                values[0] = xiEtaOver4 * (1 - xi) * (1 - eta);
                values[1] = -xiEtaOver4 * (1 + xi) * (1 - eta);
                values[2] = xiEtaOver4 * (1 + xi) * (1 + eta);
                values[3] = -xiEtaOver4 * (1 - xi) * (1 + eta);
                values[4] = -0.5 * eta * (1 - xi_2) * (1 - eta);
                values[5] = 0.5 * xi * (1 + xi) * (1 - eta_2);
                values[6] = 0.5 * eta * (1 - xi_2) * (1 + eta);
                values[7] = -0.5 * xi * (1 - xi) * (1 - eta_2);
                values[8] = (1 - xi_2) * (1 - eta_2);
                return values;
            }

            protected override sealed double[,] EvaluateDerivativesAt(double xi, double eta)
            {
                double xi2 = xi * 2.0;
                double eta2 = eta * 2.0;
                double xi_2 = xi * xi;
                double eta_2 = eta * eta;
                double xiEta = xi * eta;

                double[,] derivatives = new double[9, 2];
                derivatives[0, 0] = 0.25 * eta * (1 - xi2) * (1 - eta);
                derivatives[0, 1] = 0.25 * xi * (1 - xi) * (1 - eta2);
                derivatives[1, 0] = -0.25 * eta * (1 + xi2) * (1 - eta);
                derivatives[1, 1] = -0.25 * xi * (1 + xi) * (1 - eta2);
                derivatives[2, 0] = 0.25 * eta * (1 + xi2) * (1 + eta);
                derivatives[2, 1] = 0.25 * xi * (1 + xi) * (1 + eta2);
                derivatives[3, 0] = -0.25 * eta * (1 - xi2) * (1 + eta);
                derivatives[3, 1] = 0.25 * xi * (1 - xi) * (1 + eta2);
                derivatives[4, 1] = xiEta * (1 - eta);
                derivatives[4, 2] = -0.5 * (1 - xi_2) * (1 - eta2);
                derivatives[5, 1] = 0.5 * (1 + xi2) * (1 - eta_2);
                derivatives[5, 2] = -xiEta * (1 + xi);
                derivatives[6, 1] = -xiEta * (1 + eta);
                derivatives[6, 2] = 0.5 * (1 - xi_2) * (1 + eta2);
                derivatives[7, 1] = -0.5 * (1 - xi2) * (1 - eta_2);
                derivatives[7, 2] = xiEta * (1 - xi);
                derivatives[8, 1] = -2 * xi * (1 - eta_2);
                derivatives[8, 2] = -2 * eta * (1 - xi_2);
                return derivatives;
            }

            public override IInverseMapping2D CreateInverseMappingFor(IReadOnlyList<Node2D> nodes)
            {
                throw new NotImplementedException("Probably requires an iterative procedure");
            }
        }

        // TODO: This needs to be more efficient. The design needs to change and some of the public methods become 
        // abstract. In triangles, the natural shape function derivatives are always the same. Also the jacobian is 
        // constant for the same element.
        private class Tri3ShapeFunctions : IsoparametricInterpolation2D
        {
            public Tri3ShapeFunctions() : base(3) { }

            protected override sealed double[] EvaluateAt(double xi, double eta)
            {
                double[] values = new double[3];
                values[0] = 1 - xi - eta;
                values[1] = xi;
                values[2] = eta;
                return values;
            }

            protected override sealed double[,] EvaluateDerivativesAt(double xi, double eta)
            {
                double[,] derivatives = new double[3, 2];
                derivatives[0, 0] = -1;
                derivatives[0, 1] = -1;
                derivatives[1, 0] = 1;
                derivatives[1, 1] = 0;
                derivatives[2, 0] = 0;
                derivatives[2, 1] = 1;
                return derivatives;
            }

            public override IInverseMapping2D CreateInverseMappingFor(IReadOnlyList<Node2D> nodes)
            {
                return new InverseTri3Mapping(nodes);
            }
        }
        #endregion
    }
}
