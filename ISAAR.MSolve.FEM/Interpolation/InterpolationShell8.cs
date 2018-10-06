using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation.Inverse;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Text;

//TODO: rename IsoparametricInterpolation3DBase. It works for shells as well.
namespace ISAAR.MSolve.FEM.Interpolation
{
    public class InterpolationShell8 : IsoparametricInterpolation3DBase
    {
        private static readonly InterpolationShell8 uniqueInstance = new InterpolationShell8();

        private InterpolationShell8() : base(8)
        {
            NodalNaturalCoordinates = new NaturalPoint3D[]
            {
                //TODO: validate this
                new NaturalPoint3D(1, 1, 0),
                new NaturalPoint3D(-1, 1, 0),
                new NaturalPoint3D(-1, -1, 0),
                new NaturalPoint3D(1, -1, 0),
                new NaturalPoint3D(0, 1, 0),
                new NaturalPoint3D(-1, 0, 0),
                new NaturalPoint3D(0, -1, 0),
                new NaturalPoint3D(1, 0, 0)
            };
        }

        public override IReadOnlyList<NaturalPoint3D> NodalNaturalCoordinates { get; }

        public static InterpolationShell8 UniqueInstance => uniqueInstance;

        public override IInverseInterpolation3D CreateInverseMappingFor(IReadOnlyList<Node3D> nodes) => throw new NotImplementedException();

        protected override double[] EvaluateAt(double xi, double eta, double zeta)
        {
            var shapeFunctions = new double[8];
            shapeFunctions[4] = 0.5 * (1 - Math.Pow(xi, 2)) * (1 + eta);
            shapeFunctions[5] = 0.5 * (1 - Math.Pow(eta, 2)) * (1 - xi);
            shapeFunctions[6] = 0.5 * (1 - Math.Pow(xi, 2)) * (1 - eta);
            shapeFunctions[7] = 0.5 * (1 - Math.Pow(eta, 2)) * (1 + xi);
            shapeFunctions[0] = 0.25 * (1 + xi) * (1 + eta) - 0.5 * shapeFunctions[4] - 0.5 * shapeFunctions[7];
            shapeFunctions[1] = 0.25 * (1 - xi) * (1 + eta) - 0.5 * shapeFunctions[4] - 0.5 * shapeFunctions[5];
            shapeFunctions[2] = 0.25 * (1 - xi) * (1 - eta) - 0.5 * shapeFunctions[5] - 0.5 * shapeFunctions[6];
            shapeFunctions[3] = 0.25 * (1 + xi) * (1 - eta) - 0.5 * shapeFunctions[6] - 0.5 * shapeFunctions[7];
            return shapeFunctions;
        }

        protected override double[,] EvaluateGradientsAt(double xi, double eta, double zeta)
        {
            var naturalDerivatives = new double[8, 2];

            // Derivatives with respect to Xi
            naturalDerivatives[4, 0] = (-xi) * (1 + eta);
            naturalDerivatives[5, 0] = -0.5 * (1 - Math.Pow(eta, 2));
            naturalDerivatives[6, 0] = 0.5 * (-2 * xi) * (1 - eta);
            naturalDerivatives[7, 0] = 0.5 * (1 - Math.Pow(eta, 2));
            naturalDerivatives[0, 0] = +0.25 * (1 + eta) - 0.5 * naturalDerivatives[4, 0] - 0.5 * naturalDerivatives[7, 0];
            naturalDerivatives[1, 0] = -0.25 * (1 + eta) - 0.5 * naturalDerivatives[4, 0] - 0.5 * naturalDerivatives[5, 0];
            naturalDerivatives[2, 0] = -0.25 * (1 - eta) - 0.5 * naturalDerivatives[5, 0] - 0.5 * naturalDerivatives[6, 0];
            naturalDerivatives[3, 0] = +0.25 * (1 - eta) - 0.5 * naturalDerivatives[6, 0] - 0.5 * naturalDerivatives[7, 0];

            // Derivatives with respect to Eta
            naturalDerivatives[4, 1] = 0.5 * (1 - Math.Pow(xi, 2));
            naturalDerivatives[5, 1] = 0.5 * (-2 * eta) * (1 - xi);
            naturalDerivatives[6, 1] = 0.5 * (1 - Math.Pow(xi, 2)) * (-1);
            naturalDerivatives[7, 1] = 0.5 * (-2 * eta) * (1 + xi);
            naturalDerivatives[0, 1] = +0.25 * (1 + xi) - 0.5 * naturalDerivatives[4, 1] - 0.5 * naturalDerivatives[7, 1];
            naturalDerivatives[1, 1] = +0.25 * (1 - xi) - 0.5 * naturalDerivatives[4, 1] - 0.5 * naturalDerivatives[5, 1];
            naturalDerivatives[2, 1] = -0.25 * (1 - xi) - 0.5 * naturalDerivatives[5, 1] - 0.5 * naturalDerivatives[6, 1];
            naturalDerivatives[3, 1] = -0.25 * (1 + xi) - 0.5 * naturalDerivatives[6, 1] - 0.5 * naturalDerivatives[7, 1];

            return naturalDerivatives;
        }
    }
}
