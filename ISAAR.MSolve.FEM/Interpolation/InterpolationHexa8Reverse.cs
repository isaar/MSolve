using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation.Inverse;
using ISAAR.MSolve.Geometry.Coordinates;
using System;
using System.Collections.Generic;

namespace ISAAR.MSolve.FEM.Interpolation
{
    public class InterpolationHexa8Reverse: IsoparametricInterpolation3DBase
    {
        private static readonly InterpolationHexa8Reverse uniqueInstance = new InterpolationHexa8Reverse();

        private InterpolationHexa8Reverse() : base(8)
        {
            NodalNaturalCoordinates = new NaturalPoint3D[]
            {
                new NaturalPoint3D(1, 1, 1),
                new NaturalPoint3D(-1, 1, 1),
                new NaturalPoint3D(-1, -1, 1),
                new NaturalPoint3D(1, -1, 1),
                new NaturalPoint3D(1, 1, -1),
                new NaturalPoint3D(-1, 1, -1),
                new NaturalPoint3D(-1, -1, -1),
                new NaturalPoint3D(1, -1, -1)
            };
        }

        public override IReadOnlyList<NaturalPoint3D> NodalNaturalCoordinates { get; }

        public static InterpolationHexa8Reverse UniqueInstance => uniqueInstance;

        public override IInverseInterpolation3D CreateInverseMappingFor(IReadOnlyList<Node3D> nodes) => throw new NotImplementedException();

        protected override double[] EvaluateAt(double xi, double eta, double zeta)
        {
            throw new NotImplementedException(); //TODO: fill these
            //var shapeFunctions = new double[8];
            //shapeFunctions[0] = ;
            //shapeFunctions[1] = ;
            //shapeFunctions[2] = ;
            //shapeFunctions[3] = ;
            //shapeFunctions[4] = ;
            //shapeFunctions[5] = ;
            //shapeFunctions[6] = ;
            //shapeFunctions[7] = ;
        }

        protected override double[,] EvaluateGradientsAt(double xi, double eta, double zeta)
        {
            var naturalDerivatives = new double[3, 8]; //TODO: perhaps transpose this and the client code that uses it

            // Derivatives with respect to Xi
            naturalDerivatives[0, 0] = +0.125 * (1 + eta) * (1 + zeta);
            naturalDerivatives[0, 1] = -0.125 * (1 + eta) * (1 + zeta);
            naturalDerivatives[0, 2] = -0.125 * (1 - eta) * (1 + zeta);
            naturalDerivatives[0, 3] = +0.125 * (1 - eta) * (1 + zeta);
            naturalDerivatives[0, 4] = +0.125 * (1 + eta) * (1 - zeta);
            naturalDerivatives[0, 5] = -0.125 * (1 + eta) * (1 - zeta);
            naturalDerivatives[0, 6] = -0.125 * (1 - eta) * (1 - zeta);
            naturalDerivatives[0, 7] = +0.125 * (1 - eta) * (1 - zeta);

            // Derivatives with respect to Eta
            naturalDerivatives[1, 0] = 0.125 * (1 + xi) * (+1) * (1 + zeta);
            naturalDerivatives[1, 1] = 0.125 * (1 - xi) * (+1) * (1 + zeta);
            naturalDerivatives[1, 2] = 0.125 * (1 - xi) * (-1) * (1 + zeta);
            naturalDerivatives[1, 3] = 0.125 * (1 + xi) * (-1) * (1 + zeta);
            naturalDerivatives[1, 4] = 0.125 * (1 + xi) * (+1) * (1 - zeta);
            naturalDerivatives[1, 5] = 0.125 * (1 - xi) * (+1) * (1 - zeta);
            naturalDerivatives[1, 6] = 0.125 * (1 - xi) * (-1) * (1 - zeta);
            naturalDerivatives[1, 7] = 0.125 * (1 + xi) * (-1) * (1 - zeta);

            // Derivatives with respect to Zeta
            naturalDerivatives[2, 0] = 0.125 * (1 + xi) * (1 + eta) * (+1);
            naturalDerivatives[2, 1] = 0.125 * (1 - xi) * (1 + eta) * (+1);
            naturalDerivatives[2, 2] = 0.125 * (1 - xi) * (1 - eta) * (+1);
            naturalDerivatives[2, 3] = 0.125 * (1 + xi) * (1 - eta) * (+1);
            naturalDerivatives[2, 4] = 0.125 * (1 + xi) * (1 + eta) * (-1);
            naturalDerivatives[2, 5] = 0.125 * (1 - xi) * (1 + eta) * (-1);
            naturalDerivatives[2, 6] = 0.125 * (1 - xi) * (1 - eta) * (-1);
            naturalDerivatives[2, 7] = 0.125 * (1 + xi) * (1 - eta) * (-1);

            return naturalDerivatives;
        }
    }
}
