using ISAAR.MSolve.Discretization.Integration.Points;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation;
using System;
using System.Collections.Generic;
using System.Text;
using Xunit;

namespace ISAAR.MSolve.FEM.Tests.Interpolation.Extrapolation
{
    /// <summary>
    /// Unit testing for <see cref="ExtrapolationGaussLegendre3x3"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class GaussLegendeExtrapolation3x3
    {
        [Fact]
        private static void ExtrapolateToGaussPoints()
        {
            double tolerance = 1e-10;
            var extrapolation = ExtrapolationGaussLegendre3x3.UniqueInstance;
            var gaussPoints = GaussLegendre2D.Order3x3.IntegrationPoints;

            for (int i = 0; i < gaussPoints.Count; ++i)
            {
                var scalarsAtGPs = new double[9];
                scalarsAtGPs[i] = 1.0;
                for (int j = 0; j < gaussPoints.Count; ++j)
                {
                    double extrapolatedScalar = extrapolation.ExtrapolateScalarFromGaussPoints(scalarsAtGPs, gaussPoints[j]);
                    if (j == i) Assert.True(Utilities.AreValuesEqual(1.0, extrapolatedScalar, tolerance));
                    else Assert.True(Utilities.AreValuesEqual(0.0, extrapolatedScalar, tolerance));
                }
            }
        }

        private static void ExtrapolateToNodes()
        {
            throw new NotImplementedException("I need to calclate these myself first.");
        }
    }
}
