using ISAAR.MSolve.IGA.Entities;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.IGA.Problems.SupportiveClasses
{
    public class GaussQuadrature
    {
        #region Coordinates and Weights
        //https://pomax.github.io/bezierinfo/legendre-gauss.html
        private readonly static double[][] coordinate = new double[][]
        {
            new double[]{ 0.0 },
            new double[]{ -0.5773502691896257, 0.5773502691896257 },
            new double[]{ -0.7745966692414836, 0.000000000000000, 0.7745966692414836 },
            new double[]{ -0.8611363115940526, -0.33998104358485626, 0.33998104358485626, 0.8611363115940526 },
            new double[]{ -0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640 },
            new double[]{ -0.9324695142031521, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969, 0.6612093864662645,0.9324695142031521 },
            new double[]{ -0.9491079123427585, -0.7415311855993945, -0.4058451513773972, 0.0, 0.4058451513773972,0.7415311855993945, 0.9491079123427585 },
            new double[]{ -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363 }
        };

        private readonly static double[][] weight = new double[][]
        {
            new double[]{ 2.0 },
            new double[]{ 1.0, 1.0 },
            new double[]{ 0.555555555555556, 0.888888888888889, 0.555555555555556 },
            new double[]{ 0.34785484513745385, 0.6521451548625461, 0.6521451548625461, 0.34785484513745385 },
            new double[]{ 0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891 },
            new double[]{ 0.1713244923791704, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910, 0.3607615730481386, 0.1713244923791704 },
            new double[]{ 0.1294849661688697, 0.2797053914892766, 0.3818300505051189, 0.4179591836734694, 0.3818300505051189, 0.2797053914892766, 0.1294849661688697 },
            new double[]{ 0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763 }
        };
        #endregion

        public IList<GaussLegendrePoint3D> CalculateElementGaussPoints(int degreeKsi, int degreeHeta, IList<Knot> knotsOfElement)
        {
            IList<GaussLegendrePoint3D> gaussPointsofElement = new List<GaussLegendrePoint3D>();
            int numberOfGPKsi = degreeKsi + 1;
            int numberOfGPHeta = degreeHeta + 1;
            double[] coordinatesKsi = new double[numberOfGPKsi];
            double[] weightKsi = new double[numberOfGPKsi];
            double[] coordinatesHeta = new double[numberOfGPHeta];
            double[] weightHeta = new double[numberOfGPHeta];

            for (int i = 0; i < numberOfGPKsi; i++)
            {
                coordinatesKsi[i] = 0.5 * (knotsOfElement[0].Ksi + knotsOfElement[2].Ksi
                    + (knotsOfElement[2].Ksi - knotsOfElement[0].Ksi) * coordinate[degreeKsi][i]);
                weightKsi[i] = 0.5 * ((knotsOfElement[2].Ksi-knotsOfElement[0].Ksi) * weight[degreeKsi][i]);
            }

            for (int i = 0; i < numberOfGPHeta; i++)
            {
                coordinatesHeta[i] = 0.5 * (knotsOfElement[0].Heta + knotsOfElement[1].Heta
                    + (knotsOfElement[1].Heta - knotsOfElement[0].Heta) * coordinate[degreeHeta][i]);
                weightHeta[i] = 0.5 * ((knotsOfElement[1].Heta - knotsOfElement[0].Heta) * weight[degreeHeta][i]);
            }
            
            for (int i = 0; i < numberOfGPKsi; i++)
            {
                for (int j = 0; j < numberOfGPHeta; j++)
                {
                    gaussPointsofElement.Add(new GaussLegendrePoint3D(coordinatesKsi[i], coordinatesHeta[j], 0, null, weightKsi[i] * weightHeta[j]));
                }

            }
            return gaussPointsofElement;
        }

        internal IList<GaussLegendrePoint3D> CalculateElementGaussPoints(int degreeKsi, IList<Knot> knotsOfElement)
        {
            IList<GaussLegendrePoint3D> gaussPointsofElement = new List<GaussLegendrePoint3D>();
            int numberOfGPKsi = degreeKsi + 1;
            double[] coordinatesKsi = new double[numberOfGPKsi];
            double[] weightKsi = new double[numberOfGPKsi];

            for (int i = 0; i < numberOfGPKsi; i++)
            {
                coordinatesKsi[i] = 0.5 * (knotsOfElement[0].Ksi + knotsOfElement[1].Ksi
                    + (knotsOfElement[1].Ksi - knotsOfElement[0].Ksi) * coordinate[degreeKsi][i]);
                weightKsi[i] = 0.5 * ((knotsOfElement[1].Ksi - knotsOfElement[0].Ksi) * weight[degreeKsi][i]);
            }
            
            for (int i = 0; i < numberOfGPKsi; i++)
            {
                gaussPointsofElement.Add(new GaussLegendrePoint3D(coordinatesKsi[i], 0, 0, null, weightKsi[i]));
            }
            return gaussPointsofElement;
        }

        public IList<GaussLegendrePoint3D> CalculateElementGaussPoints(int degreeKsi, int degreeHeta, int degreeZeta, IList<Knot> knotsOfElement)
        {
            IList<GaussLegendrePoint3D> gaussPointsofElement = new List<GaussLegendrePoint3D>();
            int numberOfGPKsi = degreeKsi + 1;
            int numberOfGPHeta = degreeHeta + 1;
            int numberOfGPZeta = degreeZeta + 1;
            double[] coordinatesKsi = new double[numberOfGPKsi];
            double[] weightKsi = new double[numberOfGPKsi];
            double[] coordinatesHeta = new double[numberOfGPHeta];
            double[] weightHeta = new double[numberOfGPHeta];
            double[] coordinatesZeta = new double[numberOfGPZeta];
            double[] weightZeta = new double[numberOfGPZeta];

            for (int i = 0; i < numberOfGPKsi; i++)
            {
                coordinatesKsi[i] = 0.5 * (knotsOfElement[0].Ksi + knotsOfElement[4].Ksi
                    + (knotsOfElement[4].Ksi - knotsOfElement[0].Ksi) * coordinate[degreeKsi][i]);
                weightKsi[i] = 0.5 * ((knotsOfElement[4].Ksi - knotsOfElement[0].Ksi) * weight[degreeKsi][i]);
            }

            for (int i = 0; i < numberOfGPHeta; i++)
            {
                coordinatesHeta[i] = 0.5 * (knotsOfElement[0].Heta + knotsOfElement[2].Heta
                    + (knotsOfElement[2].Heta - knotsOfElement[0].Heta) * coordinate[degreeHeta][i]);
                weightHeta[i] = 0.5 * ((knotsOfElement[2].Heta - knotsOfElement[0].Heta) * weight[degreeHeta][i]);
            }

            for (int i = 0; i < numberOfGPZeta; i++)
            {
                coordinatesZeta[i] = 0.5 * (knotsOfElement[0].Zeta + knotsOfElement[1].Zeta
                    + (knotsOfElement[1].Zeta - knotsOfElement[0].Zeta) * coordinate[degreeZeta][i]);
                weightZeta[i] = 0.5 * ((knotsOfElement[1].Zeta - knotsOfElement[0].Zeta) * weight[degreeZeta][i]);
            }


            for (int i = 0; i < numberOfGPKsi; i++)
            {
                for (int j = 0; j < numberOfGPHeta; j++)
                {
                    for (int k = 0; k < numberOfGPZeta; k++)
                    {
                        gaussPointsofElement.Add(new GaussLegendrePoint3D(coordinatesKsi[i], coordinatesHeta[j], coordinatesZeta[k], null, weightKsi[i] * weightHeta[j]* weightZeta[k]));
                    }                    
                }
            }
            return gaussPointsofElement;
        }

    }



    public class GaussLegendrePoint3D
    {
        public double [,] DeformationMatrix { get; private set; }
        public double Ksi { get; private set; }
        public double Heta { get; private set; }
        public double Zeta { get; private set; }
        public double WeightFactor { get; private set; }

        public GaussLegendrePoint3D(double ksi, double heta, double zeta, double [,] deformationMatrix, double weightFactor)
        {
            this.Ksi = ksi;
            this.Heta = heta;
            this.Zeta = zeta;
            this.WeightFactor = weightFactor;
            this.DeformationMatrix = deformationMatrix;
        }
    }
}
