using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Globalization;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.FEM.Stochastic
{
    public class GaussianFileStochasticCoefficientsProvider : IStochasticMaterialCoefficientsProvider, IStochasticCoefficientsProvider
    {
        private readonly double[][,] evp;
        private readonly int[,] evpPos;
        private readonly double[] eigs;
        private readonly double[] translations;
        private readonly double mean;
        private readonly double covariance;
        private readonly double[] stochasticDomain;
        private static NumberFormatInfo ni = null;

        public int CurrentOrder { get; set; }
        public double[] RandomVariables { get; set; }

        public GaussianFileStochasticCoefficientsProvider(string fileName, int noOfVariables, char separator, double[] stochasticDomain, string decimalSeparator = ".")
        {
            CultureInfo ci = System.Globalization.CultureInfo.InstalledUICulture;
            ni = (System.Globalization.NumberFormatInfo)ci.NumberFormat.Clone();
            ni.NumberDecimalSeparator = decimalSeparator;
            this.stochasticDomain = stochasticDomain;

            int noOfCoords = stochasticDomain.Length;
            eigs = new double[noOfVariables];
            evpPos = new int[noOfVariables, 3];
            evp = new double[noOfCoords][,];
            translations = new double[noOfCoords];

            for (int i = 0; i < noOfCoords; i++)
            {
                translations[i] = stochasticDomain[i] / 2;
                evp[i] = new double[noOfVariables, 3];
            }

            string[] lines = File.ReadAllLines(fileName);
            string[] values;
            for (int i = 0; i < noOfVariables; i++)
            {
                for (int j = 0; j < noOfCoords; j++)
                {
                    values = lines[i + j * noOfVariables].Split(separator);
                    evp[j][i, 0] = double.Parse(values[0], ni);
                    evp[j][i, 1] = double.Parse(values[1], ni);
                    evp[j][i, 2] = double.Parse(values[2], ni);
                }

                values = lines[i + noOfCoords * noOfVariables].Split(separator);
                evpPos[i, 0] = int.Parse(values[0], ni);
                evpPos[i, 1] = int.Parse(values[1], ni);
                evpPos[i, 2] = int.Parse(values[2], ni);

                values = lines[i + (noOfCoords + 1) * noOfVariables].Split(separator);
                eigs[i] = double.Parse(values[0], ni);
            }
            values = lines[(noOfCoords + 2) * noOfVariables].Split(separator);
            mean = double.Parse(values[0], ni);
            values = lines[(noOfCoords + 2) * noOfVariables + 1].Split(separator);
            //covariance = double.Parse(values[0], ni);
            covariance = 0.15;
        }

        private double GetCoefficientInternal(double[] coordinates, int order)
        {
            if (order == -1) return 1;

            double f = covariance * Math.Sqrt(eigs[order]);
            for (int i = 0; i < coordinates.Length; i++)
            {
                int index = evpPos[order, i];
                double w = evp[i][index - 1, 0];
                double l = evp[i][index - 1, 1];
                double a = evp[i][index - 1, 2];
                double x = coordinates[i] - translations[i];

                f *= index % 2 == 0 ? a * Math.Sin(w * x) : a * Math.Cos(w * x);
            }

            return f;
        }

        #region IStochasticCoefficientsProvider Members

        double IStochasticCoefficientsProvider.GetCoefficient(double meanValue, double[] coordinates)
        {
            return GetCoefficientInternal(coordinates, CurrentOrder);
        }

        #endregion

        #region IStochasticMaterialCoefficientsProvider Members

        double IStochasticMaterialCoefficientsProvider.GetCoefficient(double meanValue, double[] coordinates)
        {
            double result = 1;
            for (int i = 0; i < RandomVariables.Length; i++)
                result += RandomVariables[i] * GetCoefficientInternal(coordinates, i);
            return result * meanValue;
        }
        
        #endregion
    }
}
