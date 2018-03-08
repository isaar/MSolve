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
    public class StochasticCoefficientsMemoizer
    {
        private readonly Dictionary<Tuple<int, double, double, double>, double> coefficientsDictionary = new Dictionary<Tuple<int, double, double, double>, double>();

        public bool Exists(double[] coordinates, int order)
        {
            double x = coordinates.Length > 0 ? coordinates[0] : 0.0;
            double y = coordinates.Length > 1 ? coordinates[1] : 0.0;
            double z = coordinates.Length > 2 ? coordinates[2] : 0.0;

            return coefficientsDictionary.ContainsKey(new Tuple<int, double, double, double>(order, x, y, z));
        }

        public double GetCoefficient(double[] coordinates, int order)
        {
            double x = coordinates.Length > 0 ? coordinates[0] : 0.0;
            double y = coordinates.Length > 1 ? coordinates[1] : 0.0;
            double z = coordinates.Length > 2 ? coordinates[2] : 0.0;

            return coefficientsDictionary[new Tuple<int, double, double, double>(order, x, y, z)];
        }

        public void SetCoefficient(double[] coordinates, int order, double coefficient)
        {
            double x = coordinates.Length > 0 ? coordinates[0] : 0.0;
            double y = coordinates.Length > 1 ? coordinates[1] : 0.0;
            double z = coordinates.Length > 2 ? coordinates[2] : 0.0;

            coefficientsDictionary.Add(new Tuple<int, double, double, double>(order, x, y, z), coefficient);
        }
    }

    public class LognormalFileStochasticCoefficientsProvider : IStochasticMaterialCoefficientsProvider, IStochasticCoefficientsProvider
    {
        private readonly StochasticCoefficientsMemoizer memoizer;
        private readonly double[][,] evp;
        private readonly int[,] evpPos;
        private readonly double[] eigs;
        private readonly double[] translations;
        private readonly double mean;
        private readonly double meanGaussian;
        private readonly double covariance;
        private readonly double covarianceGaussian;
        private readonly double[] stochasticDomain;
        private readonly int[] varianceAxes = new int[] { 0, 1, 2 };
        private static NumberFormatInfo ni = null;
        private double[] randomVariables;
        private bool changedVariables = false;

        public int CurrentOrder { get; set; }
        public double[] RandomVariables 
        {
            get { return randomVariables; }
            set 
            {
                changedVariables = true;
                randomVariables = value; 
            }
        }

        public LognormalFileStochasticCoefficientsProvider(string fileName, int noOfVariables, char separator, double[] stochasticDomain, string decimalSeparator = ".")
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
            //mean = double.Parse(values[0], ni);
            values = lines[(noOfCoords + 2) * noOfVariables + 1].Split(separator);
            //covariance = double.Parse(values[0], ni);

            mean = 1;
            covariance = 0.5;
            covarianceGaussian = Math.Sqrt(Math.Log(Math.Pow(covariance / mean, 2) + 1, Math.E));
            meanGaussian = 0;
            //meanGaussian = Math.Log(mean, Math.E) - 0.5 * Math.Pow(covarianceGaussian, 2);
        }

        public LognormalFileStochasticCoefficientsProvider(string fileName, int noOfVariables, char separator, double[] stochasticDomain, int[] varianceAxes, string decimalSeparator = ".")
            : this(fileName, noOfVariables, separator, stochasticDomain, decimalSeparator)
        {
            this.varianceAxes = varianceAxes;
        }

        public LognormalFileStochasticCoefficientsProvider(string fileName, int noOfVariables, char separator, double[] stochasticDomain, StochasticCoefficientsMemoizer memoizer, string decimalSeparator = ".")
            : this(fileName, noOfVariables, separator, stochasticDomain, decimalSeparator)
        {
            this.memoizer = memoizer;
        }

        private double GetCoefficientInternal(double[] coordinates, int order)
        {
            if (order == -1) return 1;

            double f = 1;
            if (memoizer == null || memoizer.Exists(coordinates, order) == false)
            {
                f = covarianceGaussian * Math.Sqrt(eigs[order]);
                for (int i = 0; i < varianceAxes.Length; i++)
                {
                    int index = evpPos[order, varianceAxes[i]];
                    double w = evp[i][index - 1, 0];
                    double l = evp[i][index - 1, 1];
                    double a = evp[i][index - 1, 2];
                    double x = coordinates[varianceAxes[i]] - translations[varianceAxes[i]];

                    f *= index % 2 == 0 ? a * Math.Sin(w * x) : a * Math.Cos(w * x);
                }
                if (memoizer != null)
                    memoizer.SetCoefficient(coordinates, order, f);
            }
            else
                f = memoizer.GetCoefficient(coordinates, order);
            f += meanGaussian;
            
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
            double result = Math.Log(meanValue, Math.E) - 0.5 * Math.Pow(covarianceGaussian, 2);
            for (int i = 0; i < RandomVariables.Length; i++)
                result += RandomVariables[i] * GetCoefficientInternal(coordinates, i);

            //if (changedVariables)
            //{
            //    changedVariables = false;
            //    double test = 0;
            //    for (int i = 0; i < RandomVariables.Length; i++)
            //        test += RandomVariables[i] * GetCoefficientInternal(new double[] { 0, 0, 0 }, i);
            //    using (var sw = File.AppendText(@"d:\XXX.TXT"))
            //    {
            //        sw.WriteLine(test.ToString());
            //    }
            //}

            return Math.Exp(result);
            //return Math.Exp(result) - mean + 1.0;
        }
        
        #endregion
    }
}
