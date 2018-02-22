using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using ISAAR.MSolve.FEM.Interfaces;

namespace ISAAR.MSolve.FEM.Stochastic
{
    public class PolynomialChaosCoefficientsCalculator : IPolynomialChaosCoefficients
    {
        private int mOrder, pOrder;
        private readonly bool isGaussian;
        private readonly int psiSize;
        private int[][] psiBasis;
        private int[] psiSquareNorm;
        private double[][] hermitePolynomials;
        private List<Tuple<int, int, double>>[] gaussianCoefficients;
        private Dictionary<Tuple<int, int>, double>[] lognormalCoefficients;

        public int PsiSize { get { return psiSize; } }
        public int[][] PsiBasis { get { return psiBasis; } }
        public int[] PsiSquareNorm { get { return psiSquareNorm; } }
        public double[][] HermitePolynomials { get { return hermitePolynomials; } }
        public List<Tuple<int, int, double>>[] GaussianCoefficients { get { return gaussianCoefficients; } }
        public Dictionary<Tuple<int, int>, double>[] LognormalCoefficients { get { return lognormalCoefficients; } }
        private readonly Dictionary<int, double> factorials = new Dictionary<int, double>();
        private readonly Dictionary<int, double> doubleFactorials = new Dictionary<int, double>();

        public PolynomialChaosCoefficientsCalculator(int m, int p, bool isGaussian)
        {
            this.isGaussian = isGaussian;
            this.mOrder = m;
            this.pOrder = p;
            if (mOrder < 1) throw new ArgumentException("m cannot be less than 1");
            if (pOrder < 1) throw new ArgumentException("p cannot be less than 1");

            //psiSize = 1;
            //lognormalCoefficients = new Dictionary<Tuple<int, int>, double>[] { new Dictionary<Tuple<int, int>, double>() };
            //lognormalCoefficients[0].Add(new Tuple<int, int>(0, 0), 1.0);

            psiSize = 0;
            for (int k = 0; k <= pOrder; k++)
                psiSize += (int)(Factorial(mOrder + k - 1) / (Factorial(k) * Factorial(mOrder - 1)));

            CalculatePCCoefficients();
        }

        public void WriteToDisk(string filename)
        {
            List<string> lines = new List<string>();
            lines.Add(isGaussian ? "Gaussian" : "Lognormal");
            lines.Add(mOrder.ToString());
            lines.Add(pOrder.ToString());
            lines.Add(psiSize.ToString());
            for (int i = 0; i < psiSize; i++)
                for (int j = 0; j < mOrder; j++)
                    lines.Add(psiBasis[i][j].ToString());
            for (int i = 0; i < psiSize; i++)
                lines.Add(psiSquareNorm[i].ToString());
            for (int i = 0; i < pOrder; i++)
                for (int j = 0; j < i + 2; j++)
                    lines.Add(hermitePolynomials[i][j].ToString());
            if (isGaussian)
            {
                for (int i = 0; i < mOrder + 1; i++)
                {
                    lines.Add(gaussianCoefficients[i].Count.ToString());
                    for (int j = 0; j < gaussianCoefficients[i].Count; j++)
                        lines.Add(String.Format("{0};{1};{2}", gaussianCoefficients[i][j].Item1, gaussianCoefficients[i][j].Item2, gaussianCoefficients[i][j].Item3));
                }
            }
            else
            {
                for (int i = 0; i < psiSize; i++)
                {
                    lines.Add(lognormalCoefficients[i].Count.ToString());
                    foreach (var key in lognormalCoefficients[i].Keys)
                        lines.Add(String.Format("{0};{1};{2}", key.Item1, key.Item2, lognormalCoefficients[i][key]));
                }
            }
            File.WriteAllLines(filename, lines);
        }

        //private void ReadFromDisk(string filename)
        //{
        //    var lines = File.ReadAllLines(filename);
        //    isGaussian = lines[0] == "Gaussian";
        //    mOrder = Int32.Parse(lines[1]);
        //    pOrder = Int32.Parse(lines[2]);
        //    psiSize = Int32.Parse(lines[3]);
        //    psiBasis = new int[psiSize][];
        //    int index = 4;
        //    for (int i = 0; i < psiSize; i++)
        //    {
        //        psiBasis[i] = new int[mOrder];
        //        for (int j = 0; j < mOrder; j++)
        //        {
        //            psiBasis[i][j] = Int32.Parse(lines[index]);
        //            index++;
        //        }
        //    }
        //    psiSquareNorm = new int[psiSize];
        //    for (int i = 0; i < psiSize; i++)
        //    {
        //        psiSquareNorm[i] = Int32.Parse(lines[index]);
        //        index++;
        //    }
        //    hermitePolynomials = new double[pOrder][];
        //    for (int i = 0; i < pOrder; i++)
        //    {
        //        hermitePolynomials[i] = new double[i + 2];
        //        for (int j = 0; j < i + 2; j++)
        //        {
        //            hermitePolynomials[i][j] = Double.Parse(lines[index]);
        //            index++;
        //        }
        //    }
        //    if (isGaussian)
        //    {
        //        gaussianCoefficients = new List<Tuple<int, int, double>>[mOrder + 1];
        //        for (int i = 0; i < mOrder + 1; i++)
        //        {
        //            int count = Int32.Parse(lines[index]);
        //            index++;
        //            gaussianCoefficients[i] = new List<Tuple<int, int, double>>(count);
        //            for (int j = 0; j < count; j++)
        //            {
        //                var line = lines[index].Split(';');
        //                gaussianCoefficients[i].Add(new Tuple<int, int, double>(Int32.Parse(line[0]), Int32.Parse(line[1]), Double.Parse(line[2])));
        //                index++;
        //            }
        //        }
        //    }
        //    else
        //    {
        //        lognormalCoefficients = new Dictionary<Tuple<int, int>, double>[psiSize];
        //        for (int i = 0; i < psiSize; i++)
        //        {
        //            int count = Int32.Parse(lines[index]);
        //            index++;
        //            lognormalCoefficients[i] = new Dictionary<Tuple<int, int>, double>(count);
        //            for (int j = 0; j < count; j++)
        //            {
        //                var line = lines[index].Split(';');
        //                lognormalCoefficients[i].Add(new Tuple<int, int>(Int32.Parse(line[0]), Int32.Parse(line[1])), Double.Parse(line[2]));
        //                index++;
        //            }
        //        }
        //    }
        //}

        private double Factorial(int f)
        {
            if (factorials.ContainsKey(f))
                return factorials[f];

            double result = 1;
            for (int i = 2; i <= f; i++)
                result *= i;
            factorials.Add(f, result);
            return result;
        }

        private double DoubleFactorial(int f)
        {
            if (f % 2 == 1) return 0;
            if (doubleFactorials.ContainsKey(f))
                return doubleFactorials[f];

            int f2 = f / 2;
            int result = 1;
            for (int i = 2; i <= f2; i++)
                result *= 2 * i - 1;
            doubleFactorials[f] = result;

            return result;
        }

        private double[][] CalculateHermitePolynomials(int order)
        {
            double[][] p = new double[order][];
            p[0] = new double[2] { 1, 0 };
            for (int i = 1; i < order; i++)
            {
                int i1 = i + 1;
                p[i] = new double[i + 2];
                for (int j = 0; j <= i; j++)
                    p[i][j] = p[i - 1][j] * i1 / (i1 - j);
                p[i][i1] = i1 % 2 == 1 ? 0 : Math.Pow(-1, i1 / 2) * DoubleFactorial(i1);
            }

            return p;
        }

        private double[] Convolute(double[] v1, double[] v2)
        {
            int vecSize1 = v1.Length;
            int vecSize2 = v2.Length;
            int size = vecSize1 + vecSize2 - 1;
            var result = new double[size];
            var v1Ext = new double[size];
            var v2Ext = new double[size];

            var v2Reverse = v2.Reverse().ToArray();
            Array.Clear(v2Ext, 0, size);
            Array.Copy(v2Reverse, 0, v2Ext, size - vecSize2, vecSize2);
            int index = size - 1;
            for (int i = 0; i < size - vecSize1 + 1; i++)
            {
                Array.Clear(v1Ext, 0, size);
                Array.Copy(v1, 0, v1Ext, i, vecSize1);
                for (int j = 0; j < size; j++)
                    result[index] += v1Ext[j] * v2Ext[j];
                index--;
            }

            for (int i = size - vecSize2 - 1; i >= 0; i--)
            {
                Array.Clear(v2Ext, 0, size);
                Array.Copy(v2Reverse, 0, v2Ext, i, vecSize2);
                for (int j = 0; j < size; j++)
                    result[index] += v1Ext[j] * v2Ext[j];
                index--;
            }
            return result;
        }

        private void CalculatePCCoefficients()
        {
            psiBasis = new int[psiSize][];
            for (int i = 0; i < psiSize; i++)
                psiBasis[i] = new int[mOrder];
            int mm = mOrder - 1;
            int n = 0;

            if (mOrder == 1)
                for (int i = 0; i < pOrder; i++)
                    psiBasis[i + 1][0] = i;
            else
            {
                int[] t = new int[mm];
                for (int currentOrder = 1; currentOrder <= pOrder; currentOrder++)
                {
                    bool endOfGeneration = false;
                    bool firstThisOrder = false;
                    while (!endOfGeneration)
                    {
                        n++;
                        if (!firstThisOrder)
                        {
                            for (int i = 0; i < mm; i++)
                                t[i] = i + 1;
                            firstThisOrder = true;
                        }
                        else
                        {
                            if (t[mm - 1] < mm + currentOrder)
                                t[mm - 1] += 1;
                            else
                            {
                                int j = mm;
                                while (t[j - 1] == j + currentOrder)
                                    j--;
                                t[j - 1]++;
                                for (int k = j; k < mm; k++)
                                    t[k] = t[j - 1] + k - j + 1;
                            }
                        }

                        psiBasis[n][0] = t[0] - 1;
                        for (int i = 1; i < mm; i++)
                            psiBasis[n][i] = t[i] - t[i - 1] - 1;
                        psiBasis[n][mm] = mOrder + currentOrder - t[mm - 1] - 1;
                        if (t[0] == currentOrder + 1)
                            endOfGeneration = true;
                    }
                }
            }

            psiSquareNorm = new int[psiSize];
            int norme = 1;
            for (int k = 0; k < psiSize; k++)
            {
                norme = 1;
                for (int j = 0; j < mOrder; j++)
                    norme *= (int)Factorial(psiBasis[k][j]);
                psiSquareNorm[k] = norme;
            }

            hermitePolynomials = CalculateHermitePolynomials(pOrder);
            if (isGaussian)
            {
                gaussianCoefficients = new List<Tuple<int, int, double>>[mOrder + 1];
                gaussianCoefficients[0] = new List<Tuple<int, int, double>>();
                for (int i = 0; i < psiSize; i++)
                    gaussianCoefficients[0].Add(new Tuple<int, int, double>(i, i, psiSquareNorm[i]));
                int totalMatrices = 0;

                for (int i = 0; i < mOrder; i++)
                {
                    gaussianCoefficients[i + 1] = new List<Tuple<int, int, double>>();
                    for (int j = 0; j < psiSize; j++)
                        for (int k = 0; k < psiSize; k++)
                        {
                            bool almostEqual = false;
                            for (int l = 0; l < mOrder; l++)
                                if (l != i && psiBasis[j][l] != psiBasis[k][l])
                                {
                                    almostEqual = true;
                                    break;
                                }

                            if (!almostEqual)
                            {
                                if (psiBasis[k][i] == psiBasis[j][i] - 1)
                                {
                                    gaussianCoefficients[i + 1].Add(new Tuple<int, int, double>(j, k, psiSquareNorm[j]));
                                    totalMatrices++;
                                }
                                else if (psiBasis[k][i] == psiBasis[j][i] + 1)
                                {
                                    gaussianCoefficients[i + 1].Add(new Tuple<int, int, double>(j, k, psiSquareNorm[j] * (psiBasis[j][i] + 1)));
                                    totalMatrices++;
                                }
                            }
                        }
                }

                //var consolidatedCoefficients = new Dictionary<Tuple<int, int>, List<Tuple<double, int>>>();
                //for (int i = 0; i <= mOrder; i++)
                //{
                //    foreach (var x in gaussianCoefficients[i])
                //    {
                //        List<Tuple<double, int>> c = null;
                //        var key = new Tuple<int, int>(x.Item1, x.Item2);
                //        if (!consolidatedCoefficients.ContainsKey(key))
                //        {
                //            c = new List<Tuple<double, int>>();
                //            consolidatedCoefficients.Add(key, c);
                //        }
                //        else
                //            c = consolidatedCoefficients[key];

                //        c.Add(new Tuple<double, int>(x.Item3, i));
                //    }
                //}

                //var separateSums = new Dictionary<string, List<Tuple<int, int>>>();
                //foreach (var x in consolidatedCoefficients)
                //{
                //    var s = String.Format("{0}K{1}", x.Value[0].Item1, x.Value[0].Item2);
                //    for (int i = 1; i < x.Value.Count; i++)
                //        s += String.Format("+{0}K{1}", x.Value[i].Item1, x.Value[i].Item2);
                //    if (separateSums.ContainsKey(s))
                //        separateSums[s].Add(x.Key);
                //    else
                //    {
                //        var l = new List<Tuple<int, int>>();
                //        l.Add(x.Key);
                //        separateSums.Add(s, l);
                //    }
                //}

                //var multis = consolidatedCoefficients.Where(x => x.Value.Count > 1);
                //var multis2 = separateSums.Where(x => x.Value.Count > 2);
                //totalMatrices += totalMatrices + consolidatedCoefficients.Count + multis2.Count() + separateSums.Count;
            }
            else
            {
                lognormalCoefficients = new Dictionary<Tuple<int, int>, double>[psiSize];
                for (int i = 0; i < psiSize; i++)
                    lognormalCoefficients[i] = new Dictionary<Tuple<int, int>, double>();

                for (int i = 0; i < psiSize; i++)
                {
                    for (int j = 0; j < psiSize; j++)
                        for (int k = 0; k < psiSize; k++)
                        {
                            double esp = 1;
                            double[] productPol;
                            for (int l = 0; l < mOrder; l++)
                            {
                                if (psiBasis[i][l] == 0)
                                    if (psiBasis[j][l] == 0)
                                        if (psiBasis[k][l] == 0)
                                            productPol = new double[] { 1 };
                                        else
                                            productPol = hermitePolynomials[psiBasis[k][l] - 1];
                                    else
                                        if (psiBasis[k][l] == 0)
                                            productPol = hermitePolynomials[psiBasis[j][l] - 1];
                                        else
                                            productPol = Convolute(hermitePolynomials[psiBasis[j][l] - 1], hermitePolynomials[psiBasis[k][l] - 1]);
                                else
                                    if (psiBasis[j][l] == 0)
                                        if (psiBasis[k][l] == 0)
                                            productPol = hermitePolynomials[psiBasis[i][l] - 1];
                                        else
                                            productPol = Convolute(hermitePolynomials[psiBasis[i][l] - 1], hermitePolynomials[psiBasis[k][l] - 1]);
                                    else
                                        if (psiBasis[k][l] == 0)
                                            productPol = Convolute(hermitePolynomials[psiBasis[i][l] - 1], hermitePolynomials[psiBasis[j][l] - 1]);
                                        else
                                            productPol = Convolute(Convolute(hermitePolynomials[psiBasis[i][l] - 1], hermitePolynomials[psiBasis[j][l] - 1]), hermitePolynomials[psiBasis[k][l] - 1]);

                                int degPol = productPol.Length - 1;
                                double esperanceProd = productPol[degPol];
                                for (int q = 0; q < degPol - 1; q++)
                                    esperanceProd += DoubleFactorial(degPol - q) * productPol[q];
                                esp *= esperanceProd;
                                if (esp == 0)
                                    break;
                            }

                            if (esp != 0)
                            {
                                var jk = new Tuple<int, int>(j, k);
                                var kj = new Tuple<int, int>(k, j);
                                var ik = new Tuple<int, int>(i, k);
                                var ki = new Tuple<int, int>(k, i);
                                var ij = new Tuple<int, int>(i, j);
                                var ji = new Tuple<int, int>(j, i);
                                if (lognormalCoefficients[i].ContainsKey(jk))
                                    lognormalCoefficients[i][jk] = esp;
                                else
                                    lognormalCoefficients[i].Add(jk, esp);
                                if (lognormalCoefficients[i].ContainsKey(kj))
                                    lognormalCoefficients[i][kj] = esp;
                                else
                                    lognormalCoefficients[i].Add(kj, esp);

                                if (lognormalCoefficients[j].ContainsKey(ik))
                                    lognormalCoefficients[j][ik] = esp;
                                else
                                    lognormalCoefficients[j].Add(ik, esp);
                                if (lognormalCoefficients[j].ContainsKey(ki))
                                    lognormalCoefficients[j][ki] = esp;
                                else
                                    lognormalCoefficients[j].Add(ki, esp);

                                if (lognormalCoefficients[k].ContainsKey(ij))
                                    lognormalCoefficients[k][ij] = esp;
                                else
                                    lognormalCoefficients[k].Add(ij, esp);
                                if (lognormalCoefficients[k].ContainsKey(ji))
                                    lognormalCoefficients[k][ji] = esp;
                                else
                                    lognormalCoefficients[k].Add(ji, esp);
                            }
                        }
                }
            }
        }
    }
}
