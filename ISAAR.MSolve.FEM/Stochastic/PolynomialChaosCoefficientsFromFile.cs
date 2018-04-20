using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using ISAAR.MSolve.FEM.Interfaces;

namespace ISAAR.MSolve.FEM.Stochastic
{
    public class PolynomialChaosCoefficientsFromFile : IPolynomialChaosCoefficients
    {
        private readonly int mOrder, pOrder;
        private readonly bool isGaussian;
        private readonly int psiSize;
        private readonly int[][] psiBasis;
        private readonly int[] psiSquareNorm;
        private readonly double[][] hermitePolynomials;
        private readonly List<Tuple<int, int, double>>[] gaussianCoefficients;
        private readonly Dictionary<Tuple<int, int>, double>[] lognormalCoefficients;

        public int PsiSize { get { return psiSize; } }
        public int[][] PsiBasis { get { return psiBasis; } }
        public int[] PsiSquareNorm { get { return psiSquareNorm; } }
        public double[][] HermitePolynomials { get { return hermitePolynomials; } }
        public List<Tuple<int, int, double>>[] GaussianCoefficients { get { return gaussianCoefficients; } }
        public Dictionary<Tuple<int, int>, double>[] LognormalCoefficients { get { return lognormalCoefficients; } }

        public PolynomialChaosCoefficientsFromFile(string filename)
        {
            var lines = File.ReadAllLines(filename);
            isGaussian = lines[0] == "Gaussian";
            mOrder = Int32.Parse(lines[1]);
            pOrder = Int32.Parse(lines[2]);
            psiSize = Int32.Parse(lines[3]);
            psiBasis = new int[psiSize][];
//            int index = 4;
            int index = 5;
            for (int i = 0; i < psiSize; i++)
            {
                psiBasis[i] = new int[mOrder];
                for (int j = 0; j < mOrder; j++)
                {
                    psiBasis[i][j] = Int32.Parse(lines[index]);
                    index++;
                }
            }
            psiSquareNorm = new int[psiSize];
            for (int i = 0; i < psiSize; i++)
            {
                psiSquareNorm[i] = Int32.Parse(lines[index]);
                index++;
            }
            hermitePolynomials = new double[pOrder][];
            for (int i = 0; i < pOrder; i++)
            {
                hermitePolynomials[i] = new double[i + 2];
                for (int j = 0; j < i + 2; j++)
                {
                    hermitePolynomials[i][j] = Double.Parse(lines[index]);
                    index++;
                }
            }
            if (isGaussian)
            {
                gaussianCoefficients = new List<Tuple<int, int, double>>[mOrder + 1];
                int totalMatrices = 0;
                for (int i = 0; i < mOrder + 1; i++)
                {
                    int count = Int32.Parse(lines[index]);
                    totalMatrices += count;
                    index++;
                    index++;
                    gaussianCoefficients[i] = new List<Tuple<int, int, double>>(count);
                    for (int j = 0; j < count; j++)
                    {
                        var line = lines[index].Split(';');
                        gaussianCoefficients[i].Add(new Tuple<int, int, double>(Int32.Parse(line[0]) - 1, Int32.Parse(line[1]) - 1, Double.Parse(line[2])));
                        index++;
                    }
                }
                //psiSize += totalMatrices;
            }
            else
            {
                lognormalCoefficients = new Dictionary<Tuple<int, int>, double>[psiSize];
                int totalMatrices = 0;
                for (int i = 0; i < psiSize; i++)
                {
                    int count = Int32.Parse(lines[index]);
                    totalMatrices += count;
                    index++;
                    index++;
                    lognormalCoefficients[i] = new Dictionary<Tuple<int, int>, double>(count);
                    for (int j = 0; j < count; j++)
                    {
                        var line = lines[index].Split(';');
                        lognormalCoefficients[i].Add(new Tuple<int, int>(Int32.Parse(line[0]) - 1, Int32.Parse(line[1]) - 1), Double.Parse(line[2]));
                        index++;
                    }
                }

                //var consolidatedCoefficients = new Dictionary<Tuple<int, int>, List<Tuple<double, int>>>();
                //for (int i = 0; i < psiSize; i++)
                //{
                //    foreach (var x in lognormalCoefficients[i])
                //    {
                //        List<Tuple<double, int>> c = null;
                //        if (!consolidatedCoefficients.ContainsKey(x.Key))
                //        {
                //            c = new List<Tuple<double, int>>();
                //            consolidatedCoefficients.Add(x.Key, c);
                //        }
                //        else
                //            c = consolidatedCoefficients[x.Key];

                //        c.Add(new Tuple<double, int>(x.Value, i));
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
                //psiSize += totalMatrices + consolidatedCoefficients.Count + multis2.Count() + separateSums.Count;
            }
        }
    }
}
