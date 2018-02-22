using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.FEM.Stochastic;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.PCGSkyline;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Analyzers
{
    public enum PolynomialChaosSolverPCGPreconditioner
    {
        BlockDiagonal,
        BlockSSOR,
        ILU,
        SSOR,
        Diagonal
    }
    
    public class PolynomialChaosSolverPCGMatrixCalculator : ISolverPCGMatrixCalculator
    {
        private SolverPCG<SkylineMatrix2D> solver;
        private readonly int mOrder, pOrder;
        private readonly PolynomialChaosSolverPCGPreconditioner preconditionerKind;
        //private readonly int psiSize;
        //private int[][] psiBasis;
        //private int[] psiSquareNorm;
        //private List<Tuple<int, int, double>>[] gaussianCoefficients;
        //private Dictionary<Tuple<int, int>, double>[] lognormalCoefficients;
        private readonly bool isGaussian;
        private readonly PolynomialChaosCoefficientsCalculator calculator;
        private double[] meanDiagonalCoefficients, ssorDiagonals;
        private PolynomialChaosAnalyzer analyzer;
        private int meanSize;
        private Vector tempVectorIn, tempVectorOut;
        private List<List<Tuple<int, int, int, double>>> coeffLower, coeffUpper;

        public SolverPCG<SkylineMatrix2D> Solver
        {
            get { return solver; }
            set { solver = value; }
        }

        public PolynomialChaosAnalyzer Analyzer
        {
            get { return analyzer; }
            set { analyzer = value; }
        }

        public PolynomialChaosSolverPCGMatrixCalculator(int m, int p, bool isGaussian, PolynomialChaosSolverPCGPreconditioner preconditionerKind)
        {
            this.preconditionerKind = preconditionerKind;
            this.isGaussian = isGaussian;
            this.calculator = new PolynomialChaosCoefficientsCalculator(m, p, isGaussian);
            this.mOrder = m;
            this.pOrder = p;
            //if (mOrder < 1) throw new ArgumentException("m cannot be less than 1");
            //if (pOrder < 1) throw new ArgumentException("p cannot be less than 1");

            //psiSize = 0;
            //for (int k = 0; k <= pOrder; k++)
            //    psiSize += Factorial(mOrder + k - 1) / (Factorial(k) * Factorial(mOrder - 1));

            //CalculatePCCoefficients();
            CalculateMeanDiagonalCoefficients();
            if (isGaussian)
                MakeGaussianCache();
            else
                MakeLognormalCache();
        }

        private void MakeSSORDiagonals()
        {
            if (ssorDiagonals != null) return;

            ssorDiagonals = new double[this.VectorSize];
            //foreach (ILinearSystem subdomain in solver.SubdomainsDictionary.Values)
            //foreach (ILinearSystem subdomain in solver.SubdomainsDictionary.Values)
            //{
                for (int i = 0; i < coeffLower.Count; i++)
                {
                    int blockStart = i * meanSize;
                    foreach (var c in coeffLower[i].Where(x => x.Item2 == x.Item3))
                    {
                        SkylineMatrix2D K = (SkylineMatrix2D)analyzer.Matrices[c.Item1][solver.LinearSystem.ID];
                        double[] d = K.Data as double[];
                        for (int n = 0; n < K.Rows; n++)
                            ssorDiagonals[blockStart + n] += d[K.RowIndex[n]] * c.Item4;
                    }
                //}
            }
            //for (int n = 0; n < ssorDiagonals.Length; n++)
            //    ssorDiagonals[n] = 1.0 / ssorDiagonals[n];
        }

        public int VectorSize
        {
            get
            {
                //meanSize = solver.SubdomainsDictionary.Values.First().RHS.Length;
                meanSize = solver.LinearSystem.RHS.Length;
                return meanSize * calculator.PsiSize;
            }
        }

        private void MakeGaussianCache()
        {
            var cLower = calculator.GaussianCoefficients.SelectMany((x, index) =>
            {
                var a = new List<Tuple<int, int, int, double>>();
                foreach (var item in x)
                    if (item.Item1 >= item.Item2)
                        a.Add(new Tuple<int, int, int, double>(index, item.Item1, item.Item2, item.Item3));

                return a;
            })
            .OrderBy(x => x.Item1).OrderBy(x => x.Item3).OrderBy(x => x.Item2)
            .GroupBy(x => x.Item2);

            coeffLower = new List<List<Tuple<int, int, int, double>>>();
            foreach (var elements in cLower)
            {
                var currentList = new List<Tuple<int, int, int, double>>();
                foreach (var element in elements)
                    currentList.Add(element);
                coeffLower.Add(currentList);
            }

            var cUpper = calculator.GaussianCoefficients.SelectMany((x, index) =>
            {
                var a = new List<Tuple<int, int, int, double>>();
                foreach (var item in x)
                    if (item.Item1 <= item.Item2)
                        a.Add(new Tuple<int, int, int, double>(index, item.Item1, item.Item2, item.Item3));

                return a;
            })
            .OrderBy(x => x.Item1).OrderBy(x => x.Item3).OrderBy(x => x.Item2)
            .GroupBy(x => x.Item2);

            coeffUpper = new List<List<Tuple<int, int, int, double>>>();
            foreach (var elements in cUpper)
            {
                var currentList = new List<Tuple<int, int, int, double>>();
                foreach (var element in elements)
                    currentList.Add(element);
                coeffUpper.Add(currentList);
            }
        }

        private void MakeLognormalCache()
        {
            var cLower = calculator.LognormalCoefficients.SelectMany((x, index) =>
            {
                var a = new List<Tuple<int, int, int, double>>();
                foreach (var item in x)
                    if (item.Key.Item1 >= item.Key.Item2)
                        a.Add(new Tuple<int, int, int, double>(index, item.Key.Item1, item.Key.Item2, item.Value));

                return a;
            })
            .OrderBy(x => x.Item1).OrderBy(x => x.Item3).OrderBy(x => x.Item2)
            .GroupBy(x => x.Item2);

            coeffLower = new List<List<Tuple<int, int, int, double>>>();
            foreach (var elements in cLower)
            {
                var currentList = new List<Tuple<int, int, int, double>>();
                foreach (var element in elements)
                    currentList.Add(element);
                coeffLower.Add(currentList);
            }

            var cUpper = calculator.LognormalCoefficients.SelectMany((x, index) =>
            {
                var a = new List<Tuple<int, int, int, double>>();
                foreach (var item in x)
                    if (item.Key.Item1 <= item.Key.Item2)
                        a.Add(new Tuple<int, int, int, double>(index, item.Key.Item1, item.Key.Item2, item.Value));

                return a;
            })
            .OrderBy(x => x.Item1).OrderBy(x => x.Item3).OrderBy(x => x.Item2)
            .GroupBy(x => x.Item2);

            coeffUpper = new List<List<Tuple<int, int, int, double>>>();
            foreach (var elements in cUpper)
            {
                var currentList = new List<Tuple<int, int, int, double>>();
                foreach (var element in elements)
                    currentList.Add(element);
                coeffUpper.Add(currentList);
            }
        }

        //private int Factorial(int f)
        //{
        //    int result = 1;
        //    for (int i = 2; i <= f; i++)
        //        result *= i;
        //    return result;
        //}

        //private int DoubleFactorial(int f)
        //{
        //    if (f % 2 == 1) return 0;

        //    int f2 = f / 2;
        //    int result = 1;
        //    for (int i = 2; i <= f2; i++)
        //        result *= 2 * i - 1;
        //    return result;
        //}

        //private double[][] HermitePolynomials(int order)
        //{
        //    double[][] p = new double[order][];
        //    p[0] = new double[2] { 1, 0 };
        //    for (int i = 1; i < order; i++)
        //    {
        //        int i1 = i + 1;
        //        p[i] = new double[i + 2];
        //        for (int j = 0; j <= i; j++)
        //            p[i][j] = p[i - 1][j] * i1 / (i1 - j);
        //        p[i][i1] = i1 % 2 == 1 ? 0 : Math.Pow(-1, i1 / 2) * DoubleFactorial(i1);
        //    }

        //    return p;
        //}

        //private double[] Convolute(double[] v1, double[] v2)
        //{
        //    int vecSize1 = v1.Length;
        //    int vecSize2 = v2.Length;
        //    int size = vecSize1 + vecSize2 - 1;
        //    var result = new double[size];
        //    var v1Ext = new double[size];
        //    var v2Ext = new double[size];

        //    var v2Reverse = v2.Reverse().ToArray();
        //    Array.Clear(v2Ext, 0, size);
        //    Array.Copy(v2Reverse, 0, v2Ext, size - vecSize2, vecSize2);
        //    int index = size - 1;
        //    for (int i = 0; i < size - vecSize1 + 1; i++)
        //    {
        //        Array.Clear(v1Ext, 0, size);
        //        Array.Copy(v1, 0, v1Ext, i, vecSize1);
        //        for (int j = 0; j < size; j++)
        //            result[index] += v1Ext[j] * v2Ext[j];
        //        index--;
        //    }

        //    for (int i = size - vecSize2 - 1; i >= 0; i--)
        //    {
        //        Array.Clear(v2Ext, 0, size);
        //        Array.Copy(v2Reverse, 0, v2Ext, i, vecSize2);
        //        for (int j = 0; j < size; j++)
        //            result[index] += v1Ext[j] * v2Ext[j];
        //        index--;
        //    }
        //    return result;
        //}

        //private void CalculatePCCoefficients()
        //{
        //    psiBasis = new int[psiSize][];
        //    for (int i = 0; i < psiSize; i++)
        //        psiBasis[i] = new int[mOrder];
        //    int mm = mOrder - 1;
        //    int n = 0;

        //    if (mOrder == 1)
        //        for (int i = 0; i < pOrder; i++)
        //            psiBasis[i + 1][0] = i;
        //    else
        //    {
        //        int[] t = new int[mm];
        //        for (int currentOrder = 1; currentOrder <= pOrder; currentOrder++)
        //        {
        //            bool endOfGeneration = false;
        //            bool firstThisOrder = false;
        //            while (!endOfGeneration)
        //            {
        //                n++;
        //                if (!firstThisOrder)
        //                {
        //                    for (int i = 0; i < mm; i++)
        //                        t[i] = i + 1;
        //                    firstThisOrder = true;
        //                }
        //                else
        //                {
        //                    if (t[mm - 1] < mm + currentOrder)
        //                        t[mm - 1] += 1;
        //                    else
        //                    {
        //                        int j = mm;
        //                        while (t[j - 1] == j + currentOrder)
        //                            j--;
        //                        t[j - 1]++;
        //                        for (int k = j; k < mm; k++)
        //                            t[k] = t[j - 1] + k - j + 1;
        //                    }
        //                }

        //                psiBasis[n][0] = t[0] - 1;
        //                for (int i = 1; i < mm; i++)
        //                    psiBasis[n][i] = t[i] - t[i - 1] - 1;
        //                psiBasis[n][mm] = mOrder + currentOrder - t[mm - 1] - 1;
        //                if (t[0] == currentOrder + 1)
        //                    endOfGeneration = true;
        //            }
        //        }
        //    }

        //    psiSquareNorm = new int[psiSize];
        //    int norme = 1;
        //    for (int k = 0; k < psiSize; k++)
        //    {
        //        norme = 1;
        //        for (int j = 0; j < mOrder; j++)
        //            norme *= Factorial(psiBasis[k][j]);
        //        psiSquareNorm[k] = norme;
        //    }

        //    gaussianCoefficients = new List<Tuple<int, int, double>>[mOrder];
        //    for (int i = 0; i < mOrder; i++)
        //    {
        //        gaussianCoefficients[i] = new List<Tuple<int, int, double>>();
        //        for (int j = 0; j < psiSize; j++)
        //            for (int k = 0; k < psiSize; k++)
        //            {
        //                bool almostEqual = false;
        //                for (int l = 0; l < mOrder; l++)
        //                    if (l != i && psiBasis[j][l] != psiBasis[k][l])
        //                    {
        //                        almostEqual = true;
        //                        break;
        //                    }

        //                if (!almostEqual)
        //                    if (psiBasis[k][i] == psiBasis[j][i] - 1)
        //                        gaussianCoefficients[i].Add(new Tuple<int, int, double>(j, k, psiSquareNorm[j]));
        //                    else if (psiBasis[k][i] == psiBasis[j][i] + 1)
        //                        gaussianCoefficients[i].Add(new Tuple<int, int, double>(j, k, psiSquareNorm[j] * (psiBasis[j][i] + 1)));
        //            }
        //    }

        //    var hp = HermitePolynomials(pOrder);
        //    lognormalCoefficients = new Dictionary<Tuple<int, int>, double>[psiSize];
        //    for (int i = 0; i < psiSize; i++)
        //        lognormalCoefficients[i] = new Dictionary<Tuple<int, int>, double>();

        //    for (int i = 0; i < psiSize; i++)
        //    {
        //        for (int j = 0; j < psiSize; j++)
        //            for (int k = 0; k < psiSize; k++)
        //            {
        //                double esp = 1;
        //                double[] productPol;
        //                for (int l = 0; l < mOrder; l++)
        //                {
        //                    if (psiBasis[i][l] == 0)
        //                        if (psiBasis[j][l] == 0)
        //                            if (psiBasis[k][l] == 0)
        //                                productPol = new double[] { 1 };
        //                            else
        //                                productPol = hp[psiBasis[k][l] - 1];
        //                        else
        //                            if (psiBasis[k][l] == 0)
        //                                productPol = hp[psiBasis[j][l] - 1];
        //                            else
        //                                productPol = Convolute(hp[psiBasis[j][l] - 1], hp[psiBasis[k][l] - 1]);
        //                    else
        //                        if (psiBasis[j][l] == 0)
        //                            if (psiBasis[k][l] == 0)
        //                                productPol = hp[psiBasis[i][l] - 1];
        //                            else
        //                                productPol = Convolute(hp[psiBasis[i][l] - 1], hp[psiBasis[k][l] - 1]);
        //                        else
        //                            if (psiBasis[k][l] == 0)
        //                                productPol = Convolute(hp[psiBasis[i][l] - 1], hp[psiBasis[j][l] - 1]);
        //                            else
        //                                productPol = Convolute(Convolute(hp[psiBasis[i][l] - 1], hp[psiBasis[j][l] - 1]), hp[psiBasis[k][l] - 1]);

        //                    int degPol = productPol.Length - 1;
        //                    double esperanceProd = productPol[degPol];
        //                    for (int q = 0; q < degPol - 1; q++)
        //                        esperanceProd += DoubleFactorial(degPol - q) * productPol[q];
        //                    esp *= esperanceProd;
        //                    if (esp == 0)
        //                        break;
        //                }

        //                if (esp != 0)
        //                {
        //                    var jk = new Tuple<int, int>(j, k);
        //                    var kj = new Tuple<int, int>(k, j);
        //                    var ik = new Tuple<int, int>(i, k);
        //                    var ki = new Tuple<int, int>(k, i);
        //                    var ij = new Tuple<int, int>(i, j);
        //                    var ji = new Tuple<int, int>(j, i);
        //                    if (lognormalCoefficients[i].ContainsKey(jk))
        //                        lognormalCoefficients[i][jk] = esp;
        //                    else
        //                        lognormalCoefficients[i].Add(jk, esp);
        //                    if (lognormalCoefficients[i].ContainsKey(kj))
        //                        lognormalCoefficients[i][kj] = esp;
        //                    else
        //                        lognormalCoefficients[i].Add(kj, esp);

        //                    if (lognormalCoefficients[j].ContainsKey(ik))
        //                        lognormalCoefficients[j][ik] = esp;
        //                    else
        //                        lognormalCoefficients[j].Add(ik, esp);
        //                    if (lognormalCoefficients[j].ContainsKey(ki))
        //                        lognormalCoefficients[j][ki] = esp;
        //                    else
        //                        lognormalCoefficients[j].Add(ki, esp);

        //                    if (lognormalCoefficients[k].ContainsKey(ij))
        //                        lognormalCoefficients[k][ij] = esp;
        //                    else
        //                        lognormalCoefficients[k].Add(ij, esp);
        //                    if (lognormalCoefficients[k].ContainsKey(ji))
        //                        lognormalCoefficients[k][ji] = esp;
        //                    else
        //                        lognormalCoefficients[k].Add(ji, esp);
        //                }
        //            }
        //    }
        //}

        private void CalculateMeanDiagonalCoefficients()
        {
            meanDiagonalCoefficients = new double[calculator.PsiSize];
            if (isGaussian)
                foreach (var c in calculator.GaussianCoefficients[0].Where(x => x.Item1 == x.Item2))
                    meanDiagonalCoefficients[c.Item1] = c.Item3;
            else
                foreach (var c in calculator.LognormalCoefficients[0].Where(x => x.Key.Item1 == x.Key.Item2))
                    meanDiagonalCoefficients[c.Key.Item1] = c.Value;
        }

        public void Precondition(IVector vIn, IVector vOut)
        {
            switch (preconditionerKind)
            {
                case PolynomialChaosSolverPCGPreconditioner.BlockDiagonal:
                    PreconditionBlock(vIn, vOut);
                    break;
                case PolynomialChaosSolverPCGPreconditioner.BlockSSOR:
                    PreconditionBlockSSOR(vIn, vOut);
                    break;
                case PolynomialChaosSolverPCGPreconditioner.ILU:
                    PreconditionILU(vIn, vOut);
                    break;
                case PolynomialChaosSolverPCGPreconditioner.SSOR:
                    PreconditionSSOR(vIn, vOut);
                    break;
                case PolynomialChaosSolverPCGPreconditioner.Diagonal:
                    PreconditionDiagonal(vIn, vOut);
                    break;
            }
        }

        private void PreconditionDiagonal(IVector vIn, IVector vOut)
        {
            if (analyzer.FactorizedMatrices.Count < 1)
                throw new InvalidOperationException("Cannot precondition with no factorized matrices");
            if (analyzer.FactorizedMatrices[0].Count != 1)
                throw new InvalidOperationException("Cannot precondition with more than one subdomains");

            if (tempVectorIn == null)
                tempVectorIn = new Vector(meanSize);
            if (tempVectorOut == null)
                tempVectorOut = new Vector(meanSize);

            var v = (Vector)vIn;
            var o = (Vector)vOut;
            //Array.Copy(v.Data, 0, o.Data, 0, v.Data.Length);

            //foreach (ILinearSystem subdomain in solver.SubdomainsDictionary.Values)
            //{
                SkylineMatrix2D k = (SkylineMatrix2D)analyzer.Matrices[0][solver.LinearSystem.ID];
                for (int i = 0; i < calculator.PsiSize; i++)
                {
                    Array.Copy(v.Data, i * meanSize, tempVectorIn.Data, 0, meanSize);
                    tempVectorIn.Multiply(1.0 / meanDiagonalCoefficients[i]);
                    for (int j = 0; j < k.Rows; j++)
                        tempVectorOut[j] = tempVectorIn[j] / k.Data[k.RowIndex[j]];
                    Array.Copy(tempVectorOut.Data, 0, o.Data, i * meanSize, meanSize);
                }
            //}
        }

        private void PreconditionBlock(IVector vIn, IVector vOut)
        {
            if (analyzer.FactorizedMatrices.Count < 1)
                throw new InvalidOperationException("Cannot precondition with no factorized matrices");
            if (analyzer.FactorizedMatrices[0].Count != 1)
                throw new InvalidOperationException("Cannot precondition with more than one subdomains");

            if (tempVectorIn == null)
                tempVectorIn = new Vector(meanSize);
            if (tempVectorOut == null)
                tempVectorOut = new Vector(meanSize);

            var v = vIn;
            var o = vOut;
            //var v = (Vector<double>)vIn;
            //var o = (Vector<double>)vOut;
            //Array.Copy(v.Data, 0, o.Data, 0, v.Data.Length);

            foreach (var m in analyzer.FactorizedMatrices[0].Values)
                for (int i = 0; i < calculator.PsiSize; i++)
                {
                    tempVectorIn.CopyFrom(0, meanSize, v, i * meanSize);
                    //Array.Copy(v.Data, i * meanSize, tempVectorIn.Data, 0, meanSize);
                    tempVectorIn.Multiply(1.0 / meanDiagonalCoefficients[i]);
                    m.Solve(tempVectorIn, tempVectorOut);
                    o.CopyFrom(i * meanSize, meanSize, tempVectorOut, 0);
                    //Array.Copy(tempVectorOut.Data, 0, o.Data, i * meanSize, meanSize);
                }
        }

        private void PreconditionBlockSSOR(IVector vIn, IVector vOut)
        {
            if (analyzer.FactorizedMatrices.Count < 1)
                throw new InvalidOperationException("Cannot precondition with no factorized matrices");
            if (analyzer.FactorizedMatrices[0].Count != 1)
                throw new InvalidOperationException("Cannot precondition with more than one subdomains");

            //MakeSSORDiagonals();
            //for (int i = 0; i < this.VectorSize; i++)
            //{
            //    vOut[i] = vIn[i] * ssorDiagonals[i];
            //}
            //return;

            if (tempVectorIn == null)
                tempVectorIn = new Vector(meanSize);
            if (tempVectorOut == null)
                tempVectorOut = new Vector(meanSize);

            //var v = (Vector<double>)vIn;
            //var o = (Vector<double>)vOut;
            var v = vIn;
            var o = vOut;
            var oTemp = new Vector(this.VectorSize);

            //foreach (ILinearSystem subdomain in solver.SubdomainsDictionary.Values)
            //{
                // Forward sub
                o.CopyFrom(0, v.Length, v, 0);
                //Array.Copy(v.Data, 0, o.Data, 0, v.Data.Length);
                for (int i = 0; i < coeffLower.Count; i++)
                {
                    int blockStart = i * meanSize;

                    tempVectorOut.CopyFrom(0, tempVectorOut.Length, v, blockStart);
                    //Array.Copy(v.Data, blockStart, tempVectorOut.Data, 0, meanSize);

                    //Array.Clear(tempVectorIn.Data, 0, meanSize);
                    //Array.Clear(tempVectorOut.Data, 0, meanSize);

                    foreach (var c in coeffLower[i].Where(x => x.Item2 != x.Item3))
                    {
                        SkylineMatrix2D k = (SkylineMatrix2D)analyzer.Matrices[c.Item1][solver.LinearSystem.ID];
                        int row = c.Item2 * k.Rows;
                        int col = c.Item3 * k.Rows;
                        k.Multiply(o, tempVectorOut.Data, -c.Item4, col, 0, false);
                    }

                    Array.Clear(tempVectorIn.Data, 0, meanSize);
                    foreach (var c in coeffLower[i].Where(x => x.Item2 == x.Item3 && x.Item1 == 0))
                    {
                        //Array.Copy(v.Data, i * meanSize, tempVectorIn.Data, 0, meanSize);
                        Array.Copy(tempVectorOut.Data, 0, tempVectorIn.Data, 0, meanSize);
                        tempVectorIn.Multiply(1.0 / c.Item4);
                        analyzer.FactorizedMatrices[0][1].Solve(tempVectorIn, tempVectorOut);
                        o.CopyFrom(0, meanSize, tempVectorOut, i * meanSize);
                        //Array.Copy(tempVectorOut.Data, 0, o.Data, i * meanSize, meanSize);

                    //    SkylineMatrix2D<double> K = (SkylineMatrix2D<double>)analyzer.Matrices[c.Item1][subdomain.ID];
                    //    double[] d = K.Data as double[];
                    //    for (int n = 0; n < K.Rows; n++)
                    //    {
                    //        int KL = K.RowIndex[n] + 1;
                    //        int KU = K.RowIndex[n + 1] - 1;
                    //        if (KU >= KL)
                    //        {
                    //            int k = n;
                    //            for (int KK = KL; KK <= KU; KK++)
                    //            {
                    //                k--;
                    //                tempVectorIn[n] += c.Item4 * d[KK] * tempVectorOut[k];
                    //            }
                    //        }
                    //    }
                    }
                    ////for (int n = 0; n < meanSize; n++)
                    ////    o[blockStart + n] -= tempVectorIn[n];

                    //for (int n = 0; n < meanSize; n++)
                    //    o[blockStart + n] = (o[blockStart + n] - tempVectorIn[n]) * ssorDiagonals[blockStart + n];
                }

                // Multiply with diagonal
                Array.Clear(oTemp.Data, 0, this.VectorSize);
                for (int i = 0; i < coeffLower.Count; i++)
                {
                    foreach (var c in coeffLower[i].Where(x => x.Item2 == x.Item3 && x.Item1 == 0))
                    {
                        SkylineMatrix2D k = (SkylineMatrix2D)analyzer.Matrices[c.Item1][solver.LinearSystem.ID];
                        int row = c.Item2 * k.Rows;
                        int col = c.Item3 * k.Rows;
                        k.Multiply(o, oTemp.Data, c.Item4, col, row, false);
                    }
                }

                // Backward sub
                o.CopyFrom(0, v.Length, oTemp, 0);
                //Array.Copy(oTemp.Data, 0, o.Data, 0, v.Data.Length);
                for (int i = coeffUpper.Count - 1; i >= 0; i--)
                {
                    int blockStart = i * meanSize;
                    Array.Copy(oTemp.Data, blockStart, tempVectorOut.Data, 0, meanSize);
                    //Array.Clear(tempVectorOut.Data, 0, meanSize);
                    foreach (var c in coeffUpper[i].Where(x => x.Item2 != x.Item3))
                    {
                        SkylineMatrix2D k = (SkylineMatrix2D)analyzer.Matrices[c.Item1][solver.LinearSystem.ID];
                        int row = c.Item2 * k.Rows;
                        int col = c.Item3 * k.Rows;
                        k.Multiply(o, tempVectorOut.Data, -c.Item4, col, 0, false);
                    }

                    foreach (var c in coeffUpper[i].Where(x => x.Item2 == x.Item3 && x.Item1 == 0))
                    {
                        //Array.Copy(oTemp.Data, i * meanSize, tempVectorIn.Data, 0, meanSize);
                        Array.Copy(tempVectorOut.Data, 0, tempVectorIn.Data, 0, meanSize);
                        tempVectorIn.Multiply(1.0 / c.Item4);
                        analyzer.FactorizedMatrices[0][1].Solve(tempVectorIn, tempVectorOut);
                        o.CopyFrom(0, meanSize, tempVectorOut, i * meanSize);
                        //Array.Copy(tempVectorOut.Data, 0, o.Data, i * meanSize, meanSize);

                        //SkylineMatrix2D<double> K = (SkylineMatrix2D<double>)analyzer.Matrices[c.Item1][subdomain.ID];
                        //double[] d = K.Data as double[];
                        //int n = K.Rows - 1;
                        //for (int l = 0; l < K.Rows; l++)
                        //{
                        //    int KL = K.RowIndex[n] + 1;
                        //    int KU = K.RowIndex[n + 1] - 1;
                        //    if (KU >= KL)
                        //    {
                        //        int k = n;
                        //        for (int KK = KL; KK <= KU; KK++)
                        //        {
                        //            k--;
                        //            tempVectorOut[k] += c.Item4 * d[KK] * tempVectorOut[n];
                        //            //result[k] -= d[KK] * result[n];
                        //        }
                        //    }
                        //    n--;
                        //}
                    }
                    //for (int n = 0; n < meanSize; n++)
                    //    o[blockStart + n] -= tempVectorOut[n];
                }

                //// Division
                //for (int i = 0; i < this.VectorSize; i++)
                //    o.Data[i] *= ssorDiagonals[i];
            //}
        }

        private void PreconditionSSOR(IVector vIn, IVector vOut)
        {
            if (analyzer.FactorizedMatrices.Count < 1)
                throw new InvalidOperationException("Cannot precondition with no factorized matrices");
            if (analyzer.FactorizedMatrices[0].Count != 1)
                throw new InvalidOperationException("Cannot precondition with more than one subdomains");

            MakeSSORDiagonals();
            //for (int i = 0; i < this.VectorSize; i++)
            //{
            //    vOut[i] = vIn[i] * ssorDiagonals[i];
            //}
            //return;

            if (tempVectorIn == null)
                tempVectorIn = new Vector(meanSize);
            if (tempVectorOut == null)
                tempVectorOut = new Vector(meanSize);

            var v = (Vector)vIn;
            var o = (Vector)vOut;
            var oTemp = new Vector(this.VectorSize);

            //foreach (ILinearSystem subdomain in solver.SubdomainsDictionary.Values)
            //{
                // Forward sub
                Array.Copy(v.Data, 0, o.Data, 0, v.Data.Length);
                for (int i = 0; i < coeffLower.Count; i++)
                {
                    int blockStart = i * meanSize;
                    Array.Copy(v.Data, blockStart, tempVectorOut.Data, 0, meanSize);
                    //Array.Clear(tempVectorIn.Data, 0, meanSize);
                    //Array.Clear(tempVectorOut.Data, 0, meanSize);

                    foreach (var c in coeffLower[i].Where(x => x.Item2 != x.Item3))
                    {
                        SkylineMatrix2D k = (SkylineMatrix2D)analyzer.Matrices[c.Item1][solver.LinearSystem.ID];
                        int row = c.Item2 * k.Rows;
                        int col = c.Item3 * k.Rows;
                        k.Multiply(o, tempVectorOut.Data, -c.Item4, col, 0, false);
                    }

                    Array.Clear(tempVectorIn.Data, 0, meanSize);
                    foreach (var c in coeffLower[i].Where(x => x.Item2 == x.Item3 && x.Item1 == 0))
                    {
                        ////Array.Copy(v.Data, i * meanSize, tempVectorIn.Data, 0, meanSize);
                        //Array.Copy(tempVectorOut.Data, 0, tempVectorIn.Data, 0, meanSize);
                        //tempVectorIn.Multiply(1.0 / c.Item4);
                        //analyzer.FactorizedMatrices[0][1].Solve(tempVectorIn, tempVectorOut.Data);
                        //Array.Copy(tempVectorOut.Data, 0, o.Data, i * meanSize, meanSize);

                        //double coefficientInverse = 1.0 / c.Item4;
                        //Array.Copy(tempVectorOut.Data, 0, tempVectorIn.Data, 0, meanSize);
                        //tempVectorIn.Multiply(1.0 / c.Item4);
                        tempVectorOut.Multiply(1.0 / c.Item4);
                        SkylineMatrix2D K = (SkylineMatrix2D)analyzer.Matrices[c.Item1][solver.LinearSystem.ID];
                        double[] d = K.Data as double[];
                        for (int n = 0; n < K.Rows; n++)
                        {
                            int KL = K.RowIndex[n] + 1;
                            int KU = K.RowIndex[n + 1] - 1;
                            if (KU >= KL)
                            {
                                int k = n;
                                double C = 0;
                                for (int KK = KL; KK <= KU; KK++)
                                {
                                    k--;
                                    C += d[KK] * tempVectorOut[k]; // c.Item4
                                    //tempVectorIn[n] += c.Item4 * d[KK] * tempVectorOut[k];
                                }
                                tempVectorOut[n] -= C;
                            }
                            tempVectorOut[n] /= d[KL - 1];
                        }
                        Array.Copy(tempVectorOut.Data, 0, o.Data, i * meanSize, meanSize);
                    }
                    ////for (int n = 0; n < meanSize; n++)
                    ////    o[blockStart + n] -= tempVectorIn[n];

                    //for (int n = 0; n < meanSize; n++)
                    //    o[blockStart + n] = (o[blockStart + n] - tempVectorIn[n]) * ssorDiagonals[blockStart + n];
                }

                // Multiply with diagonal
                for (int i = 0; i < v.Data.Length; i++)
                    oTemp.Data[i] = o.Data[i] * ssorDiagonals[i];
                //Array.Copy(o.Data, 0, oTemp.Data, 0, v.Data.Length);
                //Array.Clear(oTemp.Data, 0, this.VectorSize);
                //for (int i = 0; i < coeffLower.Count; i++)
                //{
                //    //foreach (var c in coeffLower[i].Where(x => x.Item2 == x.Item3 && x.Item1 == 0))
                //    //{
                //    //    SkylineMatrix2D<double> k = (SkylineMatrix2D<double>)analyzer.Matrices[c.Item1][subdomain.ID];
                //    //    int row = c.Item2 * k.Rows;
                //    //    int col = c.Item3 * k.Rows;
                //    //    k.Multiply(o, oTemp.Data, c.Item4, col, row, false);
                //    //}
                //    foreach (var c in coeffLower[i].Where(x => x.Item2 == x.Item3 && x.Item1 == 0))
                //    {
                //        SkylineMatrix2D<double> k = (SkylineMatrix2D<double>)analyzer.Matrices[c.Item1][subdomain.ID];
                //        int row = c.Item2 * k.Rows;
                //        int col = c.Item3 * k.Rows;

                //        for (int cnt = 0; cnt < k.RowIndex.Length - 1; cnt++)
                //        {
                //            oTemp.Data[row + cnt] = c.Item4 * o.Data[col + cnt] * k.Data[k.RowIndex[cnt]];
                //            //k.Multiply(o, oTemp.Data, c.Item4, col, row, false);
                //        }
                //    }
                //}

                // Backward sub
                Array.Copy(oTemp.Data, 0, o.Data, 0, v.Data.Length);
                for (int i = coeffUpper.Count - 1; i >= 0; i--)
                {
                    int blockStart = i * meanSize;
                    Array.Copy(oTemp.Data, blockStart, tempVectorOut.Data, 0, meanSize);
                    //Array.Clear(tempVectorOut.Data, 0, meanSize);
                    foreach (var c in coeffUpper[i].Where(x => x.Item2 != x.Item3))
                    {
                        SkylineMatrix2D k = (SkylineMatrix2D)analyzer.Matrices[c.Item1][solver.LinearSystem.ID];
                        int row = c.Item2 * k.Rows;
                        int col = c.Item3 * k.Rows;
                        k.Multiply(o, tempVectorOut.Data, -c.Item4, col, 0, false);
                    }

                    foreach (var c in coeffUpper[i].Where(x => x.Item2 == x.Item3 && x.Item1 == 0))
                    {
                        ////Array.Copy(oTemp.Data, i * meanSize, tempVectorIn.Data, 0, meanSize);
                        //Array.Copy(tempVectorOut.Data, 0, tempVectorIn.Data, 0, meanSize);
                        //tempVectorIn.Multiply(1.0 / c.Item4);
                        //analyzer.FactorizedMatrices[0][1].Solve(tempVectorIn, tempVectorOut.Data);
                        //Array.Copy(tempVectorOut.Data, 0, o.Data, i * meanSize, meanSize);

                        //Array.Copy(tempVectorOut.Data, 0, tempVectorIn.Data, 0, meanSize);
                        //tempVectorIn.Multiply(1.0 / c.Item4);
                        tempVectorOut.Multiply(1.0 / c.Item4);
                        SkylineMatrix2D K = (SkylineMatrix2D)analyzer.Matrices[c.Item1][solver.LinearSystem.ID];
                        //double[] d = K.Data as double[];
                        //for (int n = K.Rows - 1; n >= 0; n--)
                        //{
                        //    int KL = K.RowIndex[n] + 1;
                        //    int KU = K.RowIndex[n + 1] - 1;
                        //    if (KU >= KL)
                        //    {
                        //        int k = n;
                        //        double C = 0;
                        //        for (int KK = KL; KK <= KU; KK++)
                        //        {
                        //            k--;
                        //            C += d[KK] * tempVectorOut[k]; // c.Item4
                        //            //tempVectorIn[n] += c.Item4 * d[KK] * tempVectorOut[k];
                        //        }
                        //        tempVectorOut[n] -= C;
                        //    }
                        //    tempVectorOut[n] /= d[KL - 1];
                        //}
                        double[] d = K.Data as double[];
                        int n = K.Rows - 1;
                        for (int l = 0; l < K.Rows; l++)
                        {
                            int KL = K.RowIndex[n] + 1;
                            int KU = K.RowIndex[n + 1] - 1;
                            tempVectorOut[n] /= d[KL - 1];

                            if (KU >= KL)
                            {
                                int k = n;
                                for (int KK = KL; KK <= KU; KK++)
                                {
                                    k--;
                                    tempVectorOut[k] -= d[KK] * tempVectorOut[n];
                                    //result[k] -= d[KK] * result[n];
                                }
                            }
                            n--;
                        }
                        Array.Copy(tempVectorOut.Data, 0, o.Data, i * meanSize, meanSize);
                    }
                    //for (int n = 0; n < meanSize; n++)
                    //    o[blockStart + n] -= tempVectorOut[n];
                }

                //// Division
                //for (int i = 0; i < this.VectorSize; i++)
                //    o.Data[i] *= ssorDiagonals[i];
            //}
        }

        private void PreconditionILU(IVector vIn, IVector vOut)
        {
            if (analyzer.FactorizedMatrices.Count < 1)
                throw new InvalidOperationException("Cannot precondition with no factorized matrices");
            if (analyzer.FactorizedMatrices[0].Count != 1)
                throw new InvalidOperationException("Cannot precondition with more than one subdomains");

            //MakeSSORDiagonals();
            //for (int i = 0; i < this.VectorSize; i++)
            //{
            //    vOut[i] = vIn[i] * ssorDiagonals[i];
            //}
            //return;

            if (tempVectorIn == null)
                tempVectorIn = new Vector(meanSize);
            if (tempVectorOut == null)
                tempVectorOut = new Vector(meanSize);

            var v = (Vector)vIn;
            var o = (Vector)vOut;
            var oTemp = new Vector(this.VectorSize);
            var tOut = new Vector(meanSize);

            //foreach (ILinearSystem subdomain in solver.SubdomainsDictionary.Values)
            //{
                // Forward sub
                Array.Copy(v.Data, 0, o.Data, 0, v.Data.Length);
                for (int i = 0; i < coeffLower.Count; i++)
                {
                    int blockStart = i * meanSize;
                    Array.Copy(v.Data, blockStart, tempVectorOut.Data, 0, meanSize);
                    //Array.Clear(tempVectorIn.Data, 0, meanSize);
                    //Array.Clear(tempVectorOut.Data, 0, meanSize);

                    foreach (var c in coeffLower[i].Where(x => x.Item2 != x.Item3))
                    {
                        SkylineMatrix2D k = analyzer.FactorizedMatrices[c.Item1][solver.LinearSystem.ID];
                        int row = c.Item2 * k.Rows;
                        int col = c.Item3 * k.Rows;

                        //Array.Copy(o.Data, col, tempVectorIn.Data, 0, meanSize);
                        //tempVectorIn.Multiply(1.0 / -c.Item4);
                        //k.Solve(tempVectorIn, tOut.Data);
                        //for (int p = 0; p < meanSize; p++)
                        //    tempVectorOut[p] += tOut[p];

                        //int row = c.Item2 * k.Rows;
                        //int col = c.Item3 * k.Rows;
                        k.Multiply(o, tempVectorOut.Data, -c.Item4, col, 0, false);
                    }

                    Array.Clear(tempVectorIn.Data, 0, meanSize);
                    foreach (var c in coeffLower[i].Where(x => x.Item2 == x.Item3 && x.Item1 == 0))
                    {
                        //Array.Copy(v.Data, i * meanSize, tempVectorIn.Data, 0, meanSize);
                        Array.Copy(tempVectorOut.Data, 0, tempVectorIn.Data, 0, meanSize);
                        tempVectorIn.Multiply(1.0 / c.Item4);
                        analyzer.FactorizedMatrices[0][1].Solve(tempVectorIn, tempVectorOut);
                        Array.Copy(tempVectorOut.Data, 0, o.Data, i * meanSize, meanSize);

                        //    SkylineMatrix2D<double> K = (SkylineMatrix2D<double>)analyzer.Matrices[c.Item1][subdomain.ID];
                        //    double[] d = K.Data as double[];
                        //    for (int n = 0; n < K.Rows; n++)
                        //    {
                        //        int KL = K.RowIndex[n] + 1;
                        //        int KU = K.RowIndex[n + 1] - 1;
                        //        if (KU >= KL)
                        //        {
                        //            int k = n;
                        //            for (int KK = KL; KK <= KU; KK++)
                        //            {
                        //                k--;
                        //                tempVectorIn[n] += c.Item4 * d[KK] * tempVectorOut[k];
                        //            }
                        //        }
                        //    }
                    }
                    ////for (int n = 0; n < meanSize; n++)
                    ////    o[blockStart + n] -= tempVectorIn[n];

                    //for (int n = 0; n < meanSize; n++)
                    //    o[blockStart + n] = (o[blockStart + n] - tempVectorIn[n]) * ssorDiagonals[blockStart + n];
                }

                // Multiply with diagonal
                Array.Clear(oTemp.Data, 0, this.VectorSize);
                for (int i = 0; i < coeffLower.Count; i++)
                {
                    foreach (var c in coeffLower[i].Where(x => x.Item2 == x.Item3 && x.Item1 == 0))
                    {
                        SkylineMatrix2D k = (SkylineMatrix2D)analyzer.Matrices[c.Item1][solver.LinearSystem.ID];
                        int row = c.Item2 * k.Rows;
                        int col = c.Item3 * k.Rows;
                        k.Multiply(o, oTemp.Data, c.Item4, col, row, false);
                    }
                }

                // Backward sub
                Array.Copy(oTemp.Data, 0, o.Data, 0, v.Data.Length);
                for (int i = coeffUpper.Count - 1; i >= 0; i--)
                {
                    int blockStart = i * meanSize;
                    Array.Copy(oTemp.Data, blockStart, tempVectorOut.Data, 0, meanSize);
                    //Array.Clear(tempVectorOut.Data, 0, meanSize);
                    foreach (var c in coeffUpper[i].Where(x => x.Item2 != x.Item3))
                    {
                        SkylineMatrix2D k = (SkylineMatrix2D)analyzer.FactorizedMatrices[c.Item1][solver.LinearSystem.ID];
                        int row = c.Item2 * k.Rows;
                        int col = c.Item3 * k.Rows;
                        k.Multiply(o, tempVectorOut.Data, -c.Item4, col, 0, false);

                        //SkylineMatrix2D<double> k = analyzer.FactorizedMatrices[c.Item1][subdomain.ID];
                        //int row = c.Item2 * k.Rows;
                        //int col = c.Item3 * k.Rows;

                        //Array.Copy(o.Data, col, tempVectorIn.Data, 0, meanSize);
                        //tempVectorIn.Multiply(1.0 / -c.Item4);
                        //k.Solve(tempVectorIn, tOut.Data);
                        //for (int p = 0; p < meanSize; p++)
                        //    tempVectorOut[p] += tOut[p];
                    }

                    foreach (var c in coeffUpper[i].Where(x => x.Item2 == x.Item3 && x.Item1 == 0))
                    {
                        //Array.Copy(oTemp.Data, i * meanSize, tempVectorIn.Data, 0, meanSize);
                        Array.Copy(tempVectorOut.Data, 0, tempVectorIn.Data, 0, meanSize);
                        tempVectorIn.Multiply(1.0 / c.Item4);
                        analyzer.FactorizedMatrices[0][1].Solve(tempVectorIn, tempVectorOut);
                        Array.Copy(tempVectorOut.Data, 0, o.Data, i * meanSize, meanSize);

                        //SkylineMatrix2D<double> K = (SkylineMatrix2D<double>)analyzer.Matrices[c.Item1][subdomain.ID];
                        //double[] d = K.Data as double[];
                        //int n = K.Rows - 1;
                        //for (int l = 0; l < K.Rows; l++)
                        //{
                        //    int KL = K.RowIndex[n] + 1;
                        //    int KU = K.RowIndex[n + 1] - 1;
                        //    if (KU >= KL)
                        //    {
                        //        int k = n;
                        //        for (int KK = KL; KK <= KU; KK++)
                        //        {
                        //            k--;
                        //            tempVectorOut[k] += c.Item4 * d[KK] * tempVectorOut[n];
                        //            //result[k] -= d[KK] * result[n];
                        //        }
                        //    }
                        //    n--;
                        //}
                    }
                    //for (int n = 0; n < meanSize; n++)
                    //    o[blockStart + n] -= tempVectorOut[n];
                }

                //// Division
                //for (int i = 0; i < this.VectorSize; i++)
                //    o.Data[i] *= ssorDiagonals[i];
            //}
        }

        //private void PreconditionILU(IVector<double> vIn, IVector<double> vOut)
        //{
        //    if (analyzer.FactorizedMatrices.Count < 1)
        //        throw new InvalidOperationException("Cannot precondition with no factorized matrices");
        //    if (analyzer.FactorizedMatrices[0].Count != 1)
        //        throw new InvalidOperationException("Cannot precondition with more than one subdomains");

        //    if (tempVectorIn == null)
        //        tempVectorIn = new Vector<double>(meanSize);
        //    if (tempVectorOut == null)
        //        tempVectorOut = new Vector<double>(meanSize);

        //    var v = (Vector<double>)vIn;
        //    var o = (Vector<double>)vOut;

        //    Array.Clear(o.Data, 0, o.Length);
        //    for (int pos = 0; pos < (isGaussian ? mOrder + 1 : calculator.PsiSize); pos++)
        //    {
        //        var m = analyzer.FactorizedMatrices[pos][1];
        //        if (isGaussian)
        //            foreach (var x in calculator.GaussianCoefficients[pos])
        //            {
        //                int i = x.Item1;
        //                int j = x.Item2;

        //                Array.Copy(v.Data, j * meanSize, tempVectorIn.Data, 0, meanSize);
        //                tempVectorIn.Multiply(1.0 / x.Item3);
        //                m.Solve(tempVectorIn, tempVectorOut.Data);

        //                for (int p = 0; p < meanSize; p++)
        //                    o.Data[i * meanSize + p] += tempVectorOut.Data[p];
        //            }
        //        else
        //            foreach (var x in calculator.LognormalCoefficients[pos])
        //            {
        //                int i = x.Key.Item1;
        //                int j = x.Key.Item2;
        //                Array.Copy(v.Data, j * meanSize, tempVectorIn.Data, 0, meanSize);
        //                tempVectorIn.Multiply(1.0 / x.Value);
        //                m.Solve(tempVectorIn, tempVectorOut.Data);

        //                for (int p = 0; p < meanSize; p++)
        //                    o.Data[i * meanSize + p] += tempVectorOut.Data[p];
        //            }
        //    }
        //}

        public void MultiplyWithMatrix(IVector vIn, IVector vOut)
        {
            var data = ((Vector)vOut).Data;
            Array.Clear(data, 0, data.Length);

            //foreach (ILinearSystem subdomain in solver.SubdomainsDictionary.Values)
            //{
                for (int pos = 0; pos < (isGaussian ? mOrder + 1 : calculator.PsiSize); pos++)
                {
                    //SkylineMatrix2D<double> k = ((SkylineMatrix2D<double>)subdomain.Matrix);
                    SkylineMatrix2D k = (SkylineMatrix2D)analyzer.Matrices[pos][solver.LinearSystem.ID];
                    if (isGaussian)
                        foreach (var x in calculator.GaussianCoefficients[pos])
                        {
                            int i = x.Item1 * k.Rows;
                            int j = x.Item2 * k.Rows;
                            k.Multiply(vIn, data, x.Item3, j, i, false);
                        }
                    else
                        foreach (var x in calculator.LognormalCoefficients[pos])
                        {
                            int i = x.Key.Item1 * k.Rows;
                            int j = x.Key.Item2 * k.Rows;
                            k.Multiply(vIn, data, x.Value, j, i, false);
                        }
                }
            //}
        }
    }
}
