using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Runtime;
using Accord.Math;
using Accord.Math.Decompositions;
using ISAAR.MSolve.PreProcessor;
using MGroup.Stochastic.Interfaces;
using Troschuetz.Random.Distributions.Continuous;


namespace MGroup.Stochastic.Structural.StochasticRealizers
{
    public class KarhunenLoeveCoefficientsProvider : IUncertainParameterRealizer
    {
        public int MCsamples { get; }
        public int Partition { get; }
        public double MeanValue { get; }
        public bool MidpointMethod { get; }
        public bool IsGaussian { get; }
        public int KarLoeveTerms { get; }
        public double[] DomainBounds { get; }
        public double SigmaSquare { get; set; }
        public double CorrelationLength { get; set; }
        public double[] Xcoordinates { get; set; }
        public double[] Lambda { get; set; }
        public double[,] Eigenvectors { get; set; }
        public double[,] EigenModesAtPoint { get; set; }
        private NormalDistribution KseNormalDistribution;
        private double[] Kse;
        private bool ResetGeneration = true;
        private int PreviousIteration = -1;

        /// <summary>Initializes a new instance of the <see cref="KarhunenLoeveCoefficientsProvider"/> class.</summary>
        /// <param name="partition">The partition.</param>
        /// <param name="meanValue">The mean value.</param>
        /// <param name="midpointMethod">if set to <c>true</c> [midpoint method].</param>
        /// <param name="isGaussian">if set to <c>true</c> [is gaussian].</param>
        /// <param name="karLoeveTerms">The kar loeve terms.</param>
        /// <param name="domainBounds">The domain bounds.</param>
        /// <param name="sigmaSquare">The sigma square.</param>
        /// <param name="correlationLength">Length of the correlation.</param>
        public KarhunenLoeveCoefficientsProvider (int partition, double meanValue, bool midpointMethod, bool isGaussian, int karLoeveTerms,
            double[] domainBounds, double sigmaSquare, double correlationLength)
        {
            Partition = partition;
            MeanValue = meanValue;
            MidpointMethod = midpointMethod;
            IsGaussian = isGaussian;
            KarLoeveTerms = karLoeveTerms;
            DomainBounds = domainBounds;
            SigmaSquare = sigmaSquare;
            CorrelationLength = correlationLength;
            double[] xCoordinates = KarhunenLoeveFredholmWithFEM(KarLoeveTerms, DomainBounds, SigmaSquare, Partition, CorrelationLength).Item1;
            double[] lambda = KarhunenLoeveFredholmWithFEM(KarLoeveTerms, DomainBounds, SigmaSquare, Partition, CorrelationLength).Item2;
            double[,] eigenvectors = KarhunenLoeveFredholmWithFEM(KarLoeveTerms, DomainBounds, SigmaSquare, Partition, CorrelationLength).Item3;
            KseNormalDistribution = new NormalDistribution(27644437, 0, 1);
            Xcoordinates = xCoordinates;
            Lambda = lambda;
            Eigenvectors = eigenvectors;
        }


        /// <summary>Covariance function.</summary>
        /// <param name="x">The x.</param>
        /// <param name="y">The y.</param>
        /// <param name="sigmaSquare">The sigma square.</param>
        /// <param name="correlationLength">Length of the correlation.</param>
        /// <returns></returns>
        public double GaussianKernelCovarianceFunction(double x, double y, double sigmaSquare, double correlationLength)
        {
            //CorrelationLength = correlationLength;
            //SigmaSquare = sigmaSquare;
            double nominator = -Math.Abs(x - y) / correlationLength;
            double correlationFunction = Math.Pow(Math.E, nominator) * sigmaSquare;
            return correlationFunction;
        }

        /// <summary>Resets the sample generation.</summary>
        private void ResetSampleGeneration()
        {            
            Kse = new double[KarLoeveTerms];
            for (int i = 0; i < KarLoeveTerms; i++)
            {
                Kse[i] = KseNormalDistribution.NextDouble();
            }
        }

        /// <summary>Realizes the specified iteration based on a stochastic domain mapper and the respective domain parameters.</summary>
        /// <param name="iteration">The iteration.</param>
        /// <param name="domainMapper">The domain mapper.</param>
        /// <param name="parameters">The parameters.</param>
        /// <returns></returns>
        public double Realize(int iteration, IStochasticDomainMapper domainMapper, double[] parameters)
        {
            ResetGeneration = (PreviousIteration != iteration);
            if (ResetGeneration) ResetSampleGeneration();
            var stochasticDomainPoint = domainMapper.Map(parameters);
            double[] eigenModesAtPoint = CalculateEigenmodesAtPoint(Xcoordinates, Eigenvectors, stochasticDomainPoint[0]);
            var value = KarhunenLoeveFredholm1DSampleGenerator(stochasticDomainPoint, Lambda, eigenModesAtPoint, MeanValue, MidpointMethod, IsGaussian);
            PreviousIteration = iteration;
            return value;
        }

        /// <summary>  Calculates the Fredholm integral that provides eigenvalues and eigenfunctions of the expansion.</summary>
        /// <param name="KarLoeveTerms">The kar loeve terms.</param>
        /// <param name="domainBounds">The domain bounds.</param>
        /// <param name="sigmaSquare">The sigma square.</param>
        /// <param name="partition">The partition.</param>
        /// <param name="correlationLength">Length of the correlation.</param>
        /// <returns></returns>
        public Tuple<double[], double[], double[,]> KarhunenLoeveFredholmWithFEM(int KarLoeveTerms, double[] domainBounds, double sigmaSquare, int partition, double correlationLength)
        {
            //FEM parameters 
            int ned = 1; //number of dof per node
            int nen = 2; //number of nodes per element
            int nnp = partition; //number of nodal points
            int nfe = nnp - 1; //number of finite elements
            int neq = ned * nen; //number of element equations
            int ndof = ned * nnp; // number of degrees of freedom
            int GaussLegendreOrder = 3;

            double[] xCoordinates = new double[partition];
            int[,] IEN = new int[2, nfe];  //local to global
            int[] ID = new int[partition];

            for (int i = 0; i < partition; i++)
            {
                //LagrangianShapeFunctons();
                xCoordinates[i] = domainBounds[0] + (domainBounds[1] - domainBounds[0]) / nfe * i;
                ID[i] = i;
            }

            for (int i = 0; i < nfe; i++)
            {
                IEN[0, i] = i;
                IEN[1, i] = i + 1;
            }

            // localization matrix
            int[,] LM = new int[partition - 1, 2];
            for (int i = 0; i < partition - 1; i++)
            {
                LM[i, 0] = i;
                LM[i, 1] = i + 1;
            }

            // computing B matrix
            double[] GaussLegendreCoordinates = gauss_quad().Item1;
            double[] GaussLegendreWeights = gauss_quad().Item2;
            double[,] Bmatrix = new double[ndof, ndof];
            for (int i = 0; i < nfe; i++)
            {
                double[,] Be = new double[neq, neq];
                double det_Je = (xCoordinates[IEN[1, i]] - xCoordinates[IEN[0, i]]) / 2;
                for (int j = 0; j < GaussLegendreOrder; j++)
                {
                    double xi_gl = GaussLegendreCoordinates[j];
                    double w_gl = GaussLegendreWeights[j];
                    double[] NN = LagrangianShapeFunctions(xi_gl);
                    //Be=Be+NN'*NN*det_Je*w_gl(j)
                    Be[0, 0] = Be[0, 0] + NN[0] * NN[0] * det_Je * w_gl;
                    Be[1, 0] = Be[1, 0] + NN[0] * NN[1] * det_Je * w_gl;
                    Be[0, 1] = Be[0, 1] + NN[1] * NN[0] * det_Je * w_gl;
                    Be[1, 1] = Be[1, 1] + NN[1] * NN[1] * det_Je * w_gl;
                }
                Bmatrix[LM[i, 0], LM[i, 0]] = Bmatrix[LM[i, 0], LM[i, 0]] + Be[0, 0];
                Bmatrix[LM[i, 1], LM[i, 0]] = Bmatrix[LM[i, 1], LM[i, 0]] + Be[1, 0];
                Bmatrix[LM[i, 0], LM[i, 1]] = Bmatrix[LM[i, 0], LM[i, 1]] + Be[0, 1];
                Bmatrix[LM[i, 1], LM[i, 1]] = Bmatrix[LM[i, 1], LM[i, 1]] + Be[1, 1];
            }

            // computing C matrix
            double[,] Cmatrix = new double[ndof, ndof];
            for (int i = 0; i < nfe; i++)
            {
                double[] xe = { xCoordinates[IEN[0, i]], xCoordinates[IEN[1, i]] };
                double det_Je = (xCoordinates[IEN[1, i]] - xCoordinates[IEN[0, i]]) / 2;
                for (int j = 0; j < nfe; j++)
                {
                    double[,] Cef = new double[neq, neq];
                    double[] xf = { xCoordinates[IEN[0, j]], xCoordinates[IEN[1, j]] };
                    double det_Jf = (xCoordinates[IEN[1, j]] - xCoordinates[IEN[0, j]]) / 2;
                    for (int k = 0; k < GaussLegendreOrder; k++)
                    {
                        double xi_gl_e = GaussLegendreCoordinates[k];
                        double[] NNe = LagrangianShapeFunctions(xi_gl_e);
                        double xpk = NNe[0] * xe[0] + NNe[1] * xe[1];
                        for (int l = 0; l < GaussLegendreOrder; l++)
                        {
                            double xi_gl_f = GaussLegendreCoordinates[l];
                            double[] NNf = LagrangianShapeFunctions(xi_gl_f);
                            double xpl = NNf[0] * xf[0] + NNf[1] * xf[1];
                            //element C matrix
                            Cef[0, 0] = Cef[0, 0] + GaussianKernelCovarianceFunction(xpk, xpl, sigmaSquare, correlationLength) * NNe[0] * NNf[0] * det_Je * det_Jf * GaussLegendreWeights[k] * GaussLegendreWeights[l];
                            Cef[1, 0] = Cef[1, 0] + GaussianKernelCovarianceFunction(xpk, xpl, sigmaSquare, correlationLength) * NNe[0] * NNf[1] * det_Je * det_Jf * GaussLegendreWeights[k] * GaussLegendreWeights[l];
                            Cef[0, 1] = Cef[0, 1] + GaussianKernelCovarianceFunction(xpk, xpl, sigmaSquare, correlationLength) * NNe[1] * NNf[0] * det_Je * det_Jf * GaussLegendreWeights[k] * GaussLegendreWeights[l];
                            Cef[1, 1] = Cef[1, 1] + GaussianKernelCovarianceFunction(xpk, xpl, sigmaSquare, correlationLength) * NNe[1] * NNf[1] * det_Je * det_Jf * GaussLegendreWeights[k] * GaussLegendreWeights[l];
                        }
                    }
                    Cmatrix[LM[i, 0], LM[j, 0]] = Cmatrix[LM[i, 0], LM[j, 0]] + Cef[0, 0];
                    Cmatrix[LM[i, 0], LM[j, 1]] = Cmatrix[LM[i, 0], LM[j, 1]] + Cef[0, 1];
                    Cmatrix[LM[i, 1], LM[j, 0]] = Cmatrix[LM[i, 1], LM[j, 0]] + Cef[1, 0];
                    Cmatrix[LM[i, 1], LM[j, 1]] = Cmatrix[LM[i, 1], LM[j, 1]] + Cef[1, 1];
                }

            }
            bool sort = true;
            var gevd = new GeneralizedEigenvalueDecomposition(Cmatrix, Bmatrix, sort);
            double[] lambdaAll = gevd.RealEigenvalues;
            double[] lambda = lambdaAll.Skip(0).Take(KarLoeveTerms).ToArray();
            double[,] EigenvectorsAll = gevd.Eigenvectors;
            double[,] Eigenvectors = new double[partition, KarLoeveTerms];
            for (int i = 0; i < partition; i++)
            {
                for (int j = 0; j < KarLoeveTerms; j++)
                {
                    Eigenvectors[i, j] = EigenvectorsAll[i, j];  //each column corresponds to an eigenvector
                }
            }
            return new Tuple<double[], double[], double[,]>(xCoordinates, lambda, Eigenvectors);

        }

        /// <summary>Calculates the eigenmodes at specified point.</summary>
        /// <param name="xCoordinates">The x coordinates.</param>
        /// <param name="eigenModes">The eigen modes.</param>
        /// <param name="stochasticDomainPoint">The stochastic domain point.</param>
        /// <returns></returns>
        public double[] CalculateEigenmodesAtPoint(double[] xCoordinates, double[,] eigenModes, double stochasticDomainPoint)
        {
            List<double> xCoordinatesList = xCoordinates.ToList();
            double firstXcoordinate = xCoordinatesList.OrderBy(item => Math.Abs(stochasticDomainPoint - item)).First();
            int indexOfFirstEigenmodeValue = xCoordinatesList.IndexOf(firstXcoordinate);
            List<double> list = xCoordinatesList;
            list.Remove(firstXcoordinate);
            xCoordinatesList = xCoordinates.ToList();
            double secondXcoordinate = list.OrderBy(item => Math.Abs(stochasticDomainPoint - item)).First();
            int indexOfSecondEigenmodeValue = xCoordinatesList.IndexOf(secondXcoordinate);
            double[] eigenmodesAtPoint = new double[eigenModes.GetLength(1)];

            for (int j = 0; j < eigenModes.GetLength(1); j++)
            {
                    eigenmodesAtPoint[j] = (eigenModes[indexOfSecondEigenmodeValue, j] - eigenModes[indexOfFirstEigenmodeValue, j]) /
                                           (xCoordinates[indexOfSecondEigenmodeValue] - xCoordinates[indexOfFirstEigenmodeValue]) * 
                    (stochasticDomainPoint-xCoordinates[indexOfFirstEigenmodeValue]) + eigenModes[indexOfFirstEigenmodeValue,j];
            }
            return eigenmodesAtPoint;
        }

        /// <summary>  Field realization based on KL expansion.</summary>
        /// <param name="stochasticDomainPoint">The stochastic domain point.</param>
        /// <param name="eigenValues">The eigen values.</param>
        /// <param name="eigenmodesAtPoint">The eigenmodes at point.</param>
        /// <param name="meanValue">The mean value.</param>
        /// <param name="midpointMethod">if set to <c>true</c> [midpoint method].</param>
        /// <param name="isGaussian">if set to <c>true</c> [is gaussian].</param>
        /// <returns></returns>
        /// <exception cref="System.ArgumentException">It is not supported at the moment
        /// or
        /// It is not supported at the moment</exception>
        public double KarhunenLoeveFredholm1DSampleGenerator(double[] stochasticDomainPoint, double[] eigenValues, double[] eigenmodesAtPoint, double meanValue, bool midpointMethod, bool isGaussian)
        {
            if (midpointMethod == false) throw new ArgumentException("It is not supported at the moment");
            if (isGaussian == false) throw new ArgumentException("It is not supported at the moment");
            double fieldRealization = 0;
                      
            for (int j = 0; j < eigenmodesAtPoint.Length; j++)
            {
                fieldRealization = fieldRealization + Math.Sqrt(eigenValues[j]) * eigenmodesAtPoint[j] * Kse[j];
            }

            if (fieldRealization >= .9)
            {
                fieldRealization = .9;
            }
            else if (fieldRealization <= -.9)
            {
                fieldRealization = -.9;
            }

            fieldRealization = meanValue * (1 + fieldRealization);
            return fieldRealization;
        }

        private static double[] LagrangianShapeFunctions(double xi)
        {
            double[] NShape = { -(xi - 1) / 2, (xi + 1) / 2 };
            return NShape;
        }

        private static Tuple<double[], double[]> gauss_quad()
        {
            double[] GaussLegendreCoordinates = { -0.7746, 0, 0.7746 };
            double[] GaussLegendreWeights = { 0.5556, 0.8889, 0.5556 };
            return new Tuple<double[], double[]>(GaussLegendreCoordinates, GaussLegendreWeights);
        }
    }
}
