using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Runtime;
using Accord.Math.Decompositions;
using ISAAR.MSolve.PreProcessor;
using MGroup.Stochastic.Interfaces;
using Troschuetz.Random.Distributions.Continuous;


namespace MGroup.Stochastic.Structural.StochasticRealizers
{
    public class KarhunenLoeveCoefficientsProvider : IUncertainParameterRealizer
    {
        public static double GaussianKernelCovarianceFunction(double x, double y, double sigmaSquare, double correlationLength)
        {
            double nominator = -Math.Abs(x - y) / correlationLength;
            double correlationFunction = Math.Pow(Math.E, nominator) * sigmaSquare;
            return correlationFunction;
        }

        public double[] Realize(int iteration, double[] parameters)
        {
            throw new NotImplementedException();
        }

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
            double[] GaussLegendreCoordinates = KarhunenLoeveCoefficientsProvider.gauss_quad().Item1;
            double[] GaussLegendreWeights = KarhunenLoeveCoefficientsProvider.gauss_quad().Item2;
            double[,] Bmatrix = new double[ndof, ndof];
            for (int i = 0; i < nfe; i++)
            {
                double[,] Be = new double[neq, neq];
                double det_Je = (xCoordinates[IEN[1, i]] - xCoordinates[IEN[0, i]]) / 2;
                for (int j = 0; j < GaussLegendreOrder; j++)
                {
                    double xi_gl = GaussLegendreCoordinates[j];
                    double w_gl = GaussLegendreWeights[j];
                    double[] NN = KarhunenLoeveCoefficientsProvider.LagrangianShapeFunctions(xi_gl);
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
                        double[] NNe = KarhunenLoeveCoefficientsProvider.LagrangianShapeFunctions(xi_gl_e);
                        double xpk = NNe[0] * xe[0] + NNe[1] * xe[1];
                        for (int l = 0; l < GaussLegendreOrder; l++)
                        {
                            double xi_gl_f = GaussLegendreCoordinates[l];
                            double[] NNf = KarhunenLoeveCoefficientsProvider.LagrangianShapeFunctions(xi_gl_f);
                            double xpl = NNf[0] * xf[0] + NNf[1] * xf[1];
                            //element C matrix
                            Cef[0, 0] = Cef[0, 0] + KarhunenLoeveCoefficientsProvider.GaussianKernelCovarianceFunction(xpk, xpl, sigmaSquare, correlationLength) * NNe[0] * NNf[0] * det_Je * det_Jf * GaussLegendreWeights[k] * GaussLegendreWeights[l];
                            Cef[1, 0] = Cef[1, 0] + KarhunenLoeveCoefficientsProvider.GaussianKernelCovarianceFunction(xpk, xpl, sigmaSquare, correlationLength) * NNe[0] * NNf[1] * det_Je * det_Jf * GaussLegendreWeights[k] * GaussLegendreWeights[l];
                            Cef[0, 1] = Cef[0, 1] + KarhunenLoeveCoefficientsProvider.GaussianKernelCovarianceFunction(xpk, xpl, sigmaSquare, correlationLength) * NNe[1] * NNf[0] * det_Je * det_Jf * GaussLegendreWeights[k] * GaussLegendreWeights[l];
                            Cef[1, 1] = Cef[1, 1] + KarhunenLoeveCoefficientsProvider.GaussianKernelCovarianceFunction(xpk, xpl, sigmaSquare, correlationLength) * NNe[1] * NNf[1] * det_Je * det_Jf * GaussLegendreWeights[k] * GaussLegendreWeights[l];
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

        public double[,] KarhunenLoeveFredholm1DSampleGenerator(int MCsamples, double[] eigenValues, double[,] eigenModes, double meanValue, bool midpointMethod, bool isGaussian)
        {
            if (midpointMethod == false) throw new ArgumentException("It is not supported at the moment");
            if (isGaussian == false) throw new ArgumentException("It is not supported at the moment");
            double[,] eigenModesAtMidpoint = new double[eigenModes.GetLength(0) - 1, eigenModes.GetLength(1)];
            for (int j = 0; j < eigenModes.GetLength(1); j++)
            {
                for (int i = 0; i < eigenModes.GetLength(0) - 1; i++)
                {
                    eigenModesAtMidpoint[i, j] = eigenModes[i, j] + eigenModes[i + 1, j];
                }
            }
            double[,] fieldRealizations = new double[MCsamples, eigenModesAtMidpoint.GetLength(0)];
            for (int k = 0; k < MCsamples; k++)
            {
                double[] ksi = new double[eigenModesAtMidpoint.GetLength(1)];
                for (int j = 0; j < eigenModesAtMidpoint.GetLength(1); j++)
                {
                    var KsiNormalDistribution = new NormalDistribution(0, 1);
                    ksi[j] = KsiNormalDistribution.NextDouble();
                }
                for (int i = 0; i < eigenModesAtMidpoint.GetLength(0); i++)
                {
                    for (int j = 0; j < eigenModesAtMidpoint.GetLength(1); j++)
                    {
                        fieldRealizations[k, i] = fieldRealizations[k, i] + Math.Sqrt(eigenValues[j]) * eigenModesAtMidpoint[i, j] * ksi[j];
                    }
                    fieldRealizations[k, i] = meanValue + fieldRealizations[k, i];
                }
            }

            return fieldRealizations;
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
