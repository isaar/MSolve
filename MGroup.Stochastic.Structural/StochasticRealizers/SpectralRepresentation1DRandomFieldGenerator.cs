using ISAAR.MSolve.FEM.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Materials.Interfaces;
using MGroup.Stochastic.Interfaces;
using Troschuetz.Random.Distributions.Continuous;

namespace MGroup.Stochastic.Structural.StochasticRealizers
{
    public class SpectralRepresentation1DRandomFieldGenerator : IUncertainParameterRealizer
    {
        private const double tolerance = 1e-10;
        private readonly double b;
        private readonly double spectrumStandardDeviation, cutoffError, frequencyIncrement;
        private readonly int frequencyIntervals;
        private double[] sffTarget, omegas;
        private double wu, period;
        private int nptsSff, frequencyCounter;
        private double[] randomVariables = new double[0];
        private double[] phi;
        public double MeanValue;
        private bool ResetGeneration = true;
        private int PreviousIteration = -1;

        /// <summary>A class implementing the Spectral Respresentation methodology
        /// using Fourier series that generates 1D stochastic fields with user selected correlation structure parameters.</summary>
        /// <param name="b">The b.</param>
        /// <param name="spectrumStandardDeviation">The spectrum standard deviation.</param>
        /// <param name="meanValue">The mean value.</param>
        /// <param name="cutoffError">The cutoff error.</param>
        /// <param name="frequencyIncrement">The frequency increment.</param>
        /// <param name="frequencyIntervals">The frequency intervals.</param>
        public SpectralRepresentation1DRandomFieldGenerator(double b, double spectrumStandardDeviation, double meanValue, double cutoffError,
            double frequencyIncrement = 0.1, int frequencyIntervals = 256)
        {
            this.b = b;
            this.spectrumStandardDeviation = spectrumStandardDeviation;
            MeanValue = meanValue;
            this.cutoffError = cutoffError;
            this.frequencyIncrement = frequencyIncrement;
            this.frequencyIntervals = frequencyIntervals;
            Calculate();
        }

        public double SpectrumStandardDeviation { get { return spectrumStandardDeviation; } }
        public double[] SffTarget { get { return sffTarget; } }
        public double[] Omegas { get { return omegas; } }
        public double Wu { get { return wu; } }
        public double Period { get { return period; } }
        public int NPtsSff { get { return nptsSff; } }
        public int FrequencyIntervals { get { return frequencyIntervals; } }
        public int CurrentMCS { get; set; }
        public int CurrentFrequency { get; set; }

        /// <summary>Calculates method intrinsics.</summary>
        private void Calculate()
        {
            double integral = AutoCorrelation(0) / 2d;

            frequencyCounter = 1;
            double cumulativeSum = 0;
            double trapezoidArea = 0;

            while (cumulativeSum < (1d - cutoffError) * integral)
            {
                trapezoidArea = 0.5 * (SpectralDensity((frequencyCounter - 1) * frequencyIncrement) + SpectralDensity(frequencyCounter * frequencyIncrement)) * frequencyIncrement;
                if (trapezoidArea < tolerance)
                    break;

                cumulativeSum += trapezoidArea;
                frequencyCounter++;
            }

            wu = frequencyCounter * frequencyIncrement;
            sffTarget = new double[frequencyIntervals];
            omegas = new double[frequencyIntervals];
            double dw = wu / (double)frequencyIntervals;

            for (int i = 0; i < frequencyIntervals; i++)
            {
                omegas[i] = dw / 2 + i * dw;
                sffTarget[i] = SpectralDensity(dw / 2 + i * dw);
            }
            period = 2d * Math.PI / dw;
            nptsSff = (int)(period * wu / Math.PI);
        }

        private double AutoCorrelation(double tau)
        {
            return Math.Pow(spectrumStandardDeviation, 2) * Math.Exp(-Math.Abs(tau) / b);
        }

        private double SpectralDensity(double omega)
        {
            return .5 / Math.PI * Math.Pow(spectrumStandardDeviation, 2) * Math.Sqrt(Math.PI * b) *
                Math.Exp(-.25 * b * Math.Pow(omega, 2));
        }


        /// <summary>Realizes a sample function for respective stochastic domain mapper and parameters.</summary>
        /// <param name="iteration">The iteration.</param>
        /// <param name="domainMapper">The domain mapper.</param>
        /// <param name="parameters">The parameters.</param>
        /// <returns></returns>
        public double Realize(int iteration, IStochasticDomainMapper domainMapper, double[] parameters)
        {
            ResetGeneration = (PreviousIteration != iteration);
            if (ResetGeneration) ResetSampleGeneration();
            double[] stochasticDomainPoint = domainMapper.Map(parameters);
            var dw = wu / (double)frequencyIntervals;
            int i = 0;
            double randomCoefficient = 0;
            while (i < frequencyIntervals)
            {
                randomCoefficient += Math.Sqrt(2) * Math.Sqrt(2 * SpectralDensity(dw / 2 + i * dw) * dw) *
                                     (Math.Cos((dw / 2 + i * dw) * stochasticDomainPoint[0] +
                                     phi[i]));
                i++;
            }

            if (randomCoefficient >= .9)
            {
                randomCoefficient = .9;
            }
            else if (randomCoefficient <= -.9)
            {
                randomCoefficient = -.9;
            }

            PreviousIteration = iteration;
            return MeanValue / (1 + randomCoefficient);
        }

        public void ResetSampleGeneration()
        {
            phi = new double[frequencyIntervals];
            var Phi = new ContinuousUniformDistribution(0d, 2d * Math.PI);
            for (int i = 0; i < frequencyIntervals; i++)
            {
                phi[i] = Phi.NextDouble();
            }
        }

        //public double[] RandomVariables
        //{
        //    get
        //    {

        //        return randomVariables;
        //    }
        //    set
        //    {
        //        phi = new double[frequencyIntervals];
        //        var Phi = new ContinuousUniformDistribution(0d, 2d * Math.PI);
        //        for (int i = 0; i < frequencyIntervals; i++)
        //        {
        //            phi[i] = Phi.NextDouble();
        //        }
        //        //changedVariables = true;
        //        randomVariables = value;
        //    }
        //}
    }
}
