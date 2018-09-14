using ISAAR.MSolve.FEM.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Materials.Interfaces;
using Troschuetz.Random;
using Troschuetz.Random.Distributions.Continuous;

namespace ISAAR.MSolve.Analyzers
{
    public class PowerSpectrumTargetEvaluatorCoefficientsProvider : IStochasticMaterialCoefficientsProvider, IStochasticCoefficientsProvider
    {
        private readonly double b;
        private readonly double spectrumStandardDeviation, cutoffError, frequencyIncrement, tolerance;
        private readonly int nPhi, frequencyIntervals, nptsVRF;
        private readonly DOFType randomFieldDirection;
        private double[] sffTarget, omegas;
        private double wu, period;
        private int nptsSff, frequencyCounter;
        private double[] randomVariables = new double[0];
        private double[] phi;
        
        public PowerSpectrumTargetEvaluatorCoefficientsProvider(double b, double spectrumStandardDeviation, double cutoffError, int nPhi, int nptsVRF, DOFType randomFieldDirection,
            double frequencyIncrement = 0.1, int frequencyIntervals = 256, double tolerance = 1e-10)
        {
            this.b = b;
            this.spectrumStandardDeviation = spectrumStandardDeviation;
            this.cutoffError = cutoffError;
            this.nPhi = nPhi;
            this.randomFieldDirection = randomFieldDirection;
            this.frequencyIncrement = frequencyIncrement;
            this.frequencyIntervals = frequencyIntervals;
            this.nptsVRF = nptsVRF;
            this.tolerance = tolerance;
            Calculate();
        }

        public double SpectrumStandardDeviation { get { return spectrumStandardDeviation; } }
        public double[] SffTarget { get { return sffTarget; } }
        public double[] Omegas { get { return omegas; } }
        public double Wu { get { return wu; } }
        public double Period { get { return period; } }
        public int NPtsSff { get { return nptsSff; } }
        public int NPtsVRF { get { return nptsVRF; } }
        public int FrequencyIntervals { get { return frequencyIntervals; } }
        public int NPhi { get { return nPhi; } }
        public int CurrentMCS { get; set; }
        public int CurrentFrequency { get; set; }

        private double AutoCorrelation(double tau)
        {  
            return Math.Pow(spectrumStandardDeviation, 2) * Math.Exp(-Math.Abs(tau)/b);
        }

        private double SpectralDensity(double omega)
        {
            return .5/Math.PI * Math.Pow(spectrumStandardDeviation, 2) * Math.Sqrt(Math.PI * b) *
                Math.Exp(-.25 * b * Math.Pow(omega, 2));
            //Math.Pow(spectrumStandardDeviation, 2) * b / (Math.PI * (1 + Math.Pow(b, 2) * Math.Pow(omega, 2))));
        }
        
        private void Calculate()
        {
            //const int N = 128;
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

            //omega = dw/2:dw:wu;
            //Sff_target = feval(spectralDensityFunction,omega);
            for (int i = 0; i < frequencyIntervals; i++)
            {
                omegas[i] = dw / 2 + i * dw;
                sffTarget[i] = SpectralDensity(dw / 2 + i * dw);
            }
            period = 2d * Math.PI / dw;
            nptsSff = (int)(period * wu / Math.PI);            
        }

        public double GetCoefficient(double meanValue, double[] coordinates)
        {
            //var Sff = sffTarget;
            //var ku = wu;
            //var Std = spectrumStandardDeviation;
            //var T = period;
            //var n_wu = 20;
            //var Npts = 30;
            var dw = wu / (double)frequencyIntervals;
            int i = 0;
            double randomCoefficient = 0;
            while (i < frequencyIntervals)
            {
                randomCoefficient += Math.Sqrt(2) * Math.Sqrt(2 * SpectralDensity(dw / 2 + i * dw) * dw) *
                                     (Math.Cos((dw / 2 + i * dw) * coordinates[(int)randomFieldDirection - 1] +
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

            return meanValue / (1 + randomCoefficient);

            //w = dw/2:dw:wu;
            ////omega = dw/2:dw:wu;
            ////Sff_target = feval(spectralDensityFunction,omega);
            //for (int i = 0; i < frequencyIntervals; i++)
            //    sffTarget[i] = SpectralDensity(dw / 2 + i * dw);

            //% 2D-Plate   
            //Lx = 1;

            //num_ele_x = 20;
            //num_ele_y = 20;
            //num_ele = num_ele_x*num_ele_y;
            //Le_x = Lx/num_ele_x;
            //for i = 1:num_ele_x
            //    midpoints(i) = Le_x/2+(i-1)*Le_x;
            //end

            //for i = 1:num_ele_y
            //    mid_nodes((i-1)*num_ele_x+1:num_ele_x*i) = midpoints;
            //end

            //n_phi = 10;
            //partitionMatrix = [1:n_phi];  
            //PhaseAnglesMatrix = (partitionMatrix - 1) * 2.*pi / n_phi + 2.*pi / (2 * n_phi);
            //phaseAnglesMatrixSize = size(PhaseAnglesMatrix);
            //FrequencyIntervals = size(w,2);
  
            //// Fast Monte Carlo
            //for frequencyCounter = 1:FrequencyIntervals
            //    for c = 1:n_phi
            //        r = 0;
            //        for xCounter = 1:size(mid_nodes,2)
            //            r = r+1;   
            //            v = mid_nodes(xCounter);
            //            sampleFunctionSet(frequencyCounter,c,r) = sqrt(2.)*Std * cos(w(frequencyCounter) * v + PhaseAnglesMatrix(1, c));
            //        end
            //    end
            //end

            //savefile = 'VRF.mat';
            //save(savefile, 'sampleFunctionSet','ku','dw','Sff','w','Std','T')
        }

        public double[] RandomVariables
        {
            get
            {
                
                return randomVariables; }
            set
            {
                phi = new double[frequencyIntervals];
                var Phi = new ContinuousUniformDistribution(0d, 2d * Math.PI);
                for (int i = 0; i < frequencyIntervals; i++)
                {
                    phi[i] = Phi.NextDouble();
                }
                //changedVariables = true;
                randomVariables = value;
            }
        }
    }
}
