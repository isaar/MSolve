using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Materials.Interfaces;
using Troschuetz.Random;
using Troschuetz.Random.Distributions.Continuous;

namespace ISAAR.MSolve.Analyzers
{
    public enum RandomVariableDistributionType
    {
        Normal,
        Lognormal
    }

    public class RandomVariableTargetEvaluator:IStochasticMaterialCoefficientsProvider
    {
        private readonly double _standardDeviation;
        private readonly double _meanValue;
        private double _randomVariable;
        private readonly RandomVariableDistributionType _distributionType;


        public RandomVariableTargetEvaluator(double meanValue, double standardDeviation, RandomVariableDistributionType distributionType)
        {
            this._meanValue = meanValue;
            this._standardDeviation = standardDeviation;
            _distributionType = distributionType;
        }

        public double[] RandomVariables
        {
            get { return new double[1];}
            set { CalculateRandomVariable(); }
        }

        private void CalculateRandomVariable()
        {
            switch (_distributionType)
            {
                case RandomVariableDistributionType.Normal:
                    var normalDistribution = new NormalDistribution(_meanValue, _standardDeviation);
                    _randomVariable = normalDistribution.NextDouble();
                    break;
                case RandomVariableDistributionType.Lognormal:
                    var lognormal = new LognormalDistribution(_meanValue, _standardDeviation*_meanValue);
                    _randomVariable = lognormal.NextDouble();
                    break;
            }
        }

        public double GetCoefficient(double meanValue, double[] coordinates)
        {
            return 1/_randomVariable;
        }
    }
}
