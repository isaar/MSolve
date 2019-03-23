using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.IGA.Problems.SupportiveClasses
{
    public class BSPLines1D
    {
        public int Degree { get; private set; }

        public IVector KnotValueVector { get; private set; }

        public IVector ParametricCoordinates { get; private set; }

        public double[,] BSPLineValues { get; private set; }

        public double[,] BSPLineDerivativeValues { get; private set; }

		public double [,] BSPLineSecondDerivativeValues { get; private set; }

        public BSPLines1D(int degree, IVector knotValueVector, IVector parametricCoordinates)
        {
            if (degree <= 0)
            {
                throw new ArgumentException();
            }else if(knotValueVector == null){ throw new ArgumentNullException(); }
            this.Degree = degree;
            this.KnotValueVector = knotValueVector;
            this.ParametricCoordinates = parametricCoordinates;
        }

        public void calculateBSPLinesAndDerivatives()
        {
            int order = Degree + 1;
            int numberOfControlPoints = KnotValueVector.Length - order;
            int numberOfGaussPoints = ParametricCoordinates.Length;
            BSPLineValues = new double[numberOfControlPoints + Degree, numberOfGaussPoints];
            BSPLineDerivativeValues = new double[numberOfControlPoints + Degree, numberOfGaussPoints];
			BSPLineSecondDerivativeValues = new double[numberOfControlPoints + Degree, numberOfGaussPoints];

            for (int i = 0; i < numberOfGaussPoints; i++)
                for (int j = 0; j < numberOfControlPoints+Degree; j++)
	                if (KnotValueVector[j]<=ParametricCoordinates[i]&& ParametricCoordinates[i] <= KnotValueVector[j + 1])
		                BSPLineValues[j, i] = 1;
	                else
		                BSPLineValues[j, i] = 0;

            for (int i = 1; i <= Degree; i++)
            {
                for (int j = 0; j < numberOfControlPoints+Degree-i; j++)
                {
                    for (int k = 0; k < numberOfGaussPoints; k++)
                    {
                        double additive1 = 0;
                        double additive2 = 0;
                        double additiveDerivative1 = 0;
                        double additiveDerivative2 = 0;
	                    double additiveSecondDerivative1 = 0;
	                    double additiveSecondDerivative2 = 0;
						double denominator1 = KnotValueVector[j + i] - KnotValueVector[j];
                        double denominator2 = KnotValueVector[j + i + 1] - KnotValueVector[j + 1];
                        if (denominator1 != 0)
                        {
                            additive1 = (ParametricCoordinates[k] - KnotValueVector[j]) 
                                / denominator1 * BSPLineValues[j, k];
                            additiveDerivative1 = ((ParametricCoordinates[k] - KnotValueVector[j])
                                * BSPLineDerivativeValues[j, k] + BSPLineValues[j, k]) / denominator1;
	                        additiveSecondDerivative1 = Degree / denominator1 * BSPLineDerivativeValues[j, k];
                        }
                        if (denominator2 != 0)
                        {
                            additive2 = (KnotValueVector[j + i + 1] - ParametricCoordinates[k])
                                / denominator2 * BSPLineValues[j + 1, k];
                            additiveDerivative2 = ((KnotValueVector[j + i + 1] - ParametricCoordinates[k])
                                * BSPLineDerivativeValues[j + 1, k] - BSPLineValues[j + 1, k]) / denominator2;
	                        additiveSecondDerivative2 = -Degree / denominator2 * BSPLineDerivativeValues[j + 1, k];
                        }
                        BSPLineValues[j, k] = additive1 + additive2;
                        BSPLineDerivativeValues[j, k] = additiveDerivative1 + additiveDerivative2;
	                    BSPLineSecondDerivativeValues[j, k] = additiveSecondDerivative1 + additiveSecondDerivative2;
                    }
                }
            }
        }
    }
}
