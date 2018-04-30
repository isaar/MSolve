using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.Matrices.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISSAR.MSolve.IGAPreProcessor.Integration
{
    class BSPLines1D
    {
        public IMatrix2D<double> BSPLValues { get; }

        public IMatrix2D<double> BSPLDerivativeValues { get; }

        public BSPLines1D(int degree, IVector<double> knotValueVector, IVector<double> parametricCoordinates)
        {
            if (degree <= 0)
            {
                throw new IndexOutOfRangeException("The degree must be greater or equal to zero");
            } else if (knotValueVector==null){
                throw new ArgumentNullException("Knot Value Vector is null");
            }

            int order = degree + 1;

            int numberOfCP = knotValueVector.Length - order;

            int numberOfGP = parametricCoordinates.Length;

            BSPLValues = new Matrix2D<double>(numberOfCP+degree,numberOfGP);

            BSPLDerivativeValues = new Matrix2D<double>(numberOfCP + degree, numberOfGP);

            // BSPLine Basis Functions. Degree=0.
            for (int i = 0; i < numberOfGP; i++)
            {
                for (int j = 0; j < numberOfCP; j++)
                {
                    if (knotValueVector[j]<=parametricCoordinates[i]&&parametricCoordinates[i]<knotValueVector[j+1])
                    {
                        BSPLValues[j, i] = 1;
                    }else
                    {
                        BSPLValues[j, i] = 0;
                    }
                }
            }

            // BSPLineBasis Basis Functions. Degree=1:p
            for (int i = 0; i < degree; i++)
            {
                for (int j = 0; j < numberOfCP+degree-i; j++)
                {
                    for (int k = 0; k < numberOfGP; k++)
                    {
                        double additive1 = 0;
                        double additive2 = 0;
                        double additiveDerivative1 = 0;
                        double additiveDerivative2 = 0;

                        double denominator1 = knotValueVector[j + i] - knotValueVector[j];
                        double denominator2 = knotValueVector[j + i + 1] - knotValueVector[j + 1];

                        if (denominator1 != 0)
                        {
                            additive1 = (parametricCoordinates[k] - knotValueVector[j]) / denominator1 * BSPLValues[j, k];
                            additiveDerivative1 = ((parametricCoordinates[k]-knotValueVector[j])*BSPLDerivativeValues[j,k]
                                +BSPLValues[j,k]) / denominator1;
                        }
                        if (denominator2 != 0)
                        {
                            additive2 = (knotValueVector[j + i + 1] - parametricCoordinates[k]) / denominator2 * BSPLValues[j + 1, k];
                            additiveDerivative2 = ((knotValueVector[j + i + 1] - parametricCoordinates[k]) * BSPLDerivativeValues[j + 1, k] - BSPLValues[j + 1, k]) / denominator2;
                        }
                        BSPLValues[j, k] = additive1 + additive2;
                        BSPLDerivativeValues[j, k] = additiveDerivative1 + additiveDerivative2;

                    }
                }
            }








        }

        
    }
}
