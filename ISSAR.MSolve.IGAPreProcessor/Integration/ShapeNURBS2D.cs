using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.Matrices.Interfaces;
using ISSAR.MSolve.IGAPreProcessor.Elements;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISSAR.MSolve.IGAPreProcessor.Integration
{
    class ShapeNURBS2D
    {
        public IMatrix2D<double> nurbsValues { get; }
        public IMatrix2D<double> nurbsDerivativeValuesKsi { get; }
        public IMatrix2D<double> nurbsDerivativeValuesHeta { get; }

        public ShapeNURBS2D(NURBSElement2D element, IGAModel model, IVector<double> parametricCoordinateKsi, IVector<double> parametricCoordinateHeta, IList<ControlPoint> controlPoints)
        {
            int numberOfCPHeta = model.KnotValueVectorHeta.Length - model.DegreeHeta - 1;

            BSPLines1D bsplinesKsi = new BSPLines1D(model.DegreeKsi,model.KnotValueVectorKsi,parametricCoordinateKsi);
            BSPLines1D bsplinesHeta = new BSPLines1D(model.DegreeHeta, model.KnotValueVectorHeta, parametricCoordinateHeta);

            nurbsValues = new Matrix2D<double>((model.DegreeKsi+1)*(model.DegreeHeta+1),element.gaussPoints.Count);
            nurbsDerivativeValuesKsi = new Matrix2D<double>((model.DegreeKsi + 1) * (model.DegreeHeta + 1), element.gaussPoints.Count);
            nurbsDerivativeValuesHeta = new Matrix2D<double>((model.DegreeKsi + 1) * (model.DegreeHeta + 1), element.gaussPoints.Count);

            int supportKsi = model.DegreeKsi + 1;
            int supportHeta = model.DegreeHeta + 1;
            int numberOfElementGaussPoints = (model.DegreeKsi + 1) * (model.DegreeHeta + 1);

            for (int i = 0; i < supportKsi; i++)
            {
                for (int j = 0; j < supportHeta; j++)
                {
                    double sumKsiHeta = 0;
                    double sumdKsiHeta = 0;
                    double sumKsidHeta = 0;

                    for (int k = 0; k < numberOfElementGaussPoints; k++)
                    {
                        sumKsiHeta += bsplinesKsi.BSPLValues[element.connectivity[k] / numberOfCPHeta, i] *
                            bsplinesHeta.BSPLValues[element.connectivity[k] % numberOfCPHeta, j] *
                            controlPoints[k].Weight;

                        sumdKsiHeta = +bsplinesKsi.BSPLDerivativeValues[element.connectivity[k] / numberOfCPHeta, i] *
                            bsplinesHeta.BSPLDerivativeValues[element.connectivity[k] % numberOfCPHeta, j] *
                            controlPoints[k].Weight;

                        sumKsidHeta = +bsplinesKsi.BSPLValues[element.connectivity[k] / numberOfCPHeta, i] *
                            bsplinesHeta.BSPLDerivativeValues[element.connectivity[k] % numberOfCPHeta, j] *
                            controlPoints[k].Weight;
                    }
                    for (int k = 0; k < numberOfElementGaussPoints; k++)
                    {
                        nurbsValues[k, i * supportHeta + j] = bsplinesKsi.BSPLValues[element.connectivity[k] / numberOfCPHeta, i] *
                            bsplinesHeta.BSPLValues[element.connectivity[k] % numberOfCPHeta, j] *
                            controlPoints[k].Weight / sumKsiHeta;

                        nurbsDerivativeValuesKsi[k, i * supportHeta + j] = bsplinesHeta.BSPLValues[element.connectivity[k] % numberOfCPHeta, j] * controlPoints[k].Weight *
                            (bsplinesKsi.BSPLDerivativeValues[element.connectivity[k] / numberOfCPHeta, i] * sumKsiHeta -
                            bsplinesKsi.BSPLValues[element.connectivity[k] / numberOfCPHeta, i] * sumdKsiHeta) / Math.Pow(sumKsiHeta, 2);

                        nurbsDerivativeValuesHeta[k, i * supportHeta + j] = bsplinesKsi.BSPLValues[element.connectivity[k] / numberOfCPHeta, i] *
                            controlPoints[k].Weight *
                            (bsplinesHeta.BSPLDerivativeValues[element.connectivity[k] % numberOfCPHeta, j] * sumKsiHeta -
                            bsplinesHeta.BSPLValues[element.connectivity[k] % numberOfCPHeta, j] * sumKsidHeta) / Math.Pow(sumKsiHeta, 2);
                    }


                }
            }




        }

    }
}
