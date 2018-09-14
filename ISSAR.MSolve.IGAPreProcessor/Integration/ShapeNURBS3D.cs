using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.Matrices.Interfaces;
using ISSAR.MSolve.IGAPreProcessor.Elements;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISSAR.MSolve.IGAPreProcessor.Integration
{
    class ShapeNURBS3D
    {
        public IMatrix2D<double> nurbsValues { get; }
        public IMatrix2D<double> nurbsDerivativeValuesKsi { get; }
        public IMatrix2D<double> nurbsDerivativeValuesHeta { get; }
        public IMatrix2D<double> nurbsDerivativeValuesZeta { get; }

        public ShapeNURBS3D(NURBSElement3D element, IGAModel model, IVector<double> parametricCoordinateKsi, IVector<double> parametricCoordinateHeta, IVector<double> parametricCoordinateZeta, IList<ControlPoint> controlPoints)
        {
            int numberOfCPHeta = model.KnotValueVectorHeta.Length - model.DegreeHeta - 1;
            int numberOfCPZeta = model.KnotValueVectorZeta.Length - model.DegreeZeta - 1;

            BSPLines1D bsplinesKsi = new BSPLines1D(model.DegreeKsi, model.KnotValueVectorKsi, parametricCoordinateKsi);
            BSPLines1D bsplinesHeta = new BSPLines1D(model.DegreeHeta, model.KnotValueVectorHeta, parametricCoordinateHeta);
            BSPLines1D bsplinesZeta = new BSPLines1D(model.DegreeZeta, model.KnotValueVectorZeta, parametricCoordinateZeta);

            nurbsValues = new Matrix2D<double>((model.DegreeKsi + 1) * (model.DegreeHeta + 1)* (model.DegreeZeta + 1), element.gaussPoints.Count);
            nurbsDerivativeValuesKsi = new Matrix2D<double>((model.DegreeKsi + 1) * (model.DegreeHeta + 1)* (model.DegreeZeta + 1), element.gaussPoints.Count);
            nurbsDerivativeValuesHeta = new Matrix2D<double>((model.DegreeKsi + 1) * (model.DegreeHeta + 1)* (model.DegreeZeta + 1), element.gaussPoints.Count);
            nurbsDerivativeValuesZeta = new Matrix2D<double>((model.DegreeKsi + 1) * (model.DegreeHeta + 1) * (model.DegreeZeta + 1), element.gaussPoints.Count);

            int supportKsi = model.DegreeKsi + 1;
            int supportHeta = model.DegreeHeta + 1;
            int supportZeta = model.DegreeZeta + 1;
            int numberOfElementGaussPoints = (model.DegreeKsi + 1) * (model.DegreeHeta + 1) * (model.DegreeZeta + 1);


            for (int i = 0; i < supportKsi; i++)
            {
                for (int j = 0; j < supportHeta; j++)
                {
                    for (int k = 0; k < supportZeta; k++)
                    {
                        double sumKsiHetaZeta = 0;
                        double sumdKsiHetaZeta = 0;
                        double sumKsidHetaZeta = 0;
                        double sumKsiHetadZeta = 0;

                        for (int m = 0; m < numberOfElementGaussPoints; m++)
                        {
                            int indexKsi = element.connectivity[m] / (numberOfCPHeta * numberOfCPZeta);
                            int indexHeta = element.connectivity[m] % (numberOfCPHeta * numberOfCPZeta) / numberOfCPZeta;
                            int indexZeta = element.connectivity[m] % (numberOfCPHeta * numberOfCPZeta) % numberOfCPZeta;

                            sumKsiHetaZeta += bsplinesKsi.BSPLValues[indexKsi, i] *
                                bsplinesHeta.BSPLValues[indexHeta, j] *
                                bsplinesZeta.BSPLValues[indexZeta, k] *
                                controlPoints[m].Weight;

                            sumdKsiHetaZeta += bsplinesKsi.BSPLDerivativeValues[indexKsi,i]*
                                bsplinesHeta.BSPLValues[indexHeta, j] *
                                bsplinesZeta.BSPLValues[indexZeta, k] *
                                controlPoints[m].Weight;

                            sumKsidHetaZeta += bsplinesKsi.BSPLValues[indexKsi, i] *
                                bsplinesHeta.BSPLDerivativeValues[indexHeta,j] *
                                bsplinesZeta.BSPLValues[indexZeta, k] *
                                controlPoints[m].Weight;

                            sumKsiHetadZeta += bsplinesKsi.BSPLValues[indexKsi, i] *
                                bsplinesHeta.BSPLValues[indexHeta, j] *
                                bsplinesZeta.BSPLDerivativeValues[indexZeta,k] *
                                controlPoints[m].Weight;
                        }
                        for (int m = 0; m < numberOfElementGaussPoints; m++)
                        {
                            int indexKsi = element.connectivity[m] / (numberOfCPHeta * numberOfCPZeta);
                            int indexHeta = element.connectivity[m] % (numberOfCPHeta * numberOfCPZeta) / numberOfCPZeta;
                            int indexZeta = element.connectivity[m] % (numberOfCPHeta * numberOfCPZeta) % numberOfCPZeta;

                            nurbsValues[m, i * supportHeta * supportZeta + j * supportZeta + k] =
                                bsplinesKsi.BSPLValues[indexKsi, i] *
                                bsplinesHeta.BSPLValues[indexHeta, j] *
                                bsplinesZeta.BSPLValues[indexZeta, k] *
                                controlPoints[m].Weight / sumKsiHetaZeta;

                            nurbsDerivativeValuesKsi[m, i * supportHeta * supportZeta + j * supportZeta + k] =
                                (bsplinesKsi.BSPLDerivativeValues[indexKsi, i] * sumKsiHetaZeta -
                                bsplinesKsi.BSPLValues[indexKsi, i] * sumdKsiHetaZeta) *
                                bsplinesHeta.BSPLValues[indexHeta, j] *
                                bsplinesZeta.BSPLValues[indexZeta, k] *
                                controlPoints[m].Weight / Math.Pow(sumKsiHetaZeta, 2);

                            nurbsDerivativeValuesHeta[m, i * supportHeta * supportZeta + j * supportZeta + k] =
                                bsplinesKsi.BSPLValues[indexKsi, i] *
                                (bsplinesHeta.BSPLDerivativeValues[indexHeta, j] * sumKsiHetaZeta -
                                bsplinesHeta.BSPLValues[indexHeta, j] * sumKsidHetaZeta) *
                                bsplinesZeta.BSPLValues[indexZeta, k] *
                                controlPoints[m].Weight / Math.Pow(sumKsiHetaZeta, 2);

                            nurbsDerivativeValuesZeta[m, i * supportHeta * supportZeta + j * supportZeta + k] =
                                bsplinesKsi.BSPLValues[indexKsi, i] *
                                bsplinesHeta.BSPLValues[indexHeta, j] *
                                (bsplinesZeta.BSPLDerivativeValues[indexZeta, k] * sumKsiHetaZeta -
                                bsplinesZeta.BSPLValues[indexZeta, k] * sumKsiHetadZeta) *
                                controlPoints[m].Weight / Math.Pow(sumKsiHetaZeta, 2);
                        }

                    }
                }
            }
        }
    }
}
