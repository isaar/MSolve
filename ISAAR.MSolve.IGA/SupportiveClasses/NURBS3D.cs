using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Problems.Structural.Elements;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.IGA.Problems.SupportiveClasses
{
    public class NURBS3D
    {
        public Matrix NurbsValues { get; private set; }
        public Matrix NurbsDerivativeValuesKsi { get; private set; }
        public Matrix NurbsDerivativeValuesHeta { get; private set; }
        public Matrix NurbsDerivativeValuesZeta { get; private set; }

        public Matrix NurbsSecondDerivativeValueKsi { get; private set; }
        public Matrix NurbsSecondDerivativeValueHeta { get; private set; }
        public Matrix NurbsSecondDerivativeValueZeta { get; private set; }
        public Matrix NurbsSecondDerivativeValueKsiHeta { get; private set; }
        public Matrix NurbsSecondDerivativeValueKsiZeta { get; private set; }
        public Matrix NurbsSecondDerivativeValueHetaZeta { get; private set; }

        public NURBS3D(int numberOfControlPointsKsi, int numberOfControlPointsHeta, int numberOfControlPointsZeta,
            int degreeKsi, int degreeHeta, int degreeZeta, Vector knotValueVectorKsi,
            Vector knotValueVectorHeta, Vector knotValueVectorZeta, ControlPoint[] controlPoints,
            NaturalPoint collocationPoint)
        {
            BSPLines1D bsplinesKsi =
                new BSPLines1D(degreeKsi, knotValueVectorKsi, Vector.CreateFromArray(new double[] { collocationPoint.Xi }));
            BSPLines1D bsplinesHeta = new BSPLines1D(degreeHeta, knotValueVectorHeta,
                Vector.CreateFromArray(new double[] { collocationPoint.Eta }));
            BSPLines1D bsplinesZeta = new BSPLines1D(degreeZeta, knotValueVectorZeta,
                Vector.CreateFromArray(new double[] { collocationPoint.Zeta }));
            bsplinesKsi.calculateBSPLinesAndDerivatives();
            bsplinesHeta.calculateBSPLinesAndDerivatives();
            bsplinesZeta.calculateBSPLinesAndDerivatives();

            int supportKsi = degreeKsi + 1;
            int supportHeta = degreeHeta + 1;
            int supportZeta = degreeZeta + 1;
            int numberOfGPKsi = 1;
            int numberOfGPHeta = 1;
            int numberOfGPZeta = 1;
            int numberOfElementControlPoints = supportKsi * supportHeta * supportZeta;

            NurbsValues = Matrix.CreateZero(numberOfElementControlPoints, 1);
            NurbsDerivativeValuesKsi = Matrix.CreateZero(numberOfElementControlPoints, 1);
            NurbsDerivativeValuesHeta = Matrix.CreateZero(numberOfElementControlPoints, 1);
            NurbsDerivativeValuesZeta = Matrix.CreateZero(numberOfElementControlPoints, 1);
            NurbsSecondDerivativeValueKsi = Matrix.CreateZero(numberOfElementControlPoints, 1);
            NurbsSecondDerivativeValueHeta = Matrix.CreateZero(numberOfElementControlPoints, 1);
            NurbsSecondDerivativeValueZeta = Matrix.CreateZero(numberOfElementControlPoints, 1);
            NurbsSecondDerivativeValueKsiHeta = Matrix.CreateZero(numberOfElementControlPoints, 1);
            NurbsSecondDerivativeValueKsiZeta = Matrix.CreateZero(numberOfElementControlPoints, 1);
            NurbsSecondDerivativeValueHetaZeta = Matrix.CreateZero(numberOfElementControlPoints, 1);

            for (int i = 0; i < numberOfGPKsi; i++)
            {
                for (int j = 0; j < numberOfGPHeta; j++)
                {
                    for (int k = 0; k < numberOfGPZeta; k++)
                    {
                        double sumKsiHetaZeta = 0;
                        double sumdKsiHetaZeta = 0;
                        double sumKsidHetaZeta = 0;
                        double sumKsiHetadZeta = 0;

                        double sumdKsidKsi = 0;
                        double sumdHetadHeta = 0;
                        double sumdZetadZeta = 0;
                        double sumdKsidHeta = 0;
                        double sumdKsidZeta = 0;
                        double sumdHetadZeta = 0;

                        for (int m = 0; m < numberOfElementControlPoints; m++)
                        {
                            int indexKsi = controlPoints[m].ID /
                                           (numberOfControlPointsHeta *
                                            numberOfControlPointsZeta);
                            int indexHeta = controlPoints[m].ID %
                                            (numberOfControlPointsHeta *
                                             numberOfControlPointsZeta) /
                                            numberOfControlPointsZeta;
                            int indexZeta = controlPoints[m].ID %
                                            (numberOfControlPointsHeta *
                                             numberOfControlPointsZeta) %
                                            numberOfControlPointsZeta;

                            sumKsiHetaZeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
                                              bsplinesHeta.BSPLineValues[indexHeta, j] *
                                              bsplinesZeta.BSPLineValues[indexZeta, k] *
                                              controlPoints[m].WeightFactor;

                            sumdKsiHetaZeta += bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
                                               bsplinesHeta.BSPLineValues[indexHeta, j] *
                                               bsplinesZeta.BSPLineValues[indexZeta, k] *
                                               controlPoints[m].WeightFactor;

                            sumKsidHetaZeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
                                               bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
                                               bsplinesZeta.BSPLineValues[indexZeta, k] *
                                               controlPoints[m].WeightFactor;

                            sumKsiHetadZeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
                                               bsplinesHeta.BSPLineValues[indexHeta, j] *
                                               bsplinesZeta.BSPLineDerivativeValues[indexZeta, k] *
                                               controlPoints[m].WeightFactor;

                            sumdKsidKsi += bsplinesKsi.BSPLineSecondDerivativeValues[indexKsi, i] *
                                           bsplinesHeta.BSPLineValues[indexHeta, j] *
                                           bsplinesZeta.BSPLineValues[indexZeta, k] *
                                           controlPoints[m].WeightFactor;

                            sumdHetadHeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
                                             bsplinesHeta.BSPLineSecondDerivativeValues[indexHeta, j] *
                                             bsplinesZeta.BSPLineValues[indexZeta, k] *
                                             controlPoints[m].WeightFactor;

                            sumdZetadZeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
                                             bsplinesHeta.BSPLineValues[indexHeta, j] *
                                             bsplinesZeta.BSPLineSecondDerivativeValues[indexZeta, k] *
                                             controlPoints[m].WeightFactor;

                            sumdKsidHeta += bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
                                            bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
                                            bsplinesZeta.BSPLineValues[indexZeta, k] *
                                            controlPoints[m].WeightFactor;

                            sumdKsidZeta += bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
                                            bsplinesHeta.BSPLineValues[indexHeta, j] *
                                            bsplinesZeta.BSPLineDerivativeValues[indexZeta, k] *
                                            controlPoints[m].WeightFactor;

                            sumdHetadZeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
                                             bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
                                             bsplinesZeta.BSPLineDerivativeValues[indexZeta, k] *
                                             controlPoints[m].WeightFactor;
                        }

                        for (int m = 0; m < numberOfElementControlPoints; m++)
                        {
                            int indexKsi = controlPoints[m].ID /
                                           (numberOfControlPointsHeta *
                                            numberOfControlPointsZeta);
                            int indexHeta = controlPoints[m].ID %
                                            (numberOfControlPointsHeta *
                                             numberOfControlPointsZeta) /
                                            numberOfControlPointsZeta;
                            int indexZeta = controlPoints[m].ID %
                                            (numberOfControlPointsHeta *
                                             numberOfControlPointsZeta) %
                                            numberOfControlPointsZeta;

                            NurbsValues[m, i * supportHeta * supportZeta + j * supportZeta + k] =
                                bsplinesKsi.BSPLineValues[indexKsi, i] *
                                bsplinesHeta.BSPLineValues[indexHeta, j] *
                                bsplinesZeta.BSPLineValues[indexZeta, k] *
                                controlPoints[m].WeightFactor / sumKsiHetaZeta;

                            NurbsDerivativeValuesKsi[m, i * supportHeta * supportZeta + j * supportZeta + k] =
                                (bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * sumKsiHetaZeta -
                                 bsplinesKsi.BSPLineValues[indexKsi, i] * sumdKsiHetaZeta) *
                                bsplinesHeta.BSPLineValues[indexHeta, j] *
                                bsplinesZeta.BSPLineValues[indexZeta, k] *
                                controlPoints[m].WeightFactor / Math.Pow(sumKsiHetaZeta, 2);

                            NurbsDerivativeValuesHeta[m, i * supportHeta * supportZeta + j * supportZeta + k] =
                                bsplinesKsi.BSPLineValues[indexKsi, i] *
                                (bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] * sumKsiHetaZeta -
                                 bsplinesHeta.BSPLineValues[indexHeta, j] * sumKsidHetaZeta) *
                                bsplinesZeta.BSPLineValues[indexZeta, k] *
                                controlPoints[m].WeightFactor / Math.Pow(sumKsiHetaZeta, 2);

                            NurbsDerivativeValuesZeta[m, i * supportHeta * supportZeta + j * supportZeta + k] =
                                bsplinesKsi.BSPLineValues[indexKsi, i] *
                                bsplinesHeta.BSPLineValues[indexHeta, j] *
                                (bsplinesZeta.BSPLineDerivativeValues[indexZeta, k] * sumKsiHetaZeta -
                                 bsplinesZeta.BSPLineValues[indexZeta, k] * sumKsiHetadZeta) *
                                controlPoints[m].WeightFactor / Math.Pow(sumKsiHetaZeta, 2);

                            NurbsSecondDerivativeValueKsi[m, i * supportHeta * supportZeta + j * supportZeta + k] =
                                (bsplinesKsi.BSPLineSecondDerivativeValues[indexKsi, i] / sumKsiHetaZeta -
                                 2 * bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * sumdKsiHetaZeta /
                                 Math.Pow(sumKsiHetaZeta, 2) -
                                 bsplinesKsi.BSPLineValues[indexKsi, i] * sumdKsidKsi / Math.Pow(sumKsiHetaZeta, 2) +
                                 2 * bsplinesKsi.BSPLineValues[indexKsi, i] * Math.Pow(sumdKsiHetaZeta, 2) /
                                 Math.Pow(sumKsiHetaZeta, 3)) *
                                bsplinesHeta.BSPLineValues[indexHeta, j] *
                                bsplinesZeta.BSPLineValues[indexZeta, k] *
                                controlPoints[m].WeightFactor;

                            NurbsSecondDerivativeValueHeta[m, i * supportHeta * supportZeta + j * supportZeta + k] =
                                bsplinesKsi.BSPLineValues[indexKsi, i] *
                                (bsplinesHeta.BSPLineSecondDerivativeValues[indexHeta, j] / sumKsiHetaZeta -
                                 2 * bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] * sumKsidHetaZeta / Math.Pow(sumKsiHetaZeta, 2) -
                                 bsplinesHeta.BSPLineValues[indexHeta, j] * sumdHetadHeta / Math.Pow(sumKsiHetaZeta, 2) +
                                 2 * bsplinesHeta.BSPLineValues[indexHeta, j] * Math.Pow(sumKsidHetaZeta, 2) / Math.Pow(sumKsiHetaZeta, 3)) *
                                bsplinesZeta.BSPLineValues[indexZeta, k] *
                                controlPoints[m].WeightFactor;

                            NurbsSecondDerivativeValueZeta[m, i * supportHeta * supportZeta + j * supportZeta + k] =
                                bsplinesKsi.BSPLineValues[indexKsi, i] *
                                bsplinesHeta.BSPLineValues[indexHeta, j] *
                                (bsplinesZeta.BSPLineSecondDerivativeValues[indexZeta, k] / sumKsiHetaZeta -
                                 2 * bsplinesZeta.BSPLineDerivativeValues[indexZeta, k] * sumKsiHetadZeta /
                                 Math.Pow(sumKsiHetaZeta, 2) -
                                 bsplinesZeta.BSPLineValues[indexZeta, k] * sumdZetadZeta /
                                 Math.Pow(sumKsiHetaZeta, 2) +
                                 2 * bsplinesZeta.BSPLineValues[indexZeta, k] * Math.Pow(sumKsiHetadZeta, 2) /
                                 Math.Pow(sumKsiHetaZeta, 3)) *
                                controlPoints[m].WeightFactor;

                            NurbsSecondDerivativeValueKsiHeta[m, i * supportHeta * supportZeta + j * supportZeta + k] =
                                (bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
                                 bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] / sumKsiHetaZeta -
                                 bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
                                 bsplinesHeta.BSPLineValues[indexHeta, j] *
                                 sumKsidHetaZeta / Math.Pow(sumKsiHetaZeta, 2) -
                                 bsplinesKsi.BSPLineValues[indexKsi, i] *
                                 bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
                                 sumdKsiHetaZeta / Math.Pow(sumKsiHetaZeta, 2) -
                                 bsplinesKsi.BSPLineValues[indexKsi, i] * bsplinesHeta.BSPLineValues[indexHeta, j] *
                                 sumdKsidHeta / Math.Pow(sumKsiHetaZeta, 2) +
                                 2 * bsplinesKsi.BSPLineValues[indexKsi, i] * bsplinesHeta.BSPLineValues[indexHeta, j] *
                                 sumdKsiHetaZeta * sumKsidHetaZeta / Math.Pow(sumKsiHetaZeta, 3)) *
                                bsplinesZeta.BSPLineValues[indexZeta, k] *
                                controlPoints[m].WeightFactor;

                            NurbsSecondDerivativeValueKsiZeta[m, i * supportHeta * supportZeta + j * supportZeta + k] =
                                (bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
                                 bsplinesZeta.BSPLineDerivativeValues[indexZeta, k] / sumKsiHetaZeta -
                                 bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
                                 bsplinesZeta.BSPLineValues[indexZeta, k] *
                                 sumKsiHetadZeta / Math.Pow(sumKsiHetaZeta, 2) -
                                 bsplinesKsi.BSPLineValues[indexKsi, i] *
                                 bsplinesZeta.BSPLineDerivativeValues[indexZeta, k] *
                                 sumdKsiHetaZeta / Math.Pow(sumKsiHetaZeta, 2) -
                                 bsplinesKsi.BSPLineValues[indexKsi, i] * bsplinesZeta.BSPLineValues[indexZeta, k] *
                                 sumdKsidZeta / Math.Pow(sumKsiHetaZeta, 2) +
                                 2 * bsplinesKsi.BSPLineValues[indexKsi, i] * bsplinesZeta.BSPLineValues[indexZeta, k] *
                                 sumdKsiHetaZeta * sumKsiHetadZeta / Math.Pow(sumKsiHetaZeta, 3)) *
                                bsplinesHeta.BSPLineValues[indexHeta, j] *
                                controlPoints[m].WeightFactor;

                            NurbsSecondDerivativeValueHetaZeta[m, i * supportHeta * supportZeta + j * supportZeta + k] =
                                bsplinesKsi.BSPLineValues[indexKsi, i] *
                                (bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
                                 bsplinesZeta.BSPLineDerivativeValues[indexZeta, k] / sumKsiHetaZeta -
                                 bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
                                 bsplinesZeta.BSPLineValues[indexZeta, k] *
                                 sumKsiHetadZeta / Math.Pow(sumKsiHetaZeta, 2) -
                                 bsplinesHeta.BSPLineValues[indexHeta, j] *
                                 bsplinesZeta.BSPLineDerivativeValues[indexZeta, k] *
                                 sumKsidHetaZeta / Math.Pow(sumKsiHetaZeta, 2) -
                                 bsplinesHeta.BSPLineValues[indexHeta, j] * bsplinesZeta.BSPLineValues[indexZeta, k] *
                                 sumdHetadZeta / Math.Pow(sumKsiHetaZeta, 2) +
                                 2 * bsplinesHeta.BSPLineValues[indexHeta, j] *
                                 bsplinesZeta.BSPLineValues[indexZeta, k] *
                                 sumKsidHetaZeta * sumKsiHetadZeta / Math.Pow(sumKsiHetaZeta, 3)) *
                                controlPoints[m].WeightFactor;


                        }
                    }
                }
            }
        }

        public NURBS3D(int numberOfControlPointsKsi, int numberOfControlPointsHeta, int numberOfControlPointsZeta,
            int degreeKsi, int degreeHeta, int degreeZeta, Vector knotValueVectorKsi,
            Vector knotValueVectorHeta, Vector knotValueVectorZeta, ControlPoint[] controlPoints,
            GaussLegendrePoint3D[] gaussPoints)
        {
            var numberOfGaussPoints = gaussPoints.Length;
            var parametricGaussPointKsi = Vector.CreateZero(degreeKsi + 1);
            for (int i = 0; i < degreeKsi + 1; i++)
                parametricGaussPointKsi[i] = gaussPoints[i * (degreeZeta + 1) * (degreeHeta + 1)].Ksi;

            var parametricGaussPointHeta = Vector.CreateZero(degreeHeta + 1);
            for (int i = 0; i < degreeHeta + 1; i++)
                parametricGaussPointHeta[i] = gaussPoints[i * (degreeZeta + 1)].Heta;

            var parametricGaussPointZeta = Vector.CreateZero(degreeZeta + 1);
            for (int i = 0; i < degreeZeta + 1; i++)
                parametricGaussPointZeta[i] = gaussPoints[i].Zeta;

            BSPLines1D bsplinesKsi =
                new BSPLines1D(degreeKsi, knotValueVectorKsi, parametricGaussPointKsi);
            BSPLines1D bsplinesHeta = new BSPLines1D(degreeHeta, knotValueVectorHeta,
                parametricGaussPointHeta);
            BSPLines1D bsplinesZeta = new BSPLines1D(degreeZeta, knotValueVectorZeta,
                parametricGaussPointZeta);
            bsplinesKsi.calculateBSPLinesAndDerivatives();
            bsplinesHeta.calculateBSPLinesAndDerivatives();
            bsplinesZeta.calculateBSPLinesAndDerivatives();

            int supportKsi = degreeKsi + 1;
            int supportHeta = degreeHeta + 1;
            int supportZeta = degreeZeta + 1;
            int numberOfElementControlPoints = supportKsi * supportHeta * supportZeta;

            NurbsValues = Matrix.CreateZero(numberOfElementControlPoints, numberOfGaussPoints);
            NurbsDerivativeValuesKsi = Matrix.CreateZero(numberOfElementControlPoints, numberOfGaussPoints);
            NurbsDerivativeValuesHeta = Matrix.CreateZero(numberOfElementControlPoints, numberOfGaussPoints);
            NurbsDerivativeValuesZeta = Matrix.CreateZero(numberOfElementControlPoints, numberOfGaussPoints);

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

                        for (int m = 0; m < numberOfElementControlPoints; m++)
                        {
                            int indexKsi = controlPoints[m].ID /
                                           (numberOfControlPointsHeta *
                                            numberOfControlPointsZeta);
                            int indexHeta = controlPoints[m].ID %
                                            (numberOfControlPointsHeta *
                                             numberOfControlPointsZeta) /
                                            numberOfControlPointsZeta;
                            int indexZeta = controlPoints[m].ID %
                                            (numberOfControlPointsHeta *
                                             numberOfControlPointsZeta) %
                                            numberOfControlPointsZeta;

                            sumKsiHetaZeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
                                              bsplinesHeta.BSPLineValues[indexHeta, j] *
                                              bsplinesZeta.BSPLineValues[indexZeta, k] *
                                              controlPoints[m].WeightFactor;

                            sumdKsiHetaZeta += bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
                                               bsplinesHeta.BSPLineValues[indexHeta, j] *
                                               bsplinesZeta.BSPLineValues[indexZeta, k] *
                                               controlPoints[m].WeightFactor;

                            sumKsidHetaZeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
                                               bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
                                               bsplinesZeta.BSPLineValues[indexZeta, k] *
                                               controlPoints[m].WeightFactor;

                            sumKsiHetadZeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
                                               bsplinesHeta.BSPLineValues[indexHeta, j] *
                                               bsplinesZeta.BSPLineDerivativeValues[indexZeta, k] *
                                               controlPoints[m].WeightFactor;
                        }

                        for (int m = 0; m < numberOfElementControlPoints; m++)
                        {
                            int indexKsi = controlPoints[m].ID /
                                           (numberOfControlPointsHeta *
                                            numberOfControlPointsZeta);
                            int indexHeta = controlPoints[m].ID %
                                            (numberOfControlPointsHeta *
                                             numberOfControlPointsZeta) /
                                            numberOfControlPointsZeta;
                            int indexZeta = controlPoints[m].ID %
                                            (numberOfControlPointsHeta *
                                             numberOfControlPointsZeta) %
                                            numberOfControlPointsZeta;

                            NurbsValues[m, i * supportHeta * supportZeta + j * supportZeta + k] =
                                bsplinesKsi.BSPLineValues[indexKsi, i] *
                                bsplinesHeta.BSPLineValues[indexHeta, j] *
                                bsplinesZeta.BSPLineValues[indexZeta, k] *
                                controlPoints[m].WeightFactor / sumKsiHetaZeta;

                            NurbsDerivativeValuesKsi[m, i * supportHeta * supportZeta + j * supportZeta + k] =
                                (bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * sumKsiHetaZeta -
                                 bsplinesKsi.BSPLineValues[indexKsi, i] * sumdKsiHetaZeta) *
                                bsplinesHeta.BSPLineValues[indexHeta, j] *
                                bsplinesZeta.BSPLineValues[indexZeta, k] *
                                controlPoints[m].WeightFactor / Math.Pow(sumKsiHetaZeta, 2);

                            NurbsDerivativeValuesHeta[m, i * supportHeta * supportZeta + j * supportZeta + k] =
                                bsplinesKsi.BSPLineValues[indexKsi, i] *
                                (bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] * sumKsiHetaZeta -
                                 bsplinesHeta.BSPLineValues[indexHeta, j] * sumKsidHetaZeta) *
                                bsplinesZeta.BSPLineValues[indexZeta, k] *
                                controlPoints[m].WeightFactor / Math.Pow(sumKsiHetaZeta, 2);

                            NurbsDerivativeValuesZeta[m, i * supportHeta * supportZeta + j * supportZeta + k] =
                                bsplinesKsi.BSPLineValues[indexKsi, i] *
                                bsplinesHeta.BSPLineValues[indexHeta, j] *
                                (bsplinesZeta.BSPLineDerivativeValues[indexZeta, k] * sumKsiHetaZeta -
                                 bsplinesZeta.BSPLineValues[indexZeta, k] * sumKsiHetadZeta) *
                                controlPoints[m].WeightFactor / Math.Pow(sumKsiHetaZeta, 2);
                        }
                    }
                }
            }
        }

        public NURBS3D(Element element, IList<ControlPoint> controlPoints)
        {
            GaussQuadrature gauss = new GaussQuadrature();
            IList<GaussLegendrePoint3D> gaussPoints = gauss.CalculateElementGaussPoints(element.Patch.DegreeKsi, element.Patch.DegreeHeta, element.Patch.DegreeZeta, element.Knots);

            var parametricGaussPointKsi = Vector.CreateZero(element.Patch.DegreeKsi + 1);
            for (int i = 0; i < element.Patch.DegreeKsi + 1; i++)
            {
                parametricGaussPointKsi[i] = gaussPoints[i * (element.Patch.DegreeZeta + 1) * (element.Patch.DegreeHeta + 1)].Ksi;
            }

            var parametricGaussPointHeta = Vector.CreateZero(element.Patch.DegreeHeta + 1);
            for (int i = 0; i < element.Patch.DegreeHeta + 1; i++)
            {
                parametricGaussPointHeta[i] = gaussPoints[i * (element.Patch.DegreeZeta + 1)].Heta;
            }

            var parametricGaussPointZeta = Vector.CreateZero(element.Patch.DegreeZeta + 1);
            for (int i = 0; i < element.Patch.DegreeZeta + 1; i++)
            {
                parametricGaussPointZeta[i] = gaussPoints[i].Zeta;
            }

            BSPLines1D bsplinesKsi = new BSPLines1D(element.Patch.DegreeKsi, element.Patch.KnotValueVectorKsi, parametricGaussPointKsi);
            BSPLines1D bsplinesHeta = new BSPLines1D(element.Patch.DegreeHeta, element.Patch.KnotValueVectorHeta, parametricGaussPointHeta);
            BSPLines1D bsplinesZeta = new BSPLines1D(element.Patch.DegreeZeta, element.Patch.KnotValueVectorZeta, parametricGaussPointZeta);
            bsplinesKsi.calculateBSPLinesAndDerivatives();
            bsplinesHeta.calculateBSPLinesAndDerivatives();
            bsplinesZeta.calculateBSPLinesAndDerivatives();

            int supportKsi = element.Patch.DegreeKsi + 1;
            int supportHeta = element.Patch.DegreeHeta + 1;
            int supportZeta = element.Patch.DegreeZeta + 1;
            int numberOfElementControlPoints = supportKsi * supportHeta * supportZeta;

            NurbsValues = Matrix.CreateZero(numberOfElementControlPoints, gaussPoints.Count);
            NurbsDerivativeValuesKsi = Matrix.CreateZero(numberOfElementControlPoints, gaussPoints.Count);
            NurbsDerivativeValuesHeta = Matrix.CreateZero(numberOfElementControlPoints, gaussPoints.Count);
            NurbsDerivativeValuesZeta = Matrix.CreateZero(numberOfElementControlPoints, gaussPoints.Count);

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

                        for (int m = 0; m < numberOfElementControlPoints; m++)
                        {
                            int indexKsi = element.ControlPoints[m].ID / (element.Patch.NumberOfControlPointsHeta * element.Patch.NumberOfControlPointsZeta);
                            int indexHeta = element.ControlPoints[m].ID % (element.Patch.NumberOfControlPointsHeta * element.Patch.NumberOfControlPointsZeta) / element.Patch.NumberOfControlPointsZeta;
                            int indexZeta = element.ControlPoints[m].ID % (element.Patch.NumberOfControlPointsHeta * element.Patch.NumberOfControlPointsZeta) % element.Patch.NumberOfControlPointsZeta;

                            sumKsiHetaZeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
                                bsplinesHeta.BSPLineValues[indexHeta, j] *
                                bsplinesZeta.BSPLineValues[indexZeta, k] *
                                controlPoints[m].WeightFactor;

                            sumdKsiHetaZeta += bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
                                bsplinesHeta.BSPLineValues[indexHeta, j] *
                                bsplinesZeta.BSPLineValues[indexZeta, k] *
                                controlPoints[m].WeightFactor;

                            sumKsidHetaZeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
                                bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
                                bsplinesZeta.BSPLineValues[indexZeta, k] *
                                controlPoints[m].WeightFactor;

                            sumKsiHetadZeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
                                bsplinesHeta.BSPLineValues[indexHeta, j] *
                                bsplinesZeta.BSPLineDerivativeValues[indexZeta, k] *
                                controlPoints[m].WeightFactor;
                        }
                        for (int m = 0; m < numberOfElementControlPoints; m++)
                        {
                            int indexKsi = element.ControlPoints[m].ID / (element.Patch.NumberOfControlPointsHeta * element.Patch.NumberOfControlPointsZeta);
                            int indexHeta = element.ControlPoints[m].ID % (element.Patch.NumberOfControlPointsHeta * element.Patch.NumberOfControlPointsZeta) / element.Patch.NumberOfControlPointsZeta;
                            int indexZeta = element.ControlPoints[m].ID % (element.Patch.NumberOfControlPointsHeta * element.Patch.NumberOfControlPointsZeta) % element.Patch.NumberOfControlPointsZeta;

                            NurbsValues[m, i * supportHeta * supportZeta + j * supportZeta + k] =
                                bsplinesKsi.BSPLineValues[indexKsi, i] *
                                bsplinesHeta.BSPLineValues[indexHeta, j] *
                                bsplinesZeta.BSPLineValues[indexZeta, k] *
                                controlPoints[m].WeightFactor / sumKsiHetaZeta;

                            NurbsDerivativeValuesKsi[m, i * supportHeta * supportZeta + j * supportZeta + k] =
                                (bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * sumKsiHetaZeta -
                                bsplinesKsi.BSPLineValues[indexKsi, i] * sumdKsiHetaZeta) *
                                bsplinesHeta.BSPLineValues[indexHeta, j] *
                                bsplinesZeta.BSPLineValues[indexZeta, k] *
                                controlPoints[m].WeightFactor / Math.Pow(sumKsiHetaZeta, 2);

                            NurbsDerivativeValuesHeta[m, i * supportHeta * supportZeta + j * supportZeta + k] =
                                bsplinesKsi.BSPLineValues[indexKsi, i] *
                                (bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] * sumKsiHetaZeta -
                                bsplinesHeta.BSPLineValues[indexHeta, j] * sumKsidHetaZeta) *
                                bsplinesZeta.BSPLineValues[indexZeta, k] *
                                controlPoints[m].WeightFactor / Math.Pow(sumKsiHetaZeta, 2);

                            NurbsDerivativeValuesZeta[m, i * supportHeta * supportZeta + j * supportZeta + k] =
                                bsplinesKsi.BSPLineValues[indexKsi, i] *
                                bsplinesHeta.BSPLineValues[indexHeta, j] *
                                (bsplinesZeta.BSPLineDerivativeValues[indexZeta, k] * sumKsiHetaZeta -
                                bsplinesZeta.BSPLineValues[indexZeta, k] * sumKsiHetadZeta) *
                                controlPoints[m].WeightFactor / Math.Pow(sumKsiHetaZeta, 2);
                        }
                    }
                }
            }
        }

	    public NURBS3D(Element element, IList<ControlPoint> controlPoints, Vector parametricGaussPointKsi,
		    Vector parametricGaussPointHeta, Vector parametricGaussPointZeta)
	    {
		    var parametricPointsCount = parametricGaussPointKsi.Length * parametricGaussPointHeta.Length *
		                                parametricGaussPointZeta.Length;

			BSPLines1D bsplinesKsi = new BSPLines1D(element.Patch.DegreeKsi, element.Patch.KnotValueVectorKsi,
			    parametricGaussPointKsi);
		    BSPLines1D bsplinesHeta = new BSPLines1D(element.Patch.DegreeHeta, element.Patch.KnotValueVectorHeta,
			    parametricGaussPointHeta);
		    BSPLines1D bsplinesZeta = new BSPLines1D(element.Patch.DegreeZeta, element.Patch.KnotValueVectorZeta,
			    parametricGaussPointZeta);
		    bsplinesKsi.calculateBSPLinesAndDerivatives();
		    bsplinesHeta.calculateBSPLinesAndDerivatives();
		    bsplinesZeta.calculateBSPLinesAndDerivatives();

		    int supportKsi = parametricGaussPointKsi.Length;
		    int supportHeta = parametricGaussPointHeta.Length;
		    int supportZeta = parametricGaussPointZeta.Length;
		    int numberOfElementControlPoints = (element.Patch.DegreeKsi + 1) * (element.Patch.DegreeHeta + 1) *
		                                       (element.Patch.DegreeZeta + 1);

		    NurbsValues = Matrix.CreateZero(numberOfElementControlPoints, parametricPointsCount);
		    NurbsDerivativeValuesKsi = Matrix.CreateZero(numberOfElementControlPoints, parametricPointsCount);
		    NurbsDerivativeValuesHeta = Matrix.CreateZero(numberOfElementControlPoints, parametricPointsCount);
		    NurbsDerivativeValuesZeta = Matrix.CreateZero(numberOfElementControlPoints, parametricPointsCount);

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

					    for (int m = 0; m < numberOfElementControlPoints; m++)
					    {
						    int indexKsi = element.ControlPoints[m].ID /
						                   (element.Patch.NumberOfControlPointsHeta *
						                    element.Patch.NumberOfControlPointsZeta);
						    int indexHeta = element.ControlPoints[m].ID %
						                    (element.Patch.NumberOfControlPointsHeta *
						                     element.Patch.NumberOfControlPointsZeta) /
						                    element.Patch.NumberOfControlPointsZeta;
						    int indexZeta = element.ControlPoints[m].ID %
						                    (element.Patch.NumberOfControlPointsHeta *
						                     element.Patch.NumberOfControlPointsZeta) %
						                    element.Patch.NumberOfControlPointsZeta;

						    sumKsiHetaZeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
						                      bsplinesHeta.BSPLineValues[indexHeta, j] *
						                      bsplinesZeta.BSPLineValues[indexZeta, k] *
						                      controlPoints[m].WeightFactor;

						    sumdKsiHetaZeta += bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] *
						                       bsplinesHeta.BSPLineValues[indexHeta, j] *
						                       bsplinesZeta.BSPLineValues[indexZeta, k] *
						                       controlPoints[m].WeightFactor;

						    sumKsidHetaZeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
						                       bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] *
						                       bsplinesZeta.BSPLineValues[indexZeta, k] *
						                       controlPoints[m].WeightFactor;

						    sumKsiHetadZeta += bsplinesKsi.BSPLineValues[indexKsi, i] *
						                       bsplinesHeta.BSPLineValues[indexHeta, j] *
						                       bsplinesZeta.BSPLineDerivativeValues[indexZeta, k] *
						                       controlPoints[m].WeightFactor;
					    }

					    for (int m = 0; m < numberOfElementControlPoints; m++)
					    {
						    int indexKsi = element.ControlPoints[m].ID /
						                   (element.Patch.NumberOfControlPointsHeta *
						                    element.Patch.NumberOfControlPointsZeta);
						    int indexHeta = element.ControlPoints[m].ID %
						                    (element.Patch.NumberOfControlPointsHeta *
						                     element.Patch.NumberOfControlPointsZeta) /
						                    element.Patch.NumberOfControlPointsZeta;
						    int indexZeta = element.ControlPoints[m].ID %
						                    (element.Patch.NumberOfControlPointsHeta *
						                     element.Patch.NumberOfControlPointsZeta) %
						                    element.Patch.NumberOfControlPointsZeta;

						    NurbsValues[m, i * supportHeta * supportZeta + j * supportZeta + k] =
							    bsplinesKsi.BSPLineValues[indexKsi, i] *
							    bsplinesHeta.BSPLineValues[indexHeta, j] *
							    bsplinesZeta.BSPLineValues[indexZeta, k] *
							    controlPoints[m].WeightFactor / sumKsiHetaZeta;

						    NurbsDerivativeValuesKsi[m, i * supportHeta * supportZeta + j * supportZeta + k] =
							    (bsplinesKsi.BSPLineDerivativeValues[indexKsi, i] * sumKsiHetaZeta -
							     bsplinesKsi.BSPLineValues[indexKsi, i] * sumdKsiHetaZeta) *
							    bsplinesHeta.BSPLineValues[indexHeta, j] *
							    bsplinesZeta.BSPLineValues[indexZeta, k] *
							    controlPoints[m].WeightFactor / Math.Pow(sumKsiHetaZeta, 2);

						    NurbsDerivativeValuesHeta[m, i * supportHeta * supportZeta + j * supportZeta + k] =
							    bsplinesKsi.BSPLineValues[indexKsi, i] *
							    (bsplinesHeta.BSPLineDerivativeValues[indexHeta, j] * sumKsiHetaZeta -
							     bsplinesHeta.BSPLineValues[indexHeta, j] * sumKsidHetaZeta) *
							    bsplinesZeta.BSPLineValues[indexZeta, k] *
							    controlPoints[m].WeightFactor / Math.Pow(sumKsiHetaZeta, 2);

						    NurbsDerivativeValuesZeta[m, i * supportHeta * supportZeta + j * supportZeta + k] =
							    bsplinesKsi.BSPLineValues[indexKsi, i] *
							    bsplinesHeta.BSPLineValues[indexHeta, j] *
							    (bsplinesZeta.BSPLineDerivativeValues[indexZeta, k] * sumKsiHetaZeta -
							     bsplinesZeta.BSPLineValues[indexZeta, k] * sumKsiHetadZeta) *
							    controlPoints[m].WeightFactor / Math.Pow(sumKsiHetaZeta, 2);
					    }
				    }
			    }
		    }
	    }
    }
}
