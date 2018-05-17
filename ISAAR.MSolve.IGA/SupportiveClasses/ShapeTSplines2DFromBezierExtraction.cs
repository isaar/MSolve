using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Problems.SupportiveClasses;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.IGA.SupportiveClasses
{
    public class ShapeTSplines2DFromBezierExtraction
    {
	    public IMatrix2D TSplineValues { get; private set; }
	    public IMatrix2D TSplineDerivativeValuesKsi { get; private set; }
	    public IMatrix2D TSplineDerivativeValuesHeta { get; private set; }
		public IMatrix2D TSplineSecondDerivativesValueKsi { get; private set; }
	    public IMatrix2D TSplineSecondDerivativesValueHeta { get; private set; }
	    public IMatrix2D TSplineSecondDerivativesValueKsiHeta { get; private set; }

		public ShapeTSplines2DFromBezierExtraction(TSplineElement2D element, IList<ControlPoint> controlPoints)
	    {
			GaussQuadrature gauss = new GaussQuadrature();
		    IList<GaussLegendrePoint3D> gaussPoints = gauss.CalculateElementGaussPoints(element.DegreeKsi, element.DegreeHeta, 
			    new List<Knot>
			    {
				    new Knot(){ID=0,Ksi=-1,Heta = -1,Zeta = 0},
				    new Knot(){ID=1,Ksi=-1,Heta = 1,Zeta = 0},
				    new Knot(){ID=2,Ksi=1,Heta = -1,Zeta = 0},
				    new Knot(){ID=3,Ksi=1,Heta = 1,Zeta = 0}
				});

		    IVector parametricGaussPointKsi = new Vector(element.DegreeKsi + 1);
		    for (int i = 0; i < element.DegreeKsi + 1; i++)
		    {
			    parametricGaussPointKsi[i] = gaussPoints[i * (element.DegreeHeta + 1)].Ksi;
		    }

		    IVector parametricGaussPointHeta = new Vector(element.DegreeHeta + 1);
		    for (int i = 0; i < element.DegreeHeta + 1; i++)
		    {
			    parametricGaussPointHeta[i] = gaussPoints[i].Heta;
		    }

		    Vector knotValueVectorKsi = new Vector((element.DegreeKsi + 1) * 2);
		    Vector knotValueVectorHeta = new Vector((element.DegreeHeta + 1) * 2);
		    for (int i = 0; i < element.DegreeKsi + 1; i++)
		    {
			    knotValueVectorKsi[i]= -1;
			    knotValueVectorKsi[element.DegreeKsi + 1 + i]= 1;
		    }
		    for (int i = 0; i < element.DegreeHeta + 1; i++)
		    {
			    knotValueVectorHeta[i]= -1;
			    knotValueVectorHeta[element.DegreeHeta + 1 + i]= 1;
		    }

		    BSPLines1D bernsteinKsi = new BSPLines1D(element.DegreeKsi, knotValueVectorKsi, parametricGaussPointKsi);
		    BSPLines1D bernsteinHeta = new BSPLines1D(element.DegreeHeta, knotValueVectorHeta, parametricGaussPointHeta);
		    bernsteinKsi.calculateBSPLinesAndDerivatives();
		    bernsteinHeta.calculateBSPLinesAndDerivatives();

		    int supportKsi = element.DegreeKsi + 1;
		    int supportHeta = element.DegreeHeta + 1;

			Matrix2D bernsteinShapeFunctions = new Matrix2D(supportKsi*supportHeta,supportKsi*supportHeta);
		    Matrix2D bernsteinShapeFunctionDerivativesKsi = new Matrix2D(supportKsi * supportHeta, supportKsi * supportHeta);
		    Matrix2D bernsteinShapeFunctionDerivativesHeta = new Matrix2D(supportKsi * supportHeta, supportKsi * supportHeta);

		    for (int i = 0; i < supportKsi; i++)
		    {
			    for (int j = 0; j < supportHeta; j++)
			    {
				    for (int k = 0; k < supportKsi; k++)
				    {
					    for (int l = 0; l < supportHeta; l++)
					    {
						    bernsteinShapeFunctions[k * supportHeta + l, i * supportHeta + j] =
							    bernsteinKsi.BSPLineValues[k, i] * bernsteinHeta.BSPLineValues[l, j];
						    bernsteinShapeFunctionDerivativesKsi[k * supportHeta + l, i * supportHeta + j] =
							    bernsteinKsi.BSPLineDerivativeValues[k, i] * bernsteinHeta.BSPLineValues[l, j];
						    bernsteinShapeFunctionDerivativesHeta[k * supportHeta + l, i * supportHeta + j] =
							    bernsteinKsi.BSPLineValues[k, i] * bernsteinHeta.BSPLineDerivativeValues[l, j];
					    }
				    }
			    }   
		    }

			Matrix2D rationalTSplines = element.ExtractionOperator * bernsteinShapeFunctions.Transpose();
		    Matrix2D rationalTSplineDerivativesKsi = element.ExtractionOperator * bernsteinShapeFunctionDerivativesKsi.Transpose();
		    Matrix2D rationalTSplineDerivativesHeta = element.ExtractionOperator * bernsteinShapeFunctionDerivativesHeta.Transpose();

			TSplineValues= new Matrix2D(element.ControlPoints.Count, supportKsi*supportHeta);
			TSplineDerivativeValuesKsi = new Matrix2D(element.ControlPoints.Count, supportKsi * supportHeta);
			TSplineDerivativeValuesHeta = new Matrix2D(element.ControlPoints.Count, supportKsi * supportHeta);

		    for (int i = 0; i < supportKsi; i++)
		    {
			    for (int j = 0; j < supportHeta; j++)
			    {
				    double sumKsiHeta = 0;
				    double sumdKsiHeta = 0;
				    double sumKsidHeta = 0;
				    var index = i * supportHeta + j;

					for (int k = 0; k < element.ControlPoints.Count; k++)
				    {
					    sumKsiHeta += rationalTSplines[k, index] * element.ControlPoints[k].WeightFactor;
					    sumdKsiHeta += rationalTSplineDerivativesKsi[k, index] * element.ControlPoints[k].WeightFactor;
					    sumKsidHeta += rationalTSplineDerivativesHeta[k, index] * element.ControlPoints[k].WeightFactor;
					}

				    for (int k = 0; k < element.ControlPoints.Count; k++)
				    {
					    TSplineValues[k, index] = rationalTSplines[k, index] * element.ControlPoints[k].WeightFactor / sumKsiHeta;
					    TSplineDerivativeValuesKsi[k, index] = (rationalTSplineDerivativesKsi[k, index] * sumKsiHeta -
					                                            rationalTSplines[k, index] * sumdKsiHeta) /
					                                           Math.Pow(sumKsiHeta, 2) * element.ControlPoints[k].WeightFactor;
					    TSplineDerivativeValuesHeta[k, index] = (rationalTSplineDerivativesHeta[k, index] * sumKsiHeta -
					                                             rationalTSplines[k, index] * sumKsidHeta) /
					                                            Math.Pow(sumKsiHeta, 2) * element.ControlPoints[k].WeightFactor;
				    }
			    }
		    }
		}

		public ShapeTSplines2DFromBezierExtraction(TSplineKirchhoffLoveShellElement element, IList<ControlPoint> controlPoints)
		{
			GaussQuadrature gauss = new GaussQuadrature();
			IList<GaussLegendrePoint3D> gaussPoints = gauss.CalculateElementGaussPoints(element.DegreeKsi, element.DegreeHeta,
				new List<Knot>
				{
					new Knot(){ID=0,Ksi=-1,Heta = -1,Zeta = 0},
					new Knot(){ID=1,Ksi=-1,Heta = 1,Zeta = 0},
					new Knot(){ID=2,Ksi=1,Heta = -1,Zeta = 0},
					new Knot(){ID=3,Ksi=1,Heta = 1,Zeta = 0}
				});

			IVector parametricGaussPointKsi = new Vector(element.DegreeKsi + 1);
			for (int i = 0; i < element.DegreeKsi + 1; i++)
			{
				parametricGaussPointKsi[i] = gaussPoints[i * (element.DegreeHeta + 1)].Ksi;
			}

			IVector parametricGaussPointHeta = new Vector(element.DegreeHeta + 1);
			for (int i = 0; i < element.DegreeHeta + 1; i++)
			{
				parametricGaussPointHeta[i] = gaussPoints[i].Heta;
			}

			Vector knotValueVectorKsi = new Vector((element.DegreeKsi + 1) * 2);
			Vector knotValueVectorHeta = new Vector((element.DegreeHeta + 1) * 2);
			for (int i = 0; i < element.DegreeKsi + 1; i++)
			{
				knotValueVectorKsi[i] = -1;
				knotValueVectorKsi[element.DegreeKsi + 1 + i] = 1;
			}
			for (int i = 0; i < element.DegreeHeta + 1; i++)
			{
				knotValueVectorHeta[i] = -1;
				knotValueVectorHeta[element.DegreeHeta + 1 + i] = 1;
			}

			BSPLines1D bernsteinKsi = new BSPLines1D(element.DegreeKsi, knotValueVectorKsi, parametricGaussPointKsi);
			BSPLines1D bernsteinHeta = new BSPLines1D(element.DegreeHeta, knotValueVectorHeta, parametricGaussPointHeta);
			bernsteinKsi.calculateBSPLinesAndDerivatives();
			bernsteinHeta.calculateBSPLinesAndDerivatives();

			int supportKsi = element.DegreeKsi + 1;
			int supportHeta = element.DegreeHeta + 1;

			Matrix2D bernsteinShapeFunctions = new Matrix2D(supportKsi * supportHeta, supportKsi * supportHeta);
			Matrix2D bernsteinShapeFunctionDerivativesKsi = new Matrix2D(supportKsi * supportHeta, supportKsi * supportHeta);
			Matrix2D bernsteinShapeFunctionDerivativesHeta = new Matrix2D(supportKsi * supportHeta, supportKsi * supportHeta);
			Matrix2D bernsteinShapeFunctionSecondDerivativesKsi = new Matrix2D(supportKsi * supportHeta, supportKsi * supportHeta);
			Matrix2D bernsteinShapeFunctionSecondDerivativesHeta = new Matrix2D(supportKsi * supportHeta, supportKsi * supportHeta);
			Matrix2D bernsteinShapeFunctionSecondDerivativesKsiHeta = new Matrix2D(supportKsi * supportHeta, supportKsi * supportHeta);

			for (int i = 0; i < supportKsi; i++)
			{
				for (int j = 0; j < supportHeta; j++)
				{
					for (int k = 0; k < supportKsi; k++)
					{
						for (int l = 0; l < supportHeta; l++)
						{
							bernsteinShapeFunctions[k * supportHeta + l, i * supportHeta + j] =
								bernsteinKsi.BSPLineValues[k, i] * bernsteinHeta.BSPLineValues[l, j];
							bernsteinShapeFunctionDerivativesKsi[k * supportHeta + l, i * supportHeta + j] =
								bernsteinKsi.BSPLineDerivativeValues[k, i] * bernsteinHeta.BSPLineValues[l, j];
							bernsteinShapeFunctionDerivativesHeta[k * supportHeta + l, i * supportHeta + j] =
								bernsteinKsi.BSPLineValues[k, i] * bernsteinHeta.BSPLineDerivativeValues[l, j];
							bernsteinShapeFunctionSecondDerivativesKsi[k * supportHeta + l, i * supportHeta + j] =
								bernsteinKsi.BSPLineSecondDerivativeValues[k, i] * bernsteinHeta.BSPLineValues[l, j];
							bernsteinShapeFunctionSecondDerivativesHeta[k * supportHeta + l, i * supportHeta + j] =
								bernsteinKsi.BSPLineValues[k, i] * bernsteinHeta.BSPLineSecondDerivativeValues[l, j];
							bernsteinShapeFunctionSecondDerivativesKsiHeta[k * supportHeta + l, i * supportHeta + j] =
								bernsteinKsi.BSPLineDerivativeValues[k, i] * bernsteinHeta.BSPLineDerivativeValues[l, j];
						}
					}
				}
			}

			Matrix2D rationalTSplines = element.ExtractionOperator * bernsteinShapeFunctions.Transpose();
			Matrix2D rationalTSplineDerivativesKsi = element.ExtractionOperator * bernsteinShapeFunctionDerivativesKsi.Transpose();
			Matrix2D rationalTSplineDerivativesHeta = element.ExtractionOperator * bernsteinShapeFunctionDerivativesHeta.Transpose();
			Matrix2D rationalTSplineSecondDerivativesKsi = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesKsi.Transpose();
			Matrix2D rationalTSplineSecondDerivativesHeta = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesHeta.Transpose();
			Matrix2D rationalTSplineSecondDerivativesKsiHeta = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesKsiHeta.Transpose();

			TSplineValues = new Matrix2D(element.ControlPoints.Count, supportKsi * supportHeta);
			TSplineDerivativeValuesKsi = new Matrix2D(element.ControlPoints.Count, supportKsi * supportHeta);
			TSplineDerivativeValuesHeta = new Matrix2D(element.ControlPoints.Count, supportKsi * supportHeta);
			TSplineSecondDerivativesValueKsi = new Matrix2D(element.ControlPoints.Count, supportKsi * supportHeta);
			TSplineSecondDerivativesValueHeta = new Matrix2D(element.ControlPoints.Count, supportKsi * supportHeta);
			TSplineSecondDerivativesValueKsiHeta = new Matrix2D(element.ControlPoints.Count, supportKsi * supportHeta);

			for (int i = 0; i < supportKsi; i++)
			{
				for (int j = 0; j < supportHeta; j++)
				{
					double sumKsiHeta = 0;
					double sumdKsiHeta = 0;
					double sumKsidHeta = 0;
					double sumdKsidKsi = 0;
					double sumdHetadHeta = 0;
					double sumdKsidHeta = 0;

					var index = i * supportHeta + j;

					for (int k = 0; k < element.ControlPoints.Count; k++)
					{
						sumKsiHeta += rationalTSplines[k, index] * element.ControlPoints[k].WeightFactor;
						sumdKsiHeta += rationalTSplineDerivativesKsi[k, index] * element.ControlPoints[k].WeightFactor;
						sumKsidHeta += rationalTSplineDerivativesHeta[k, index] * element.ControlPoints[k].WeightFactor;
						sumdKsidKsi += rationalTSplineSecondDerivativesKsi[k, index] * element.ControlPoints[k].WeightFactor;
						sumdHetadHeta += rationalTSplineSecondDerivativesHeta[k, index] * element.ControlPoints[k].WeightFactor;
						sumdKsidHeta += rationalTSplineSecondDerivativesKsiHeta[k, index] * element.ControlPoints[k].WeightFactor;
					}

					for (int k = 0; k < element.ControlPoints.Count; k++)
					{
						TSplineValues[k, index] = rationalTSplines[k, index] * element.ControlPoints[k].WeightFactor / sumKsiHeta;
						TSplineDerivativeValuesKsi[k, index] = (rationalTSplineDerivativesKsi[k, index] * sumKsiHeta -
																rationalTSplines[k, index] * sumdKsiHeta) /
															   Math.Pow(sumKsiHeta, 2) * element.ControlPoints[k].WeightFactor;
						TSplineDerivativeValuesHeta[k, index] = (rationalTSplineDerivativesHeta[k, index] * sumKsiHeta -
																 rationalTSplines[k, index] * sumKsidHeta) /
																Math.Pow(sumKsiHeta, 2) * element.ControlPoints[k].WeightFactor;
						TSplineSecondDerivativesValueKsi[k, index] = (rationalTSplineSecondDerivativesKsi[k, index] / sumKsiHeta -
						                                              2 * rationalTSplineDerivativesKsi[k, index] * sumdKsiHeta /
						                                              Math.Pow(sumKsiHeta, 2) -
						                                              rationalTSplines[k, index] * sumdKsidKsi / Math.Pow(sumKsiHeta, 2) +
						                                              2 * rationalTSplines[k, index] * Math.Pow(sumdKsiHeta, 2) /
						                                              Math.Pow(sumKsiHeta, 3))* element.ControlPoints[k].WeightFactor;
						TSplineSecondDerivativesValueHeta[k, index] = (rationalTSplineSecondDerivativesHeta[k, index] / sumKsiHeta -
						                                               2 * rationalTSplineDerivativesHeta[k, index] * sumKsidHeta /
						                                               Math.Pow(sumKsiHeta, 2) -
						                                               rationalTSplines[k, index] * sumdHetadHeta /
						                                               Math.Pow(sumKsiHeta, 2) +
						                                               2 * rationalTSplines[k, index] * Math.Pow(sumKsidHeta, 2) /
						                                               Math.Pow(sumKsiHeta, 3)) * element.ControlPoints[k].WeightFactor;
						TSplineSecondDerivativesValueKsiHeta[k, index] = (rationalTSplineSecondDerivativesKsiHeta[k, index] / sumKsiHeta -
						                                                  rationalTSplineDerivativesKsi[k, index] * sumKsidHeta /
						                                                  Math.Pow(sumKsiHeta, 2) -
						                                                  rationalTSplineDerivativesHeta[k, index] * sumdKsiHeta /
						                                                  Math.Pow(sumKsiHeta, 2) -
						                                                  rationalTSplines[k, index] * sumdKsidHeta /
						                                                  Math.Pow(sumKsiHeta, 2) +
						                                                  2 * rationalTSplines[k, index] * sumdKsiHeta * sumKsidHeta /
						                                                  Math.Pow(sumKsiHeta, 3)) *
						                                                 element.ControlPoints[k].WeightFactor;
					}
				}
			}
		}
	}
}
