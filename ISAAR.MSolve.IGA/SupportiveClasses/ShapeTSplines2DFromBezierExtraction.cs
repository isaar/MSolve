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

		    var bKsi = MatrixPart(supportKsi, bernsteinKsi.BSPLineValues);
		    var bdKsi = MatrixPart(supportKsi, bernsteinKsi.BSPLineDerivativeValues);

		    var bheta = MatrixPart(supportHeta, bernsteinHeta.BSPLineValues);
		    var bdheta = MatrixPart(supportHeta, bernsteinHeta.BSPLineDerivativeValues);

		    var bernsteinShapeFunctions = KroneckerProduct(bKsi, bheta);
		    Matrix2D bernsteinShapeFunctionDerivativesKsi = KroneckerProduct(bheta, bdKsi);
		    Matrix2D bernsteinShapeFunctionDerivativesHeta = KroneckerProduct(bdheta, bKsi);

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

			var bKsi = MatrixPart(supportKsi, bernsteinKsi.BSPLineValues);
			var bdKsi = MatrixPart(supportKsi, bernsteinKsi.BSPLineDerivativeValues);
			var bddKsi = MatrixPart(supportKsi, bernsteinKsi.BSPLineSecondDerivativeValues);

			var bheta = MatrixPart(supportHeta, bernsteinHeta.BSPLineValues);
			var bdheta = MatrixPart(supportHeta, bernsteinHeta.BSPLineDerivativeValues);
			var bddheta = MatrixPart(supportHeta, bernsteinHeta.BSPLineSecondDerivativeValues);

			var bernsteinShapeFunctions = KroneckerProduct(bKsi, bheta);
			Matrix2D bernsteinShapeFunctionDerivativesKsi = KroneckerProduct(bheta, bdKsi);
			Matrix2D bernsteinShapeFunctionDerivativesHeta = KroneckerProduct(bdheta, bKsi);
			Matrix2D bernsteinShapeFunctionSecondDerivativesKsi = KroneckerProduct(bheta, bddKsi);
			Matrix2D bernsteinShapeFunctionSecondDerivativesHeta = KroneckerProduct(bddheta, bKsi);
			Matrix2D bernsteinShapeFunctionSecondDerivativesKsiHeta = KroneckerProduct(bdheta, bdKsi);

			
			Matrix2D rationalTSplines = element.ExtractionOperator * bernsteinShapeFunctions;
			Matrix2D rationalTSplineDerivativesKsi = element.ExtractionOperator * bernsteinShapeFunctionDerivativesKsi;
			Matrix2D rationalTSplineDerivativesHeta = element.ExtractionOperator * bernsteinShapeFunctionDerivativesHeta;
			Matrix2D rationalTSplineSecondDerivativesKsi = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesKsi;
			Matrix2D rationalTSplineSecondDerivativesHeta = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesHeta;
			Matrix2D rationalTSplineSecondDerivativesKsiHeta = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesKsiHeta;

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


		public ShapeTSplines2DFromBezierExtraction(TSplineKirchhoffLoveShellElement element, IList<ControlPoint> controlPoints, Vector parametricGaussPointKsi, Vector parametricGaussPointHeta)
		{

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

			int supportElementKsi = element.DegreeKsi + 1;
			int supportElementHeta = element.DegreeHeta + 1;
			int supportKsi = parametricGaussPointKsi.Length;
			int supportHeta = parametricGaussPointHeta.Length;

			var bKsi = MatrixPart(supportElementKsi, supportKsi, bernsteinKsi.BSPLineValues);
			var bdKsi = MatrixPart(supportElementKsi, supportKsi, bernsteinKsi.BSPLineDerivativeValues);
			var bddKsi = MatrixPart(supportElementKsi, supportKsi, bernsteinKsi.BSPLineSecondDerivativeValues);

			var bheta = MatrixPart(supportElementHeta, supportHeta, bernsteinHeta.BSPLineValues);
			var bdheta = MatrixPart(supportElementHeta, supportHeta, bernsteinHeta.BSPLineDerivativeValues);
			var bddheta = MatrixPart(supportElementHeta, supportHeta, bernsteinHeta.BSPLineSecondDerivativeValues);

			var bernsteinShapeFunctions = KroneckerProduct(bKsi, bheta);
			Matrix2D bernsteinShapeFunctionDerivativesKsi = KroneckerProduct(bheta, bdKsi);
			Matrix2D bernsteinShapeFunctionDerivativesHeta = KroneckerProduct(bdheta, bKsi);
			Matrix2D bernsteinShapeFunctionSecondDerivativesKsi = KroneckerProduct(bheta, bddKsi);
			Matrix2D bernsteinShapeFunctionSecondDerivativesHeta = KroneckerProduct(bddheta, bKsi);
			Matrix2D bernsteinShapeFunctionSecondDerivativesKsiHeta = KroneckerProduct(bdheta, bdKsi);


			Matrix2D rationalTSplines = element.ExtractionOperator * bernsteinShapeFunctions;
			Matrix2D rationalTSplineDerivativesKsi = element.ExtractionOperator * bernsteinShapeFunctionDerivativesKsi;
			Matrix2D rationalTSplineDerivativesHeta = element.ExtractionOperator * bernsteinShapeFunctionDerivativesHeta;
			Matrix2D rationalTSplineSecondDerivativesKsi = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesKsi;
			Matrix2D rationalTSplineSecondDerivativesHeta = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesHeta;
			Matrix2D rationalTSplineSecondDerivativesKsiHeta = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesKsiHeta;

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
																	  Math.Pow(sumKsiHeta, 3)) * element.ControlPoints[k].WeightFactor;
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

		public ShapeTSplines2DFromBezierExtraction(TSplineKirchhoffLoveShellElementMaterial element, IList<ControlPoint> controlPoints, Vector parametricGaussPointKsi, Vector parametricGaussPointHeta)
		{

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

			int supportElementKsi = element.DegreeKsi + 1;
			int supportElementHeta = element.DegreeHeta + 1;
			int supportKsi = parametricGaussPointKsi.Length;
			int supportHeta = parametricGaussPointHeta.Length;

			var bKsi = MatrixPart(supportElementKsi, supportKsi, bernsteinKsi.BSPLineValues);
			var bdKsi = MatrixPart(supportElementKsi, supportKsi, bernsteinKsi.BSPLineDerivativeValues);
			var bddKsi = MatrixPart(supportElementKsi, supportKsi, bernsteinKsi.BSPLineSecondDerivativeValues);

			var bheta = MatrixPart(supportElementHeta, supportHeta, bernsteinHeta.BSPLineValues);
			var bdheta = MatrixPart(supportElementHeta, supportHeta, bernsteinHeta.BSPLineDerivativeValues);
			var bddheta = MatrixPart(supportElementHeta, supportHeta, bernsteinHeta.BSPLineSecondDerivativeValues);

			var bernsteinShapeFunctions = KroneckerProduct(bKsi, bheta);
			Matrix2D bernsteinShapeFunctionDerivativesKsi = KroneckerProduct(bheta, bdKsi);
			Matrix2D bernsteinShapeFunctionDerivativesHeta = KroneckerProduct(bdheta, bKsi);
			Matrix2D bernsteinShapeFunctionSecondDerivativesKsi = KroneckerProduct(bheta, bddKsi);
			Matrix2D bernsteinShapeFunctionSecondDerivativesHeta = KroneckerProduct(bddheta, bKsi);
			Matrix2D bernsteinShapeFunctionSecondDerivativesKsiHeta = KroneckerProduct(bdheta, bdKsi);


			Matrix2D rationalTSplines = element.ExtractionOperator * bernsteinShapeFunctions;
			Matrix2D rationalTSplineDerivativesKsi = element.ExtractionOperator * bernsteinShapeFunctionDerivativesKsi;
			Matrix2D rationalTSplineDerivativesHeta = element.ExtractionOperator * bernsteinShapeFunctionDerivativesHeta;
			Matrix2D rationalTSplineSecondDerivativesKsi = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesKsi;
			Matrix2D rationalTSplineSecondDerivativesHeta = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesHeta;
			Matrix2D rationalTSplineSecondDerivativesKsiHeta = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesKsiHeta;

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
																	  Math.Pow(sumKsiHeta, 3)) * element.ControlPoints[k].WeightFactor;
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


		public ShapeTSplines2DFromBezierExtraction(TSplineElement2D element, IList<ControlPoint> controlPoints, Vector parametricGaussPointKsi, Vector parametricGaussPointHeta)
		{

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

			int supportElementKsi = element.DegreeKsi + 1;
			int supportElementHeta = element.DegreeHeta + 1;
			int supportKsi = parametricGaussPointKsi.Length;
			int supportHeta = parametricGaussPointHeta.Length;

			var bKsi = MatrixPart(supportElementKsi, supportKsi, bernsteinKsi.BSPLineValues);
			var bdKsi = MatrixPart(supportElementKsi, supportKsi, bernsteinKsi.BSPLineDerivativeValues);
			var bddKsi = MatrixPart(supportElementKsi, supportKsi, bernsteinKsi.BSPLineSecondDerivativeValues);

			var bheta = MatrixPart(supportElementHeta, supportHeta, bernsteinHeta.BSPLineValues);
			var bdheta = MatrixPart(supportElementHeta, supportHeta, bernsteinHeta.BSPLineDerivativeValues);
			var bddheta = MatrixPart(supportElementHeta, supportHeta, bernsteinHeta.BSPLineSecondDerivativeValues);

			var bernsteinShapeFunctions = KroneckerProduct(bKsi, bheta);
			Matrix2D bernsteinShapeFunctionDerivativesKsi = KroneckerProduct(bheta, bdKsi);
			Matrix2D bernsteinShapeFunctionDerivativesHeta = KroneckerProduct(bdheta, bKsi);
			Matrix2D bernsteinShapeFunctionSecondDerivativesKsi = KroneckerProduct(bheta, bddKsi);
			Matrix2D bernsteinShapeFunctionSecondDerivativesHeta = KroneckerProduct(bddheta, bKsi);
			Matrix2D bernsteinShapeFunctionSecondDerivativesKsiHeta = KroneckerProduct(bdheta, bdKsi);


			Matrix2D rationalTSplines = element.ExtractionOperator * bernsteinShapeFunctions;
			Matrix2D rationalTSplineDerivativesKsi = element.ExtractionOperator * bernsteinShapeFunctionDerivativesKsi;
			Matrix2D rationalTSplineDerivativesHeta = element.ExtractionOperator * bernsteinShapeFunctionDerivativesHeta;
			Matrix2D rationalTSplineSecondDerivativesKsi = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesKsi;
			Matrix2D rationalTSplineSecondDerivativesHeta = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesHeta;
			Matrix2D rationalTSplineSecondDerivativesKsiHeta = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesKsiHeta;

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
																	  Math.Pow(sumKsiHeta, 3)) * element.ControlPoints[k].WeightFactor;
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


		public ShapeTSplines2DFromBezierExtraction(TSplineKirchhoffLoveShellElementMaterial element, IList<ControlPoint> controlPoints)
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

			var bKsi = MatrixPart(supportKsi, bernsteinKsi.BSPLineValues);
			var bdKsi = MatrixPart(supportKsi, bernsteinKsi.BSPLineDerivativeValues);
			var bddKsi = MatrixPart(supportKsi, bernsteinKsi.BSPLineSecondDerivativeValues);

			var bheta = MatrixPart(supportHeta, bernsteinHeta.BSPLineValues);
			var bdheta = MatrixPart(supportHeta, bernsteinHeta.BSPLineDerivativeValues);
			var bddheta = MatrixPart(supportHeta, bernsteinHeta.BSPLineSecondDerivativeValues);

			var bernsteinShapeFunctions = KroneckerProduct(bKsi, bheta);
			Matrix2D bernsteinShapeFunctionDerivativesKsi = KroneckerProduct(bheta, bdKsi);
			Matrix2D bernsteinShapeFunctionDerivativesHeta = KroneckerProduct(bdheta, bKsi);
			Matrix2D bernsteinShapeFunctionSecondDerivativesKsi = KroneckerProduct(bheta, bddKsi);
			Matrix2D bernsteinShapeFunctionSecondDerivativesHeta = KroneckerProduct(bddheta, bKsi);
			Matrix2D bernsteinShapeFunctionSecondDerivativesKsiHeta = KroneckerProduct(bdheta, bdKsi);


			Matrix2D rationalTSplines = element.ExtractionOperator * bernsteinShapeFunctions;
			Matrix2D rationalTSplineDerivativesKsi = element.ExtractionOperator * bernsteinShapeFunctionDerivativesKsi;
			Matrix2D rationalTSplineDerivativesHeta = element.ExtractionOperator * bernsteinShapeFunctionDerivativesHeta;
			Matrix2D rationalTSplineSecondDerivativesKsi = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesKsi;
			Matrix2D rationalTSplineSecondDerivativesHeta = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesHeta;
			Matrix2D rationalTSplineSecondDerivativesKsiHeta = element.ExtractionOperator * bernsteinShapeFunctionSecondDerivativesKsiHeta;

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
																	  Math.Pow(sumKsiHeta, 3)) * element.ControlPoints[k].WeightFactor;
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


		private static Matrix2D MatrixPart(int support, double[,] matrix)
	    {
		    var A = new Matrix2D(support, support);
		    for (int i = 0; i < support; i++)
		    {
			    for (int j = 0; j < support; j++)
			    {
				    A[i, j] = matrix[i, j];
			    }
		    }

		    return A;
	    }

	    private static Matrix2D MatrixPart(int support1, int support2, double[,] matrix)
	    {
		    var A = new Matrix2D(support1, support2);
		    for (int i = 0; i < support1; i++)
		    {
			    for (int j = 0; j < support2; j++)
			    {
				    A[i, j] = matrix[i, j];
			    }
		    }

		    return A;
	    }


		private static Matrix2D KroneckerProduct(Matrix2D A, Matrix2D B)
	    {
			Matrix2D C=new Matrix2D(A.Rows*B.Rows,A.Columns*B.Columns);
		    for (int rowAIndex = 0; rowAIndex < A.Rows; rowAIndex++)
		    {
			    for (int rowBIndex = 0; rowBIndex < B.Rows; rowBIndex++)
			    {
				    for (int columnAIndex = 0; columnAIndex < A.Columns; columnAIndex++)
				    {
					    for (int columnBIndex = 0; columnBIndex < B.Columns; columnBIndex++)
					    {
							C[rowAIndex*B.Rows + rowBIndex, columnAIndex*B.Columns + columnBIndex] =
							    A[rowAIndex, columnAIndex] * B[rowBIndex, columnBIndex];
					    }
				    }
			    }
		    }

		    return C;
	    }
	}
}
