using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Interfaces;
using IVector = ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces.IVector;
using Vector = ISAAR.MSolve.Numerical.LinearAlgebra.Vector;

namespace ISAAR.MSolve.IGA.Postprocessing
{
	public class ParaviewNurbs3D
	{
		private Model _model;
		private IVectorView _solution;
		private string _filename;
		public ParaviewNurbs3D(Model model, IVectorView solution, string filename)
		{
			_model = model;
			_solution = solution;
			_filename = filename;
		}

		public void CreateParaviewFile()
		{
			var uniqueKnotsKsi = _model.PatchesDictionary[0].KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity();
			var uniqueKnotsHeta = _model.PatchesDictionary[0].KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity();
			var uniqueKnotsZeta = _model.PatchesDictionary[0].KnotValueVectorZeta.RemoveDuplicatesFindMultiplicity();

			var numberOfKnotsKsi = uniqueKnotsKsi[0].Length;
			var numberOfKnotsHeta = uniqueKnotsHeta[0].Length;
			var numberOfKnotsZeta = uniqueKnotsZeta[0].Length;
			
			var knots = new double[numberOfKnotsKsi * numberOfKnotsHeta*numberOfKnotsZeta, 3];
			var count = 0;
			var patch = _model.PatchesDictionary[0];

			var projectiveControlPoints = CalculateProjectiveControlPoints();

			for (int knotKsiIndex = 0; knotKsiIndex < numberOfKnotsKsi; knotKsiIndex++)
			{
				for (int knotHetaIndex = 0; knotHetaIndex < numberOfKnotsHeta; knotHetaIndex++)
				{
					for (int knotZetaIndex = 0; knotZetaIndex < numberOfKnotsZeta; knotZetaIndex++)
					{
						var ksiCoordinate = uniqueKnotsKsi[0][knotKsiIndex];
						var hetaCoordinate = uniqueKnotsHeta[0][knotHetaIndex];
						var zetaCoordinate = uniqueKnotsZeta[0][knotZetaIndex];

						var point3D = SolidPoint3D(
							patch.NumberOfControlPointsKsi - 1, patch.DegreeKsi, patch.KnotValueVectorKsi,
							patch.NumberOfControlPointsHeta - 1, patch.DegreeHeta, patch.KnotValueVectorHeta,
							patch.NumberOfControlPointsZeta-1, patch.DegreeZeta,patch.KnotValueVectorZeta,
							projectiveControlPoints, ksiCoordinate, hetaCoordinate, zetaCoordinate);

						knots[count, 0] = point3D[0] / point3D[3];
						knots[count, 1] = point3D[1] / point3D[3];
						knots[count++, 2] = point3D[2] / point3D[3];
					}
				}
			}

			var numberOfElements = (numberOfKnotsKsi - 1) * (numberOfKnotsHeta - 1) * (numberOfKnotsZeta - 1);
			var elementConnectivity = new int[numberOfElements, 8];

			foreach (var element in _model.Elements)
				for (int i = 0; i < element.Knots.Count; i++)
					elementConnectivity[element.ID, i] = element.Knots[i].ID;

			var knotDisplacements = new double[knots.GetLength(0), 3];
			foreach (var element in _model.Elements)
			{
				var localDisplacements = new double[element.ControlPoints.Count, 3];
				var counterCP = 0;
				foreach (var controlPoint in element.ControlPoints)
				{
					//var dofX = _model.ControlPointDOFsDictionary[controlPoint.ID][DOFType.X];
					//var dofY = _model.ControlPointDOFsDictionary[controlPoint.ID][DOFType.Y];
					//var dofZ = _model.ControlPointDOFsDictionary[controlPoint.ID][DOFType.Z];


					var dofX = _model.GlobalDofOrdering.GlobalFreeDofs[controlPoint, DOFType.X];
					var dofY = _model.GlobalDofOrdering.GlobalFreeDofs[controlPoint, DOFType.Y];
					var dofZ = _model.GlobalDofOrdering.GlobalFreeDofs[controlPoint, DOFType.Y];
					localDisplacements[counterCP, 0] = (dofX == -1) ? 0.0 : _solution[dofX];
					localDisplacements[counterCP, 1] = (dofY == -1) ? 0.0 : _solution[dofY];
					localDisplacements[counterCP++, 2] = (dofZ == -1) ? 0.0 : _solution[dofZ];
				}
				var elementKnotDisplacements = element.ElementType.CalculateDisplacementsForPostProcessing(element, localDisplacements);
				for (int i = 0; i < elementConnectivity.GetLength(1); i++)
				{
					var knotConnectivity = elementConnectivity[element.ID, i];
					knotDisplacements[knotConnectivity, 0] = elementKnotDisplacements[i, 0];
					knotDisplacements[knotConnectivity, 1] = elementKnotDisplacements[i, 1];
					knotDisplacements[knotConnectivity, 2] = elementKnotDisplacements[i, 2];
				}
			}

			WriteParaviewFile3D(knots, knotDisplacements);
		}

		public void WriteParaviewFile3D(double[,] nodeCoordinates, double[,] displacements)
		{
			var patch = _model.PatchesDictionary[0];

			var numberOfKnotsKsi = patch.KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[0].Length;
			var numberOfKnotsHeta = patch.KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[0].Length;
			var numberOfKnotZeta = patch.KnotValueVectorZeta.RemoveDuplicatesFindMultiplicity()[0].Length;

			var numberOfNodes = nodeCoordinates.GetLength(0);

			using (StreamWriter outputFile = new StreamWriter($"..\\..\\..\\OutputFiles\\{_filename}Paraview.vts"))
			{
				outputFile.WriteLine("<?xml version=\"1.0\"?>");
				outputFile.WriteLine("<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\" >");
				outputFile.WriteLine($"<StructuredGrid  WholeExtent=\"{0} {numberOfKnotZeta - 1} {0} {numberOfKnotsHeta-1} {0} {numberOfKnotsKsi - 1}\">");
				outputFile.WriteLine($"<Piece Extent=\"{0} {numberOfKnotZeta - 1} {0} {numberOfKnotsHeta-1} {0} {numberOfKnotsKsi - 1}\">");

				outputFile.WriteLine("<PointData Vectors=\"Disp\"  >");
				outputFile.WriteLine("<DataArray type=\"Float32\" Name=\"Displacement\" NumberOfComponents=\"3\" format=\"ascii\">");
				for (int i = 0; i < numberOfNodes; i++)
					outputFile.WriteLine($"{displacements[i, 0]} {displacements[i, 1]} {displacements[i, 2]}");
				outputFile.WriteLine("</DataArray>");
				outputFile.WriteLine("</PointData>");

				outputFile.WriteLine("<Celldata>");
				outputFile.WriteLine("</Celldata>");
				outputFile.WriteLine("<Points>");
				outputFile.WriteLine("<DataArray type=\"Float32\" Name=\"Array\" NumberOfComponents=\"3\" format=\"ascii\">");
				for (int i = 0; i < numberOfNodes; i++)
					outputFile.WriteLine($"{nodeCoordinates[i, 0]} {nodeCoordinates[i, 1]} {nodeCoordinates[i, 2]}");
				outputFile.WriteLine("</DataArray>");
				outputFile.WriteLine("</Points>");

				outputFile.WriteLine("</Piece>");
				outputFile.WriteLine("</StructuredGrid>");
				outputFile.WriteLine("</VTKFile>");

			}

		}

		public static Vector SolidPoint3D(int numberOfCPKsi, int degreeKsi, IVector knotValueVectorKsi,
			int numberOfCPHeta, int degreeHeta, IVector knotValueVectorHeta, int numberOfCPZeta, int degreeZeta,
			IVector knotValueVectorZeta, double[,] projectiveControlPointCoordinates, double ksiCoordinate, double hetaCoordinate,
			double zetaCoordinate)
		{
			var spanKsi = ParaviewNurbs2D.FindSpan(numberOfCPKsi, degreeKsi, ksiCoordinate, knotValueVectorKsi);
			var spanHeta = ParaviewNurbs2D.FindSpan(numberOfCPHeta, degreeHeta, hetaCoordinate, knotValueVectorHeta);
			var spanZeta = ParaviewNurbs2D.FindSpan(numberOfCPZeta, degreeZeta, zetaCoordinate, knotValueVectorZeta);

			var pointFunctionsKsi = ParaviewNurbs2D.BasisFunctions(spanKsi, ksiCoordinate, degreeKsi, knotValueVectorKsi);
			var pointFunctionsHeta = ParaviewNurbs2D.BasisFunctions(spanHeta, hetaCoordinate, degreeHeta, knotValueVectorHeta);
			var pointFunctionsZeta = ParaviewNurbs2D.BasisFunctions(spanZeta, zetaCoordinate, degreeZeta, knotValueVectorZeta);

			var cartesianPoint = new Vector(4);
			var indexKsi = spanKsi - degreeKsi;

			for (int k = 0; k <= degreeZeta; k++)
			{
				var indexZeta = spanZeta - degreeZeta + k;
				for (int j = 0; j <= degreeHeta; j++)
				{
					var indexHeta = spanHeta - degreeHeta + j;
					for (int i = 0; i <= degreeKsi; i++)
					{
						var cpIndex = (indexKsi + i) * (numberOfCPHeta + 1) * (numberOfCPZeta + 1) +
						              indexHeta * (numberOfCPZeta + 1) + indexZeta;

						var cpCoordinates = new Vector(new double[]
						{
							projectiveControlPointCoordinates[cpIndex, 0],
							projectiveControlPointCoordinates[cpIndex, 1],
							projectiveControlPointCoordinates[cpIndex, 2],
							projectiveControlPointCoordinates[cpIndex, 3]
						});
						cpCoordinates.Scale(pointFunctionsKsi[i]);
						cpCoordinates.Scale(pointFunctionsHeta[j]);
						cpCoordinates.Scale(pointFunctionsZeta[k]);
						cartesianPoint= new Vector(cartesianPoint+cpCoordinates);
					}
				}
			}

			return cartesianPoint;
		}

		private double[,] CalculateProjectiveControlPoints()
		{
			var projectiveCPs = new double[_model.PatchesDictionary[0].ControlPoints.Count, 4];
			foreach (var controlPoint in _model.PatchesDictionary[0].ControlPoints)
			{
				var weight = controlPoint.WeightFactor;
				projectiveCPs[controlPoint.ID, 0] = controlPoint.X * weight;
				projectiveCPs[controlPoint.ID, 1] = controlPoint.Y * weight;
				projectiveCPs[controlPoint.ID, 2] = controlPoint.Z * weight;
				projectiveCPs[controlPoint.ID, 3] = weight;
			}

			return projectiveCPs;
		}
	}
}
