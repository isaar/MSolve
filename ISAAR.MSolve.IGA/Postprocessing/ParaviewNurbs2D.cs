using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Threading;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Interfaces;
using IVector = ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces.IVector;
using Vector = ISAAR.MSolve.Numerical.LinearAlgebra.Vector;

namespace ISAAR.MSolve.IGA.Postprocessing
{
	public class ParaviewNurbs2D
	{
		private Model _model;
		private IVectorView _solution;
		private string _filename;

		public ParaviewNurbs2D(Model model, IVectorView solution, string filename)
		{
			_model = model;
			_solution = solution;
			_filename = filename;
		}

		public void CreateParaview2DFile()
		{
			var uniqueKnotsKsi = _model.PatchesDictionary[0].KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity();
			var uniqueKnotsHeta = _model.PatchesDictionary[0].KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity();

			var numberOfKnotsKsi = uniqueKnotsKsi[0].Length;
			var numberOfKnotsHeta = uniqueKnotsHeta[0].Length;

			var knots = new double[numberOfKnotsKsi * numberOfKnotsHeta, 3];
			var count = 0;
			var patch = _model.PatchesDictionary[0];

			var projectiveControlPoints = CalculateProjectiveControlPoints();

			for (var knotKsiIndex = 0; knotKsiIndex < numberOfKnotsKsi; knotKsiIndex++)
			{
				for (var knotHetaIndex = 0; knotHetaIndex < numberOfKnotsHeta; knotHetaIndex++)
				{
					var hetaCoordinate = uniqueKnotsHeta[0][knotHetaIndex];
					var ksiCoordinate = uniqueKnotsKsi[0][knotKsiIndex];
					var point3D = SurfacePoint2D(patch.NumberOfControlPointsKsi - 1, patch.DegreeKsi,
						patch.KnotValueVectorKsi, patch.NumberOfControlPointsHeta-1, patch.DegreeHeta,
						patch.KnotValueVectorHeta, projectiveControlPoints, ksiCoordinate, hetaCoordinate);
					knots[count, 0] = point3D[0] / point3D[3];
					knots[count, 1] = point3D[1] / point3D[3];
					knots[count++, 2] = point3D[2] / point3D[3];
				}
			}

			var incrementKsi = numberOfKnotsHeta;
			var incrementHeta = 1;
			var nodePattern = new int[] {0, incrementKsi, incrementKsi + 1, 1};
			var elementConnectivity = CreateElement2DConnectivity(nodePattern, uniqueKnotsKsi[0].Length-1,
				uniqueKnotsHeta[0].Length-1, incrementKsi, incrementHeta);
			var knotDisplacements= new double[knots.GetLength(0),2];


			foreach (var element in _model.Elements)
			{
				var localDisplacements = new double[element.ControlPoints.Count, 2];
				var counterCP = 0;
				foreach (var controlPoint in element.ControlPoints)
				{
					localDisplacements[counterCP, 0] =
						(!_model.GlobalDofOrdering.GlobalFreeDofs.Contains(controlPoint, DOFType.X))
							? 0.0
							: _solution[_model.GlobalDofOrdering.GlobalFreeDofs[controlPoint, DOFType.X]];
					localDisplacements[counterCP++, 1] = 
						(!_model.GlobalDofOrdering.GlobalFreeDofs.Contains(controlPoint, DOFType.Y))
						? 0.0
						: _solution[_model.GlobalDofOrdering.GlobalFreeDofs[controlPoint, DOFType.Y]];
				}
				var elementKnotDisplacements=element.ElementType.CalculateDisplacementsForPostProcessing(element, localDisplacements);
				for (int i = 0; i < elementConnectivity.GetLength(1); i++)
				{
					var knotConnectivity = elementConnectivity[element.ID, i];
					knotDisplacements[knotConnectivity, 0] = elementKnotDisplacements[i, 0];
					knotDisplacements[knotConnectivity, 1] = elementKnotDisplacements[i, 1];
				}
			}

			Write2DNurbsFile(knots, elementConnectivity,"Quad4",knotDisplacements);
		}

		public void Write2DNurbsFile(double[,] nodeCoordinates, int[,] elementConnectivity, string elementType,double[,] displacements )
		{
			var dimensions = 2;
			var numberOfNodes = nodeCoordinates.GetLength(0);
			var numberOfCells = elementConnectivity.GetLength(0);

			int numberOfVerticesPerCell=0;
			int paraviewCellCode=0;

			if (elementType == "Quad4")
			{
				numberOfVerticesPerCell = 4;
				paraviewCellCode = 9;
			}

			var dofPerVertex = 2;
			using (StreamWriter outputFile = new StreamWriter($"..\\..\\..\\OutputFiles\\{_filename}Paraview.vtu"))
			{
				outputFile.WriteLine("<VTKFile type=\"UnstructuredGrid\"  version=\"0.1\"   >");
				outputFile.WriteLine("<UnstructuredGrid>");
				outputFile.WriteLine($"<Piece  NumberOfPoints=\"{numberOfNodes}\" NumberOfCells=\"{numberOfCells}\">");

				outputFile.WriteLine("<Points>");
				outputFile.WriteLine("<DataArray  type=\"Float64\"  NumberOfComponents=\"3\"  format=\"ascii\" >");
				for (int i = 0; i < numberOfNodes; i++)
					outputFile.WriteLine($"{nodeCoordinates[i,0]} {nodeCoordinates[i, 1]} {nodeCoordinates[i, 2]}");

				outputFile.WriteLine("</DataArray>");
				outputFile.WriteLine("</Points>");

				outputFile.WriteLine("<Cells>");
				outputFile.WriteLine("<DataArray  type=\"Int32\"  Name=\"connectivity\"  format=\"ascii\">");
				for (int i = 0; i < numberOfCells; i++)
				{
					for (int j = 0; j < elementConnectivity.GetLength(1); j++)
						outputFile.Write($"{elementConnectivity[i,j]} ");
					outputFile.WriteLine("");
				}

				outputFile.WriteLine("</DataArray>");
				outputFile.WriteLine("<DataArray  type=\"Int32\"  Name=\"offsets\"  format=\"ascii\">");

				var offset = 0;
				for (int i = 0; i < numberOfCells; i++)
				{
					offset += numberOfVerticesPerCell;
					outputFile.WriteLine(offset);
				}

				outputFile.WriteLine("</DataArray>");
				outputFile.WriteLine("<DataArray  type=\"UInt8\"  Name=\"types\"  format=\"ascii\">");
				for (int i = 0; i < numberOfCells; i++)
					outputFile.WriteLine(paraviewCellCode);

				outputFile.WriteLine("</DataArray>");
				outputFile.WriteLine("</Cells>");

				outputFile.WriteLine("<PointData  Vectors=\"U\">");

				outputFile.WriteLine("<DataArray  type=\"Float64\"  Name=\"U\" NumberOfComponents=\"3\" format=\"ascii\">");

				for (int i = 0; i < numberOfNodes; i++)
					outputFile.WriteLine($"{displacements[i,0]} {displacements[i, 1]} 0.0");

				outputFile.WriteLine("</DataArray>");
				outputFile.WriteLine("</PointData>");
				outputFile.WriteLine("</Piece>");
				outputFile.WriteLine("</UnstructuredGrid>");
				outputFile.WriteLine("</VTKFile>");
			}
		}

		private int[,] CreateElement2DConnectivity(int[] nodePattern, int numberOfElementsKsi, int numberOfElementsHeta,
			int incrementKsi, int incrementHeta)
		{
			var increment = 0;
			var elementConnectivity = new int[numberOfElementsKsi * numberOfElementsHeta, nodePattern.Length];
			var elementCounter = 0;
			for (int elementKsi = 0; elementKsi < numberOfElementsKsi; elementKsi++) 
			{
				increment = elementKsi * incrementKsi;
				for (int elementHeta = 0; elementHeta < numberOfElementsHeta; elementHeta++)
				{
					for (int i = 0; i < nodePattern.Length; i++)
						elementConnectivity[elementCounter, i] = nodePattern[i] + increment;
					increment += incrementHeta;
					elementCounter++;
				}
			}

			return elementConnectivity;
		}

		/// <summary>
		/// Creates a control point coordinates matrix in projective coordinates
		/// </summary>
		/// <returns></returns>
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

		/// <summary>
		/// NURBS Book Algorithm A2.1
		/// Algorithm that performs binary search to locate the knot span
		/// in which a point is located in the knot value vector
		/// </summary>
		/// <param name="numberOfBasisFunctions"></param>
		/// <param name="degree"></param>
		/// <param name="pointCoordinate"></param>
		/// <param name="knotValueVector"></param>
		/// <returns></returns>
		public static int FindSpan(int numberOfBasisFunctions, int degree, double pointCoordinate, IVector knotValueVector)
		{
			if (pointCoordinate == knotValueVector[numberOfBasisFunctions + 1]) return numberOfBasisFunctions;
			int minimum = degree;
			int maximum = numberOfBasisFunctions + 1;
			int mid = (minimum + maximum) / 2;
			while (pointCoordinate < knotValueVector[mid] || pointCoordinate >= knotValueVector[mid + 1])
			{
				if (pointCoordinate < knotValueVector[mid])
					maximum = mid;
				else
					minimum = mid;
				mid = (minimum + maximum) / 2;
			}

			return mid;
		}

		/// <summary>
		/// NURBS Book Algorithm A2.2
		/// Optimized Algorithm that calculates the basis functions.
		/// </summary>
		/// <param name="spanId"></param>
		/// <param name="pointCoordinate"></param>
		/// <param name="degree"></param>
		/// <param name="knotValueVector"></param>
		/// <returns></returns>
		public static Vector BasisFunctions(int spanId, double pointCoordinate, int degree, IVector knotValueVector)
		{
			var basisFunctions = new Vector(degree + 1);
			var left = new Vector(degree + 1);
			var right = new Vector(degree + 1);
			basisFunctions[0] = 1;
			for (int j = 1; j <= degree; j++)
			{
				left[j] = pointCoordinate - knotValueVector[spanId + 1 - j];
				right[j] = knotValueVector[spanId + j] - pointCoordinate;
				var saved = 0.0;
				for (int r = 0; r < j; r++)
				{
					var temp = basisFunctions[r] / (right[r + 1] + left[j - r]);
					basisFunctions[r] = saved + right[r + 1] * temp;
					saved = left[j - r] * temp;
				}

				basisFunctions[j] = saved;
			}

			return basisFunctions;
		}


		/// <summary>
		/// NURBS Book Algoritmh A3.5.
		/// Based on IGAFEM code.
		/// </summary>
		/// <param name="numberOfCPKsi"></param>
		/// <param name="degreeKsi"></param>
		/// <param name="knotValueVectorKsi"></param>
		/// <param name="numberOfCPHeta"></param>
		/// <param name="degreeHeta"></param>
		/// <param name="knotValueVectorHeta"></param>
		/// <param name="controlPoints"></param>
		/// <param name="coordinateKsi"></param>
		/// <param name="coordinateHeta"></param>
		/// <returns></returns>
		public static Vector SurfacePoint2D(int numberOfCPKsi, int degreeKsi, IVector knotValueVectorKsi, int numberOfCPHeta,
			int degreeHeta, IVector knotValueVectorHeta, double[,] projectiveControlPointCoordinates,
			double coordinateKsi, double coordinateHeta)
		{
			var spanKsi = FindSpan(numberOfCPKsi, degreeKsi, coordinateKsi, knotValueVectorKsi);
			var spanHeta = FindSpan(numberOfCPHeta, degreeHeta, coordinateHeta, knotValueVectorHeta);

			var pointFunctionsKsi = BasisFunctions(spanKsi, coordinateKsi, degreeKsi, knotValueVectorKsi);
			var pointFunctionsHeta = BasisFunctions(spanHeta, coordinateHeta, degreeHeta, knotValueVectorHeta);

			var cartesianPoint = new Vector(4);
			var indexKsi = spanKsi - degreeKsi;

			for (var j = 0; j <= degreeHeta; j++)
			{
				var temp = new Vector(4);
				var indexHeta = spanHeta - degreeHeta + j;

				for (int i = 0; i <= degreeKsi; i++)
				{
					var cpIndex = (indexKsi+i)*(numberOfCPHeta+1)+indexHeta; 
					var cpCoordinates = new Vector(new double[]
					{
						projectiveControlPointCoordinates[cpIndex, 0],
						projectiveControlPointCoordinates[cpIndex, 1],
						projectiveControlPointCoordinates[cpIndex, 2],
						projectiveControlPointCoordinates[cpIndex, 3]
					});
					cpCoordinates.Scale(pointFunctionsKsi[i]);
					temp = new Vector(temp + cpCoordinates);
				}

				temp.Scale(pointFunctionsHeta[j]);
				cartesianPoint = new Vector(cartesianPoint + temp);
			}

			return cartesianPoint;
		}
	}
}
