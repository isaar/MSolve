using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Solvers.Interfaces;

namespace ISAAR.MSolve.IGA.Postprocessing
{
	public class ParaviewTsplines2D
	{
		private Model _model;
		private ILinearSystem _linearSystem;
		private string _filename;

		public ParaviewTsplines2D(Model model, ILinearSystem linearSystem, string filename)
		{
			_model = model;
			_linearSystem = linearSystem;
			_filename = filename;
		}

		public void CreateParaviewFile()
		{
			var projectiveControlPoints = CalculateProjectiveControlPoints();
			var numberOfPointsPerElement = 4;
			var nodes = new double[_model.Elements.Count * numberOfPointsPerElement, 2];
			var pointIndex = 0;
			foreach (var element in _model.Elements)
			{
				var tsplineElement = element as TSplineKirchhoffLoveShellElement;
				var elementPoints = tsplineElement.CalculatePointsForPostProcessing(tsplineElement);
				for (int i = 0; i < elementPoints.GetLength(0); i++)
				{
					nodes[pointIndex, 0] = elementPoints[i, 0];
					nodes[pointIndex++, 1] = elementPoints[i, 1];
				}
			}

			var elementConnectivity = CreateTsplineConnectivity();

			var pointDisplacements = new double[nodes.GetLength(0), 2];
			foreach (var element in _model.Elements)
			{
				var localDisplacements = new double[element.ControlPoints.Count, 2];
				var counterCP = 0;
				foreach (var controlPoint in element.ControlPoints)
				{
					//var dofX = _model.ControlPointDOFsDictionary[controlPoint.ID][DOFType.X];
					//var dofY = _model.ControlPointDOFsDictionary[controlPoint.ID][DOFType.Y];

					var dofX = _model.GlobalDofOrdering.GlobalFreeDofs[controlPoint, DOFType.X];
					var dofY = _model.GlobalDofOrdering.GlobalFreeDofs[controlPoint, DOFType.Y];

					localDisplacements[counterCP, 0] = (dofX == -1) ? 0.0 : _linearSystem.Solution[dofX];
					localDisplacements[counterCP++, 1] = (dofY == -1) ? 0.0 : _linearSystem.Solution[dofY];
				}
				var elementKnotDisplacements = element.ElementType.CalculateDisplacementsForPostProcessing(element, localDisplacements);
				for (int i = 0; i < elementConnectivity.GetLength(1); i++)
				{
					var knotConnectivity = elementConnectivity[element.ID, i];
					pointDisplacements[knotConnectivity, 0] = elementKnotDisplacements[i, 0];
					pointDisplacements[knotConnectivity, 1] = elementKnotDisplacements[i, 1];
				}
			}

			Write2DTSplinesFile(nodes, elementConnectivity, pointDisplacements);
		}

		public void Write2DTSplinesFile(double[,] nodeCoordinates, int[,] elementConnectivity, double[,] displacements)
		{
			var numberOfPoints = nodeCoordinates.GetLength(0);
			var numberOfCells = elementConnectivity.GetLength(0);

			int numberOfVerticesPerCell = 0;
			int paraviewCellCode = 0;

			numberOfVerticesPerCell = 4;
			paraviewCellCode = 9;

			using (StreamWriter outputFile = new StreamWriter($"..\\..\\..\\OutputFiles\\{_filename}Paraview.vtu"))
			{
				outputFile.WriteLine("<?xml version=\"1.0\"?>");
				outputFile.WriteLine("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">");
				outputFile.WriteLine("<UnstructuredGrid>");

				outputFile.WriteLine($"<Piece NumberOfPoints=\"{numberOfPoints}\" NumberOfCells=\"{numberOfCells}\">");
				outputFile.WriteLine($"<PointData Vectors=\"U\">");
				outputFile.WriteLine($"<DataArray type=\"Float32\" Name=\"U\" format=\"ascii\" NumberOfComponents=\"3\">");

				for (int i = 0; i < numberOfPoints; i++)
					outputFile.WriteLine($"{displacements[i, 0]} {displacements[i, 1]} 0.0");

				outputFile.WriteLine("</DataArray>");
				outputFile.WriteLine("</PointData>");
				outputFile.WriteLine("<Points>");
				outputFile.WriteLine("<DataArray type=\"Float32\" NumberOfComponents=\"3\">");

				for (int i = 0; i < numberOfPoints; i++)
					outputFile.WriteLine($"{nodeCoordinates[i, 0]} {nodeCoordinates[i, 1]} 0.0");
				

				outputFile.WriteLine("</DataArray>");
				outputFile.WriteLine("</Points>");
				outputFile.WriteLine("<Cells>");
				outputFile.WriteLine("<DataArray type=\"Int32\" Name=\"connectivity\">");

				for (int i = 0; i < numberOfCells; i++)
				{
					for (int j = 0; j < elementConnectivity.GetLength(1); j++)
						outputFile.Write($"{elementConnectivity[i, j]} ");
					outputFile.WriteLine("");
				}

				outputFile.WriteLine("</DataArray>");
				outputFile.WriteLine("<DataArray type=\"Int32\" Name=\"offsets\">");

				var offset = 0;
				for (int i = 0; i < numberOfCells; i++)
				{
					offset += numberOfVerticesPerCell;
					outputFile.WriteLine(offset);
				}

				outputFile.WriteLine("</DataArray>");
				outputFile.WriteLine("<DataArray type=\"Int32\" Name=\"types\">");

				for (int i = 0; i < numberOfCells; i++)
					outputFile.WriteLine(paraviewCellCode);

				outputFile.WriteLine("</DataArray>");
				outputFile.WriteLine("</Cells>");
				outputFile.WriteLine("</Piece>");
				outputFile.WriteLine("</UnstructuredGrid>");
				outputFile.WriteLine("</VTKFile>");
			}
		}

		private int[,] CreateTsplineConnectivity()
		{
			var elementConnectivity = new int[_model.Elements.Count, 4];
			var pointCounter = 0;
			for (int i = 0; i < _model.Elements.Count; i++)
			{
				elementConnectivity[i, 0] = pointCounter++;
				elementConnectivity[i, 1] = pointCounter++;
				elementConnectivity[i, 2] = pointCounter++;
				elementConnectivity[i, 3] = pointCounter++;
			}

			return elementConnectivity;
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
