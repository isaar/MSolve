using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.FEM.Postprocessing
{
    public class ParaviewEmbedded3D
    {
        private Model_v2 _model;
        private IVectorView _solution;
        private string _filename;

        public ParaviewEmbedded3D(Model_v2 model, IVectorView solution, string filename)
        {
            _model = model;
            _solution = solution;
            _filename = filename;
        }

        public void CreateParaviewFile()
        {
            WriteParaviewFile3D();
        }

        //private int[] conn = new int[] { 1, 2, 4, 3, 5, 6, 8, 7 };
        private int[] conn = new int[] { 6, 7, 4, 5, 2, 3, 0, 1 };
        private void WriteParaviewFile3D()
        {
            var elements = _model.Elements.Where(e=>e.ElementType is Hexa8NonLinear).ToList();
            var nodes = new List<INode>();
            elements.ForEach(e => nodes.AddRange(e.Nodes));
            nodes = nodes.Distinct().ToList();

            var numberOfPoints = nodes.Count;
            var numberOfCells = elements.Count;

            int numberOfVerticesPerCell = 0;
            int paraviewCellCode = 0;

            numberOfVerticesPerCell = 8;

            paraviewCellCode = 13;

            using (StreamWriter outputFile= new StreamWriter($"E:\\GEORGE_DATA\\DESKTOP\\output files\\{_filename}Paraview.vtu"))
            {
                outputFile.WriteLine("<?xml version=\"1.0\"?>");
                outputFile.WriteLine("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">");
                outputFile.WriteLine("<UnstructuredGrid>");

                outputFile.WriteLine($"<Piece NumberOfPoints=\"{numberOfPoints}\" NumberOfCells=\"{numberOfCells}\">");
                outputFile.WriteLine($"<PointData Vectors=\"U\">");
                outputFile.WriteLine($"<DataArray type=\"Float32\" Name=\"U\" format=\"ascii\" NumberOfComponents=\"3\">");

                for (int i = 0; i < numberOfPoints; i++)
                {
                    var node = _model.Nodes[i];
                    var dx = (!_model.GlobalDofOrdering.GlobalFreeDofs.Contains(node, DOFType.X))
                            ? 0.0
                            : _solution[_model.GlobalDofOrdering.GlobalFreeDofs[node, DOFType.X]];
                    var dy = (!_model.GlobalDofOrdering.GlobalFreeDofs.Contains(node, DOFType.Y))
                            ? 0.0
                            : _solution[_model.GlobalDofOrdering.GlobalFreeDofs[node, DOFType.Y]];
                    var dz = (!_model.GlobalDofOrdering.GlobalFreeDofs.Contains(node, DOFType.Z))
                            ? 0.0
                            : _solution[_model.GlobalDofOrdering.GlobalFreeDofs[node, DOFType.Z]];

                    outputFile.WriteLine($"{dx} {dy} {dz}");
                }
                    


                outputFile.WriteLine("</DataArray>");

                outputFile.WriteLine("</PointData>");
                outputFile.WriteLine("<Points>");
                outputFile.WriteLine("<DataArray type=\"Float32\" NumberOfComponents=\"3\">");

                for (int i = 0; i < numberOfPoints; i++)
                    outputFile.WriteLine($"{_model.Nodes[i].X} {_model.Nodes[i].Y} {_model.Nodes[i].Z}");

                outputFile.WriteLine("</DataArray>");
                outputFile.WriteLine("</Points>");
                outputFile.WriteLine("<Cells>");
                outputFile.WriteLine("<DataArray type=\"Int32\" Name=\"connectivity\">");

                for (int i = 0; i < numberOfCells; i++)
                {
                    for (int j = 0; j < _model.Elements[i].Nodes.Count; j++)
                        outputFile.Write($"{_model.Elements[i].Nodes[conn[j]].ID} ");
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
    }
}
