using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Tensors;

namespace ISAAR.MSolve.XFEM.Output.VTK
{
    class DiscontinuousMeshVTKWriter
    {
        public static string vtkReaderVersion = "4.1";
        private static readonly string directory =
            Directory.GetParent(Directory.GetCurrentDirectory()).Parent.FullName + "\\Resources\\";
        

        private readonly Model2D originalModel;
        private readonly VTKMesh2D auxiliaryMesh;
        private StreamWriter writer;

        public DiscontinuousMeshVTKWriter(Model2D model)
        {
            this.originalModel = model;
            this.auxiliaryMesh = CreateAuxiliaryMesh(model);
        }

        private static VTKMesh2D CreateAuxiliaryMesh(Model2D originalModel)
        {
            var allCells = new List<VTKCell>();
            var allPoints = new List<VtkPoint2D>();
            int pointID = 0;

            foreach (var element in originalModel.Elements)
            {
                var cellVertices = new VtkPoint2D[element.Nodes.Count];
                for (int i = 0; i < element.Nodes.Count; ++i)
                {
                    var point = new VtkPoint2D(pointID++, element.Nodes[i]);
                    cellVertices[i] = point;
                    allPoints.Add(point);
                }
                int cellType = VTKCell.CellTypeCodes[element.ElementType];
                allCells.Add(new VTKCell(element, cellType, cellVertices));
            }

            return new VTKMesh2D(allPoints, allCells);
        }

        public void InitializeFile(string filename)
        {
            // Header
            string path = directory + filename + ".vtk";
            writer = new StreamWriter(path);
            writer.Write("# vtk DataFile Version ");
            writer.WriteLine(vtkReaderVersion);
            writer.WriteLine(filename);
            writer.Write("ASCII\n\n");

            WriteMesh();

            // Fields
            writer.Write("\n\n");
            writer.WriteLine("POINT_DATA " + auxiliaryMesh.Points.Count);
        }

        private void WriteMesh()
        {
            // Nodes 
            writer.WriteLine("DATASET UNSTRUCTURED_GRID");
            writer.WriteLine(String.Format("POINTS {0} double", auxiliaryMesh.Points.Count));
            for (int i = 0; i < auxiliaryMesh.Points.Count; ++i) // Their indices in Model.Nodes are equal to their IDs
            {
                VtkPoint2D point = auxiliaryMesh.Points[i];
                writer.Write(String.Format("{0} {1} 0.0\n", point.X, point.Y));
            }

            // Element connectivity
            int elementDataCount = 0;
            foreach (var cell in auxiliaryMesh.Cells) elementDataCount += 1 + cell.Vertices.Count;
            writer.WriteLine(String.Format("\nCELLS {0} {1}", auxiliaryMesh.Cells.Count, elementDataCount));
            foreach (var cell in auxiliaryMesh.Cells)
            {
                writer.Write(cell.Vertices.Count);
                foreach (var point in cell.Vertices)
                {
                    writer.Write(' ');
                    writer.Write(point.ID);
                }
                writer.WriteLine();
            }

            // Element types
            writer.WriteLine("\nCELL_TYPES " + auxiliaryMesh.Cells.Count);
            foreach (var cell in auxiliaryMesh.Cells)
            {
                writer.WriteLine(cell.TypeCode);
            }
        }

        public void WriteScalarField(string fieldName,
            IReadOnlyDictionary<XContinuumElement2D, IReadOnlyList<double>> valuesAtElementNodes)
        {
            var nodalValues = new double[auxiliaryMesh.Points.Count];
            foreach (var cell in auxiliaryMesh.Cells)
            {
                IReadOnlyList<double> valuesAtCellVertices = valuesAtElementNodes[cell.OriginalElement];
                for (int i = 0; i < cell.Vertices.Count; ++i)
                {
                    nodalValues[cell.Vertices[i].ID] = valuesAtCellVertices[i];
                }
            }

            writer.WriteLine(String.Format("SCALARS {0} double 1", fieldName));
            writer.WriteLine("LOOKUP_TABLE default");
            for (int i = 0; i < auxiliaryMesh.Points.Count; ++i)
            {
                writer.WriteLine(nodalValues[i]);
            }
            writer.WriteLine();
        }

        public void WriteVector2DField(string fieldName,
            IReadOnlyDictionary<XContinuumElement2D, IReadOnlyList<Vector2>> vectorsAtElementNodes)
        {
            var nodalVectors = new double[auxiliaryMesh.Points.Count, 2];
            foreach (var cell in auxiliaryMesh.Cells)
            {
                IReadOnlyList<Vector2> valuesAtCellVertices = vectorsAtElementNodes[cell.OriginalElement];
                for (int i = 0; i < cell.Vertices.Count; ++i)
                {
                    nodalVectors[cell.Vertices[i].ID, 0] = valuesAtCellVertices[i][0];
                    nodalVectors[cell.Vertices[i].ID, 1] = valuesAtCellVertices[i][1];
                }
            }

            writer.WriteLine(String.Format("VECTORS {0} double", fieldName));
            for (int i = 0; i < auxiliaryMesh.Points.Count; ++i)
            {
                writer.WriteLine(String.Format("{0} {1} 0.0", nodalVectors[i, 0], nodalVectors[i, 1]));
            }
            writer.WriteLine();
        }

        public void WriteTensor2DField(string fieldName, 
            IReadOnlyDictionary<XContinuumElement2D, IReadOnlyList<Tensor2D>> tensorsAtElementNodes)
        {
            var nodalTensors = new Tensor2D[auxiliaryMesh.Points.Count];
            foreach (var cell in auxiliaryMesh.Cells)
            {
                IReadOnlyList<Tensor2D> valuesAtCellVertices = tensorsAtElementNodes[cell.OriginalElement];
                for (int i = 0; i < cell.Vertices.Count; ++i)
                {
                    nodalTensors[cell.Vertices[i].ID] = valuesAtCellVertices[i];
                }
            }

            writer.WriteLine(String.Format("TENSORS {0} double", fieldName));
            for (int i = 0; i < auxiliaryMesh.Points.Count; ++i)
            {
                writer.WriteLine(String.Format("{0} {1} {2} 0.0 0.0 0.0 0.0 0.0 0.0",
                    nodalTensors[i].XX, nodalTensors[i].YY, nodalTensors[i].XY));
            }
            writer.WriteLine();
        }

        public void CloseCurrentFile()
        {
            writer.Flush();
            writer.Close();
        }
    }
}
