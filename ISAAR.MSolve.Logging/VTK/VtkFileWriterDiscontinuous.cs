using System;
using System.Collections.Generic;
using System.IO;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging.VTK;

//TODO: Extend this to FEM, but only if all elements are continuum (or at least 2D) without embedding or other fancy stuff.
namespace ISAAR.MSolve.XFEM.Output.VTK
{
    public class VtkFileWriterDiscontinuous<TNode> : IDisposable where TNode : INode
    {
        public static string vtkReaderVersion = "4.1";
        private readonly VtkMeshDiscontinuous<TNode> mesh;
        private readonly StreamWriter writer;
        private bool writeFieldsNext;

        public VtkFileWriterDiscontinuous(VtkMeshDiscontinuous<TNode> mesh, string filePath)
        {
            this.mesh = mesh;

            this.writer = new StreamWriter(filePath);
            writer.Write("# vtk DataFile Version ");
            writer.WriteLine(vtkReaderVersion);
            writer.WriteLine(filePath);
            writer.Write("ASCII\n\n");
            writeFieldsNext = false;
        }

        public void Dispose()
        {
            if (writer != null) writer.Dispose();
        }

        public void WriteMesh()
        {
            if (writeFieldsNext) throw new InvalidOperationException("A mesh has already been written.");

            // Nodes 
            writer.WriteLine("DATASET UNSTRUCTURED_GRID");
            writer.WriteLine(String.Format("POINTS {0} double", mesh.VtkPoints.Count));
            for (int i = 0; i < mesh.VtkPoints.Count; ++i) // Their indices in Model.Nodes are equal to their IDs
            {
                VtkPoint point = mesh.VtkPoints[i];
                writer.Write(String.Format("{0} {1} 0.0\n", point.X, point.Y));
            }

            // Element connectivity
            int elementDataCount = 0;
            foreach (VtkCell cell in mesh.VtkCells) elementDataCount += 1 + cell.Vertices.Count;
            writer.WriteLine(String.Format("\nCELLS {0} {1}", mesh.VtkCells.Count, elementDataCount));
            foreach (VtkCell cell in mesh.VtkCells)
            {
                writer.Write(cell.Vertices.Count);
                foreach (VtkPoint point in cell.Vertices)
                {
                    writer.Write(' ');
                    writer.Write(point.ID);
                }
                writer.WriteLine();
            }

            // Element types
            writer.WriteLine("\nCELL_TYPES " + mesh.VtkCells.Count);
            foreach (VtkCell cell in mesh.VtkCells) writer.WriteLine(cell.Code);
        }

        public void WriteScalarField(string fieldName,
            IReadOnlyDictionary<ICell<TNode>, IReadOnlyList<double>> valuesAtElementNodes)
        {
            var nodalValues = new double[mesh.VtkPoints.Count];
            for (int e = 0; e < mesh.VtkCells.Count; ++e)
            {
                VtkCell cell = mesh.VtkCells[e];
                IReadOnlyList<double> valuesAtCellVertices = valuesAtElementNodes[mesh.OriginalElements[e]];
                for (int i = 0; i < cell.Vertices.Count; ++i)
                {
                    nodalValues[cell.Vertices[i].ID] = valuesAtCellVertices[i];
                }
            }

            WriteFieldsHeader(nodalValues.Length);
            writer.WriteLine(String.Format("SCALARS {0} double 1", fieldName));
            writer.WriteLine("LOOKUP_TABLE default");
            for (int i = 0; i < mesh.VtkPoints.Count; ++i)
            {
                writer.WriteLine(nodalValues[i]);
            }
            writer.WriteLine();
        }

        public void WriteTensor2DField(string fieldName,
            IReadOnlyDictionary<ICell<TNode>, IReadOnlyList<Tensor2D>> tensorsAtElementNodes)
        {
            var nodalTensors = new Tensor2D[mesh.VtkPoints.Count];
            for (int e = 0; e < mesh.VtkCells.Count; ++e)
            {
                VtkCell cell = mesh.VtkCells[e];
                IReadOnlyList<Tensor2D> valuesAtCellVertices = tensorsAtElementNodes[mesh.OriginalElements[e]];
                for (int i = 0; i < cell.Vertices.Count; ++i)
                {
                    nodalTensors[cell.Vertices[i].ID] = valuesAtCellVertices[i];
                }
            }

            WriteFieldsHeader(nodalTensors.Length);
            writer.WriteLine(String.Format("TENSORS {0} double", fieldName));
            for (int i = 0; i < mesh.VtkPoints.Count; ++i)
            {
                writer.WriteLine(String.Format("{0} {1} {2} 0.0 0.0 0.0 0.0 0.0 0.0",
                    nodalTensors[i].XX, nodalTensors[i].YY, nodalTensors[i].XY));
            }
            writer.WriteLine();
        }

        public void WriteVector2DField(string fieldName,
            IReadOnlyDictionary<ICell<TNode>, IReadOnlyList<Vector2>> vectorsAtElementNodes)
        {
            var nodalVectors = new double[mesh.VtkPoints.Count, 2];
            for (int e = 0; e < mesh.VtkCells.Count; ++e)
            {
                VtkCell cell = mesh.VtkCells[e];
                IReadOnlyList<Vector2> valuesAtCellVertices = vectorsAtElementNodes[mesh.OriginalElements[e]];
                for (int i = 0; i < cell.Vertices.Count; ++i)
                {
                    nodalVectors[cell.Vertices[i].ID, 0] = valuesAtCellVertices[i][0];
                    nodalVectors[cell.Vertices[i].ID, 1] = valuesAtCellVertices[i][1];
                }
            }

            WriteFieldsHeader(nodalVectors.Length);
            writer.WriteLine(String.Format("VECTORS {0} double", fieldName));
            for (int i = 0; i < mesh.VtkPoints.Count; ++i)
            {
                writer.WriteLine(String.Format("{0} {1} 0.0", nodalVectors[i, 0], nodalVectors[i, 1]));
            }
            writer.WriteLine();
        }

        private void WriteFieldsHeader(int numPoints)
        {
            if (!writeFieldsNext) // Fields header
            {
                writer.Write("\n\n");
                writer.WriteLine("POINT_DATA " + numPoints);
                writeFieldsNext = true;
            }
        }
    }
}
