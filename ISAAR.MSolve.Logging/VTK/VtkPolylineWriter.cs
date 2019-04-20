using System;
using System.Collections.Generic;
using System.IO;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.XFEM.Output.VTK
{
    public class VtkPolylineWriter : IDisposable
    {
        public const string vtkReaderVersion = "4.1";
        private readonly StreamWriter writer;

        public VtkPolylineWriter(string filePath)
        {
            this.writer = new StreamWriter(filePath);
            writer.Write("# vtk DataFile Version ");
            writer.WriteLine(vtkReaderVersion);
            writer.WriteLine(filePath);
            writer.Write("ASCII\n\n");
        }

        public void Dispose()
        {
            if (writer != null) writer.Dispose();
        }

        public void WritePolyline(IReadOnlyList<CartesianPoint> vertices)
        {
            // Vertices 
            writer.WriteLine("DATASET POLYDATA");
            writer.WriteLine($"POINTS {vertices.Count} double");
            for (int i = 0; i < vertices.Count; ++i)
            {
                var vertex = vertices[i];
                writer.WriteLine($"{vertex.X} {vertex.Y} 0.0");
            }
            writer.WriteLine();

            // Line segments
            int numLines = vertices.Count - 1;
            writer.WriteLine($"LINES {numLines} {3* numLines}");
            for (int i = 0; i < numLines; ++i)
            {
                writer.WriteLine($"2 {i} {i+1}"); 
            }
            writer.WriteLine();
        }

        public void WritePolylineAsUnstructured(IReadOnlyList<CartesianPoint> vertices)
        {
            // Nodes 
            writer.WriteLine("DATASET UNSTRUCTURED_GRID");
            writer.WriteLine($"POINTS {vertices.Count} double");
            for (int i = 0; i < vertices.Count; ++i)
            {
                var vertex = vertices[i];
                writer.WriteLine($"{vertex.X} {vertex.Y} 0.0");
            }
            writer.WriteLine();

            // Cell connectivity
            int numPolylines = 1;
            int cellDataCount = vertices.Count + 1;
            writer.WriteLine($"\nCELLS {numPolylines} {cellDataCount}");

            writer.Write(vertices.Count);
            for (int i = 0; i < vertices.Count; ++i)
            {
                writer.Write($" {i}");
            }
            writer.WriteLine();

            // Cell types
            writer.WriteLine("\nCELL_TYPES " + numPolylines);
            writer.WriteLine("4");
        }
    }
}
