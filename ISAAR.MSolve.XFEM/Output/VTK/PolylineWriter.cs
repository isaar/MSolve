using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Output.VTK
{
    class PolylineWriter
    {
        public static string vtkReaderVersion = "4.1";
        private static readonly string directory =
            Directory.GetParent(Directory.GetCurrentDirectory()).Parent.FullName + "\\Resources\\";
        private StreamWriter writer;

        public PolylineWriter()
        {
        }

        public void CloseCurrentFile()
        {
            writer.Flush();
            writer.Close();
        }

        public void InitializeFile(string filename, bool isAbsolute)
        {
            // Header
            string path = isAbsolute ? filename + ".vtk" : directory + filename + ".vtk";
            writer = new StreamWriter(path);
            writer.Write("# vtk DataFile Version ");
            writer.WriteLine(vtkReaderVersion);
            writer.WriteLine(filename);
            writer.Write("ASCII\n\n");
        }

        public void WritePolyline(IReadOnlyList<ICartesianPoint2D> vertices)
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

        public void WritePolylineAsUnstructured(IReadOnlyList<ICartesianPoint2D> vertices)
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
