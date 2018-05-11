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
            for (int i = 0; i < vertices.Count; ++i) // Their indices in Model.Nodes are equal to their IDs
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
    }
}
