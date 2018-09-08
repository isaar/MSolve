using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Output.VTK
{
    class VTKPointWriter
    {
        public static string vtkReaderVersion = "4.1";
        private static readonly string directory =
            Directory.GetParent(Directory.GetCurrentDirectory()).Parent.FullName + "\\Resources\\";

        private StreamWriter writer;

        public VTKPointWriter()
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

        public void WriteScalarField(string fieldName, IReadOnlyDictionary<ICartesianPoint2D, double> pointValues)
        {
            // Points
            writer.WriteLine("DATASET UNSTRUCTURED_GRID");
            writer.WriteLine($"POINTS {pointValues.Count} double");
            foreach (var point in pointValues.Keys)
            {
                writer.WriteLine($"{point.X} {point.Y} 0.0");
            }

            // Values
            writer.Write("\n\n");
            writer.WriteLine($"POINT_DATA {pointValues.Count}");
            writer.WriteLine($"SCALARS {fieldName} double 1");
            writer.WriteLine("LOOKUP_TABLE default");
            foreach (var value in pointValues.Values)
            {
                writer.WriteLine(value);
            }
            writer.WriteLine();
        }

        
    }
}
