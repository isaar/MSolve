using System;
using System.Collections.Generic;
using System.IO;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.XFEM.Output.VTK
{
    public class VtkPointWriter : IDisposable
    {
        public const string vtkReaderVersion = "4.1";
        private readonly StreamWriter writer;

        public VtkPointWriter(string filePath)
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

        public void WriteScalarField(string fieldName, IReadOnlyDictionary<CartesianPoint, double> pointValues)
        {
            // Points
            writer.WriteLine("DATASET UNSTRUCTURED_GRID");
            writer.WriteLine($"POINTS {pointValues.Count} double");
            foreach (var point in pointValues.Keys)
            {
                writer.WriteLine($"{point.X} {point.Y} {point.Z}");
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

        //TODO: Avoid the duplication.
        public void WriteScalarField(string fieldName, IReadOnlyDictionary<INode, double> nodalValues)
        {
            // Points
            writer.WriteLine("DATASET UNSTRUCTURED_GRID");
            writer.WriteLine($"POINTS {nodalValues.Count} double");
            foreach (var point in nodalValues.Keys)
            {
                writer.WriteLine($"{point.X} {point.Y} {point.Z}");
            }

            // Values
            writer.Write("\n\n");
            writer.WriteLine($"POINT_DATA {nodalValues.Count}");
            writer.WriteLine($"SCALARS {fieldName} double 1");
            writer.WriteLine("LOOKUP_TABLE default");
            foreach (var value in nodalValues.Values)
            {
                writer.WriteLine(value);
            }
            writer.WriteLine();
        }

    }
}
