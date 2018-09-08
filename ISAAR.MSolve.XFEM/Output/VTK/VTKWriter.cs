using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Tensors;

//TODO: This should be IDisposable
//TODO: It should work with Mesh, not Model. What about data that is interpolated over specific regions (e.g. narrow band LSM)?
//TODO: Passing dictionaries might be preferrable to passing arrays, since the user does not have to guess the entries order and
//      he doesn't have to create the array from the dictionary he has. Even better, provide both functionalitites.
namespace ISAAR.MSolve.XFEM.Output.VTK
{
    class VTKWriter
    {
        public static string vtkReaderVersion = "4.1";
        private static readonly string directory = 
            Directory.GetParent(Directory.GetCurrentDirectory()).Parent.FullName + "\\Resources\\";
        private static readonly Dictionary<IsoparametricElementType2D, int> cellTypeCodes = 
            new Dictionary<IsoparametricElementType2D, int>()
            {
                { IsoparametricElementType2D.Quad4, 9 }
            };

        private StreamWriter writer;
        private readonly Model2D model;

        public VTKWriter(Model2D model)
        {
            this.model = model;
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

            WriteMesh();

            // Fields
            writer.Write("\n\n");
            writer.WriteLine("POINT_DATA " + model.Nodes.Count);
        }

        private void WriteMesh()
        {
            // Nodes 
            writer.WriteLine("DATASET UNSTRUCTURED_GRID");
            writer.WriteLine($"POINTS {model.Nodes.Count} double");
            for (int n = 0; n < model.Nodes.Count; ++n) // Their indices in Model.Nodes are equal to their IDs
            {
                XNode2D node = model.Nodes[n];
                writer.WriteLine($"{node.X} {node.Y} 0.0");
            }

            // Element connectivity
            int elementDataCount = 0;
            foreach (var element in model.Elements) elementDataCount += 1 + element.Nodes.Count;
            writer.WriteLine($"\nCELLS {model.Elements.Count} {elementDataCount}");
            foreach (var element in model.Elements)
            {
                writer.Write(element.Nodes.Count); 
                foreach (var node in element.Nodes)
                {
                    writer.Write(' ');
                    writer.Write(node.ID);
                }
                writer.WriteLine();
            }

            // Element types
            writer.WriteLine("\nCELL_TYPES " + model.Elements.Count);
            foreach (var element in model.Elements)
            {
                writer.WriteLine(cellTypeCodes[element.ElementType]);
            }
        }

        public void WriteScalarField(string fieldName, double[] nodalValues)
        {
            writer.WriteLine($"SCALARS {fieldName} double 1");
            writer.WriteLine("LOOKUP_TABLE default");
            for (int i = 0; i < model.Nodes.Count; ++i)
            {
                writer.WriteLine(nodalValues[i]);
            }
            writer.WriteLine();
        }

        public void WriteVector2DField(string fieldName, double[,] nodalValues)
        {
            writer.WriteLine($"VECTORS {fieldName} double");
            for (int i = 0; i < model.Nodes.Count; ++i)
            {
                writer.WriteLine($"{nodalValues[i, 0]} {nodalValues[i, 1]} 0.0");
            }
            writer.WriteLine();
        }
        
        public void WriteTensor2DField(string fieldName, IReadOnlyList<Tensor2D> nodalTensors)
        {
            writer.WriteLine($"TENSORS {fieldName} double");
            for (int i = 0; i < model.Nodes.Count; ++i)
            {
                writer.WriteLine($"{nodalTensors[i].XX} {nodalTensors[i].YY} {nodalTensors[i].XY} 0.0 0.0 0.0 0.0 0.0 0.0");
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
