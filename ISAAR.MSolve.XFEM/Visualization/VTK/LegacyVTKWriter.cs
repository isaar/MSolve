using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Visualization.VTK
{
    class LegacyVTKWriter
    {
        public static string vtkReaderVersion = "4.1";
        private static readonly string directory = 
            Directory.GetParent(Directory.GetCurrentDirectory()).Parent.FullName + "\\Output\\";
        private static readonly Dictionary<IsoparametricElementType2D, int> cellTypeCodes = 
            new Dictionary<IsoparametricElementType2D, int>()
            {
                { IsoparametricElementType2D.Quad4, 9 }
            };

        private StreamWriter writer;
        public Model2D Model { get; set; }

        public void InitializeFile(string filename)
        {
            string path = directory + filename + ".vtk";
            writer = new StreamWriter(path);
            writer.Write("# vtk DataFile Version ");
            writer.WriteLine(vtkReaderVersion);
            writer.WriteLine(filename);
            writer.Write("ASCII\n\n");
        }

        public void WriteMesh()
        {
            // Nodes 
            writer.WriteLine("DATASET UNSTRUCTURED_GRID");
            writer.WriteLine(String.Format("POINTS {0} double", Model.Nodes.Count));
            for (int n = 0; n < Model.Nodes.Count; ++n) // Their indices in Model.Nodes are equal to their IDs
            {
                XNode2D node = Model.Nodes[n];
                writer.Write(String.Format("{0} {1} 0.0\n", node.X, node.Y));
            }
            writer.WriteLine();

            // Element connectivity
            int elementDataCount = 0;
            foreach (var element in Model.Elements) elementDataCount += 1 + element.Nodes.Count;
            writer.WriteLine(String.Format("CELLS {0} {1}", Model.Elements.Count, elementDataCount));
            foreach (var element in Model.Elements)
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
            writer.WriteLine("CELL_TYPES " + Model.Elements.Count);
            foreach (var element in Model.Elements)
            {
                writer.WriteLine(cellTypeCodes[element.ElementType]);
            }
            writer.Write("\n\n");
        }

        public void WriteScalarField(string fieldName, double[] nodalValues)
        {
            writer.WriteLine("POINT_DATA " + Model.Nodes.Count);
            writer.WriteLine(String.Format("SCALARS {0} float 1", fieldName));
            writer.WriteLine("LOOKUP_TABLE default");
            for (int n = 0; n < Model.Nodes.Count; ++n)
            {
                writer.WriteLine(nodalValues[n]);
            }
            writer.WriteLine();
        }

        public void CloseCurrentFile()
        {
            writer.Flush();
        }
    }
}
