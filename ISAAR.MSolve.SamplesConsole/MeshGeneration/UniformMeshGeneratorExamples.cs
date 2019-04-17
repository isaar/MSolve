using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Mesh.Generation;
using ISAAR.MSolve.Discretization.Mesh.Generation.Custom;
using ISAAR.MSolve.FEM.Entities;

namespace ISAAR.MSolve.SamplesConsole.MeshGeneration
{
    public class UniformMeshGeneratorExamples
    {
        public static void Generate2DMesh()
        {
        }

        public static void Generate3DMesh()
        {
            double minX = 0.0, minY = 0.0, minZ = 0.0;
            double maxX = 100.0, maxY = 100.0, maxZ = 100.0;
            int cellsPerX = 4, cellsPerY = 4, cellsPerZ = 4;

            var generator = new UniformMeshGenerator3D<Node>(minX, minY, minZ, maxX, maxY, maxZ, cellsPerX, cellsPerY, cellsPerZ);
            generator.StartIDsAt0 = false;
            (IReadOnlyList<Node> vertices, IReadOnlyList<CellConnectivity<Node>> cells) = 
                generator.CreateMesh((id, x, y, z) => new Node(id: id, x: x, y:  y, z: z ));

            Console.WriteLine($"{vertices.Count} vertices:");
            foreach (var vertex in vertices)
            {
                Console.WriteLine($"{vertex.ID} {vertex.X} {vertex.Y} {vertex.Z}");
            }

            Console.WriteLine();
            Console.WriteLine($"{cells.Count} cells:");
            int cellID = 1;
            string gmshHexa8Tag = "5 2 1 1";
            foreach (var cell in cells)
            {
                Console.Write($"{cellID++} {gmshHexa8Tag} ");
                foreach (var vertex in cell.Vertices) Console.Write($"{vertex.ID} ");
                Console.WriteLine();
            }
        }
    }
}
