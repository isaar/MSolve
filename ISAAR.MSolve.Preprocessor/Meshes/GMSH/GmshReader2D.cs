using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using ISAAR.MSolve.FEM.Entities;

namespace ISAAR.MSolve.Preprocessor.Meshes.GMSH
{
    /// <summary>
    /// Creates meshes by reading GMSH output files (.msh). Unrecognized GMSH cell types will be ignored along with any 1D cells
    /// in the .msh file, therefore some care is needed.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class GmshReader2D: IDisposable, IMeshProvider2D<Node2D>
    {
        private readonly StreamReader reader;

        /// <summary>
        /// Opens the .msh file but doesn't read it.
        /// </summary>
        /// <param name="mshFilePath">The absolute path of the .msh file where GMSH has written the mesh data. The .msh file 
        ///     will not be modified.</param>
        public GmshReader2D(string mshFilePath)
        {
            reader = new StreamReader(mshFilePath);
        }

        /// <summary>
        /// Reads the whole .msh file and converts it to MSolve mesh data.
        /// </summary>
        /// <returns></returns>
        public (IReadOnlyList<Node2D> vertices, IReadOnlyList<CellConnectivity2D> cells) CreateMesh()
        {
            // Vertices must be listed before cells
            Node2D[] vertices = ReadVertices();
            IReadOnlyList<CellConnectivity2D> cells = ReadCells(vertices);
            return (vertices, cells);
        }

        public void Dispose()
        {
            if (reader != null) reader.Dispose();
        }

        private Node2D[] ReadVertices()
        {
            string line;

            // Find node segment
            while (true)
            {
                line = reader.ReadLine();
                if (line[0] == '$')
                    if (line.Equals("$Nodes")) break; // Next line will be the nodes count.
            }

            // Read the vertices
            int numVertices = int.Parse(reader.ReadLine());
            Node2D[] vertices = new Node2D[numVertices];
            while (true)
            {
                line = reader.ReadLine();
                if (line[0] == '$') break; // This line is "$EndNodes". Next line will be the next segment.
                else
                {
                    // Line = nodeID x y z
                    string[] words = line.Split(new char[] { ' ' });
                    int id = int.Parse(words[0]) - 1; // MSolve uses 0-based indexing
                    double x = double.Parse(words[1]);
                    double y = double.Parse(words[2]);
                    vertices[id] = new Node2D(id, x, y);
                }
            }

            return vertices;
        }

        // It must be called after vertices are read.
        private IReadOnlyList<CellConnectivity2D> ReadCells(Node2D[] vertices)
        {
            string line;

            // Find cell segment
            while (true)
            {
                line = reader.ReadLine();
                if (line[0] == '$')
                    if (line.Equals("$Elements")) break; // Next line will be the number of cells.
            }

            // Read the elements
            int numFauxCells = int.Parse(reader.ReadLine()); // not all of them are actual 2D cells
            var cells = new List<CellConnectivity2D>();
            var cellFactory = new GmshCell2DFactory(vertices);
            while (true)
            {
                line = reader.ReadLine();
                if (line[0] == '$') break; // This line is "$EndElements". Next line will be the next segment.
                else
                {
                    // Line = cellID cellCode tagsCount <tags>(0 or more) vertexIds(2 or more)
                    string[] words = line.Split(new char[] { ' ' });
                    int code = int.Parse(words[1]);
                    int numTags = int.Parse(words[2]);

                    int firstVertexPos = 3 + numTags;
                    int numVertices = words.Length - firstVertexPos;
                    int[] vertexIDs = new int[numVertices];
                    for (int i = 0; i < numVertices; ++i)
                    {
                        vertexIDs[i] = int.Parse(words[firstVertexPos + i]) - 1; // MSolve uses 0-based indexing
                    }

                    bool validCell = cellFactory.TryCreateCell(code, vertexIDs, out CellConnectivity2D cell);
                    if (validCell) cells.Add(cell);
                }
            }

            return cells;
        }
    }
}
