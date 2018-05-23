using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Geometry.Mesh.Gmsh
{
    // Nodes must be listed before elements
    class GmshReader
    {
        private static readonly string directory =
            Directory.GetParent(Directory.GetCurrentDirectory()).Parent.FullName + "\\Resources\\";

        private readonly StreamReader reader;

        public GmshReader(string filename, bool absolute = false)
        {
            if (absolute) reader = new StreamReader(filename);
            else
            {
                string path = directory + filename;
                reader = new StreamReader(path);
            }
        }

        ~GmshReader()
        {
            reader.Close();
        }

        public Tuple<IReadOnlyList<XNode2D>, IReadOnlyList<GmshElement>> ReadMesh()
        {
            XNode2D[] nodes = ReadNodes();
            IReadOnlyList<GmshElement> elements = ReadElements(nodes);
            return new Tuple<IReadOnlyList<XNode2D>, IReadOnlyList<GmshElement>>(nodes, elements);
        }

        private XNode2D[] ReadNodes()
        {
            string line;

            // Find node segment
            while (true)
            {
                line = reader.ReadLine();
                if (line[0] == '$')
                    if (line.Equals("$Nodes")) break; // Next line will be the nodes count.
            }

            // Read the nodes
            int nodesCount = int.Parse(reader.ReadLine());
            XNode2D[] nodes = new XNode2D[nodesCount];
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
                    nodes[id] = new XNode2D(id, x, y);
                }
            }

            return nodes;
        }

        // It must be called after nodes are read.
        private IReadOnlyList<GmshElement> ReadElements(XNode2D[] nodes)
        {
            string line;

            // Find element segment
            while (true)
            {
                line = reader.ReadLine();
                if (line[0] == '$')
                    if (line.Equals("$Elements")) break; // Next line will be the elements count.
            }

            // Read the elements
            int fauxElementsCount = int.Parse(reader.ReadLine());
            var elements = new List<GmshElement>();
            var elementFactory = new GmshElementFactory(nodes);
            while (true)
            {
                line = reader.ReadLine();
                if (line[0] == '$') break; // This line is "$EndElements". Next line will be the next segment.
                else
                {
                    // Line = elementID elementCode tagsCount <tags>(0 or more) nodeIds(2 or more)
                    string[] words = line.Split(new char[] { ' ' });
                    int code = int.Parse(words[1]);
                    int tagsCount = int.Parse(words[2]);

                    int firstNodePos = 3 + tagsCount;
                    int nodesCount = words.Length - firstNodePos;
                    int[] nodeIDs = new int[nodesCount];
                    for (int i = 0; i < nodesCount; ++i)
                    {
                        nodeIDs[i] = int.Parse(words[firstNodePos + i]) - 1; // MSolve uses 0-based indexing
                    }

                    GmshElement element;
                    bool validElement = elementFactory.TryCreateElement(code, nodeIDs, out element);
                    if (validElement) elements.Add(element);
                }
            }

            return elements;
        }
            
    }
}
