using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.Mesh.Gmsh;
using ISAAR.MSolve.XFEM.Geometry.Mesh.Providers;

namespace ISAAR.MSolve.XFEM.Tests.GRACM
{
    interface IMeshProvider
    {
        (XNode2D[] nodes, List<XNode2D[]> elementConnectivity) CreateMesh();
    }

    class DCBUniformMeshProvider : IMeshProvider
    {
        private readonly double elementSize;

        /// <summary>
        /// Creates a uniform rectilinear mesh.
        /// </summary>
        /// <param name="elementSize">The approximate size of the elements.</param>
        public DCBUniformMeshProvider(double elementSize)
        {
            this.elementSize = elementSize;
        }

        public (XNode2D[] nodes, List<XNode2D[]> elementConnectivity) CreateMesh()
        {
            int elementsPerXAxis = (int)(DCB.L / elementSize) + 1;
            int elementsPerYAxis = (int)(DCB.h / elementSize) + 1;
            var meshGenerator = new UniformRectilinearMeshGenerator(DCB.L, DCB.h, elementsPerXAxis, elementsPerYAxis);
            return meshGenerator.CreateMesh();
        }
    }

    class DCBRefinedMeshProvider : IMeshProvider
    {
        private readonly double coarseElementSize;
        private readonly double fineElementSize;

        /// <summary>
        /// Creates a rectilinear mesh that is refined around the expected crack path.
        /// </summary>
        /// <param name="fineElementSize">The approximate size of the elements in the fine region.</param>
        /// <param name="coarseElementSize">The approximate size of the elements in the coarse region.</param>
        public DCBRefinedMeshProvider(double fineElementSize, double coarseElementSize)
        {
            this.fineElementSize = fineElementSize;
            this.coarseElementSize = coarseElementSize;
        }

        public (XNode2D[] nodes, List<XNode2D[]> elementConnectivity) CreateMesh()
        {
            double[,] meshSizeAlongX = new double[,]
            {
                { 0.0, coarseElementSize },
                { DCB.a, fineElementSize},
                //{ 3.0 / 5.0* L, fineElementSize},
                { 6.4, fineElementSize},
                { DCB.L, coarseElementSize}
            };
            double[,] meshSizeAlongY = new double[,]
            {
                { 0.0, fineElementSize },
                //{ 0.0, fineElementSize},
                //{ 0.5, fineElementSize},
                { DCB.h/2 + DCB.h/10, fineElementSize}, //TODO: I think this was coarse
                { DCB.h, coarseElementSize}
            };

            var meshGenerator = new RectilinearMeshGenerator(meshSizeAlongX, meshSizeAlongY);
            return meshGenerator.CreateMesh();
        }
    }

    class GmshMeshProvider : IMeshProvider
    {
        private readonly string meshFilePath;

        public GmshMeshProvider(string meshFilePath)
        {
            this.meshFilePath = meshFilePath;
        }

        public (XNode2D[] nodes, List<XNode2D[]> elementConnectivity) CreateMesh()
        {
            var meshReader = new GmshReader(meshFilePath);
            Tuple<IReadOnlyList<XNode2D>, IReadOnlyList<GmshElement>> meshEntities = meshReader.ReadMesh();

            //TODO: standardize the output of all mesh generators, readers, etc
            var nodes = new XNode2D[meshEntities.Item1.Count];
            var elementConnectivity = new List<XNode2D[]>();
            nodes = meshEntities.Item1.ToArray();
            foreach (var element in meshEntities.Item2)
            {
                elementConnectivity.Add(element.Nodes.ToArray());
            }

            return (nodes, elementConnectivity);
            
        }
    }
}
