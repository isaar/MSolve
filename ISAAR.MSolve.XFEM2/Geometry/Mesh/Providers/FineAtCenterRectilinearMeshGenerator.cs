using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Geometry.Mesh.SuiteSparse
{
    class FineAtCenterRectilinearMeshGenerator
    {
        public double[] domainLowerBounds, domainUpperBounds;
        public double fineElementSize;
        public int coarseElementCountPerRegion;
        public double[] interestAreaDimensions;

        public FineAtCenterRectilinearMeshGenerator()
        { }

        public (XNode2D[] nodes, List<XNode2D[]> elementConnectivity) CreateMesh()
        {
            double[] coordinatesX = FindNodalCoordinates(0);
            double[] coordinatesY = FindNodalCoordinates(1);
            var baseGenerator = new RectilinearMeshGenerator(coordinatesX, coordinatesY);
            return baseGenerator.CreateMesh();
        }

        public double[] FindNodalCoordinates(int dimension)
        {
            double p1 = domainLowerBounds[dimension];
            double p4 = domainUpperBounds[dimension];
            double center = 0.5 * (p1 + p4);

            // The mesh is finer in a region ~2.5*interestAreaDimension around the center
            double p2 = center - 1.25 * interestAreaDimensions[dimension];
            double p3 = center + 1.25 * interestAreaDimensions[dimension];

            // Adjust the fine region bounds
            int fineElementsCount = (int)(Math.Ceiling((p3 - p2) / fineElementSize));
            p2 = center - fineElementsCount * fineElementSize / 2.0;
            p3 = center + fineElementsCount * fineElementSize / 2.0;

            // The coordinate array
            int totalNodesCount = 2 * coarseElementCountPerRegion + fineElementsCount + 1;
            double[] coordinates = new double[totalNodesCount];

            // Region 1: coarse. Does not include the common node with region 2.
            int start1 = 0;
            int end1 = coarseElementCountPerRegion;
            double spacing1 = (p2 - p1) / coarseElementCountPerRegion;
            for (int i = start1; i < end1; ++i) coordinates[i] = i * spacing1;

            //Region 2: fine. Includes the common node with regions 1, but not the common node with region 3.
            int start2 = coarseElementCountPerRegion;
            int end2 = coarseElementCountPerRegion + fineElementsCount;
            double spacing2 = fineElementSize;
            for (int i = start2; i < end2; ++i) coordinates[i] = p2 + (i - start2) * spacing2;

            // Region 3: coarse. Includes the common node with region 2.
            int start3 = coarseElementCountPerRegion + fineElementsCount;
            int end3 = totalNodesCount;
            double spacing3 = (p4 - p3) / coarseElementCountPerRegion;
            for (int i = start3; i < end3; ++i) coordinates[i] = p3 + (i - start3) * spacing3;

            return coordinates;
        }
    }
}
