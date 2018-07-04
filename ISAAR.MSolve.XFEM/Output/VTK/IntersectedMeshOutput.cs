using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Triangulation;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Interpolation.GaussPointSystems;
using ISAAR.MSolve.XFEM.Interpolation.InverseMappings;
using System;
using System.Collections.Generic;
using System.Text;

//TODO: the extrapolations from Gauss points to nodes of the subtriangles does not work near the tip, where the displacement 
//      field contains interpolated tip functions. This could be fixed (for the subtriangle nodes only) by basing the 
//      extrapolation on the XFEM interpolation instead of the standard Tri3.
//TODO: Can I cheat and assign H(x)= +-1 to points on the crack, depending on which side of the crack the triangle lies? That 
//      way I could avoid the need to extrapolate the displacements form Gauss points. Would it also suffice for stresses or do
//      I have to extrapolate?
namespace ISAAR.MSolve.XFEM.Output.VTK
{
    class IntersectedMeshOutput
    {
        private const int triangleVtkCode = 5;
        private readonly TrackingExteriorCrackLSM crack;
        private readonly Model2D model;
        private readonly string pathNoExtension;
        private readonly CartesianTriangulator triangulator;

        //private IReadOnlyList<VtkPoint2D> points;
        //private IReadOnlyList<VtkCell2D> cells;

        public IntersectedMeshOutput(Model2D model, TrackingExteriorCrackLSM crack, string pathNoExtension)
        {
            this.crack = crack;
            this.model = model;
            this.pathNoExtension = pathNoExtension;
            this.triangulator = new CartesianTriangulator();
        }

        public void WriteOutputData(IDofOrderer dofOrderer, Vector freeDisplacements, Vector constrainedDisplacements, int step)
        {
            // TODO: guess initial capacities from previous steps or from the model
            var allPoints = new List<VtkPoint2D>();
            var allCells = new List<VtkCell2D>();
            var displacements = new List<Vector2>();
            int pointCounter = 0;

            foreach (XContinuumElement2D element in model.Elements)
            {
                Vector elementDisplacements = dofOrderer.ExtractDisplacementVectorOfElementFromGlobal(element,
                    freeDisplacements, constrainedDisplacements);
                bool isEnriched = IsEnriched(element);
                if (!isEnriched)
                {
                    // Mesh
                    var cellPoints = new VtkPoint2D[element.Nodes.Count];
                    for (int p = 0; p < cellPoints.Length; ++p)
                    {
                        cellPoints[p] = new VtkPoint2D(pointCounter++, element.Nodes[p]);
                        allPoints.Add(cellPoints[p]);
                    }
                    allCells.Add(new VtkCell2D(VtkCell2D.CellTypeCodes[element.ElementType], cellPoints));

                    // Displacements
                    for (int p = 0; p < cellPoints.Length; ++p)
                    {
                        displacements.Add(Vector2.Create(elementDisplacements[2 * p], elementDisplacements[2 * p + 1]));
                    }
                }
                else
                {
                    IInverseMapping2D inverseMapping = element.Interpolation.CreateInverseMappingFor(element.Nodes);
                    Vector standardDisplacements = dofOrderer.ExtractDisplacementVectorOfElementFromGlobal(element,
                        freeDisplacements, constrainedDisplacements);
                    Vector enrichedDisplacements =
                        dofOrderer.ExtractEnrichedDisplacementsOfElementFromGlobal(element, freeDisplacements);

                    

                    // Triangulate and then operate on each triangle
                    SortedSet<ICartesianPoint2D> triangleVertices = crack.FindTriangleVertices(element);
                    IReadOnlyList<TriangleCartesian2D> triangles = triangulator.CreateMesh(triangleVertices);

                    foreach (TriangleCartesian2D triangle in triangles)
                    {
                        // Mesh
                        int numCellPoints = 3;
                        var cellPoints = new VtkPoint2D[numCellPoints];
                        for (int p = 0; p < numCellPoints; ++p)
                        {
                            cellPoints[p] = 
                                new VtkPoint2D(pointCounter++, new CartesianPoint2D(triangle.Vertices[p].Coordinates));
                            allPoints.Add(cellPoints[p]);
                        }
                        allCells.Add(new VtkCell2D(triangleVtkCode, cellPoints));

                        // Displacements, strains and stresses are not defined on the crack, thus they must be evaluated at GPs   
                        // and extrapolated to each point of interest. However how should I choose the Gauss points? Here I take 
                        // the Gauss points of the subtriangles.
                        IGaussPointSystem extrapolation = new Tri3GPSystem();

                        // Find the Gauss points of the triangle in the natural system of the element
                        var triangleNodesNatural = new INaturalPoint2D[numCellPoints];
                        for (int p = 0; p < numCellPoints; ++p)
                        {
                            triangleNodesNatural[p] = inverseMapping.TransformCartesianToNatural(cellPoints[p]);
                        }
                        NaturalPoint2D[] triangleGPsNatural = 
                            FindTriangleGPsNatural(triangleNodesNatural, extrapolation.GaussPoints);

                        // Find the field values at the Gauss points of the triangle (their coordinates are in the natural 
                        // system of the element)
                        var displacementsAtGPs = new Vector2[triangleGPsNatural.Length];
                        for (int gp = 0; gp < triangleGPsNatural.Length; ++gp)
                        {
                            EvaluatedInterpolation2D evalInterpol =
                                    element.Interpolation.EvaluateAt(element.Nodes, triangleGPsNatural[gp]);
                            displacementsAtGPs[gp] = element.CalculateDisplacementField(triangleGPsNatural[gp],
                                evalInterpol, standardDisplacements, enrichedDisplacements);
                        }

                        // Extrapolate the field values to the triangle nodes. We need their coordinates in the auxiliary 
                        // system of the triangle. We could use the inverse interpolation of the triangle to map the natural 
                        // (element local) coordinates of the nodes to the auxiliary system of the triangle. Fortunately they 
                        // can be accessed by the extrapolation object directly.
                        IReadOnlyList<Vector2> displacementsAtTriangleNodes = 
                            extrapolation.ExtrapolateVectorFromGaussPointsToNodes(displacementsAtGPs);
                        for (int p = 0; p < numCellPoints; ++p)
                        {
                            displacements.Add(displacementsAtTriangleNodes[p]);
                        }
                    }
                }
            }

            using (var writer = new VtkFileWriter2D($"{pathNoExtension}_{step}.vtk"))
            {
                writer.WriteMesh(allPoints, allCells);
                writer.WriteVector2DField("displacements", displacements);
            }
        }

        //TODO: this should be available to all XFEM classes.
        private bool IsEnriched(XContinuumElement2D element)
        {
            foreach (XNode2D node in element.Nodes)
            {
                if (node.EnrichmentItems.Count != 0)
                {
                    return true;
                }
            }
            return false;
        }

        private NaturalPoint2D[] FindTriangleGPsNatural(IReadOnlyList<INaturalPoint2D> triangleNodesNatural, 
            IReadOnlyList<GaussPoint2D> triangleGPsAuxiliary)
        {
            //TODO: consider removing all this type safety in the coordinate systems or use generics
            // Copy the triangle nodes' natural coordinates to pseudo-cartesian nodes, so that interpolations can be used. 
            var nodesPseudoCartesian = new Node2D[3];
            for (int i = 0; i < 3; ++i)
            {
                nodesPseudoCartesian[i] = new Node2D(i, triangleNodesNatural[i].Xi, triangleNodesNatural[i].Eta);
            }

            var triangleGPsNatural = new NaturalPoint2D[triangleGPsAuxiliary.Count];
            for (int i = 0; i < triangleGPsAuxiliary.Count; ++i)
            {
                // Map the triangle's Gauss points from the auxiliary system (natural system of the triangle) to the natural system
                // (pseudo Cartesian system of the triangle)
                ICartesianPoint2D pseudoCartesian = IsoparametricInterpolation2D.Tri3.TransformNaturalToCartesian(
                    nodesPseudoCartesian, triangleGPsAuxiliary[i]);

                // Copy the pseudo cartesian coordinates to natural system
                triangleGPsNatural[i] = new NaturalPoint2D(pseudoCartesian.X, pseudoCartesian.Y);
            }

            return triangleGPsNatural;
        }
    }
}
