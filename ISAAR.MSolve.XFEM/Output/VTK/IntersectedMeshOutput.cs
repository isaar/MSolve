using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.Integration;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation;
using ISAAR.MSolve.FEM.Interpolation.Inverse;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging.VTK;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

//TODO: the extrapolations from Gauss points to nodes of the subtriangles does not work near the tip, where the displacement 
//      field contains interpolated tip functions. This could be fixed (for the subtriangle nodes only) by basing the 
//      extrapolation on the XFEM interpolation instead of the standard Tri3.
//TODO: Can I cheat and assign H(x)= +-1 to points on the crack, depending on which side of the crack the triangle lies? That 
//      way I could avoid the need to extrapolate the displacements form Gauss points. Would it also suffice for stresses or do
//      I have to extrapolate?
//TODO: refactor the huge method and break it down to smaller ones.
//TODO: apply averaging at standard nodes
//TODO: plot the principal stresses, not only S11, S12, S22
namespace ISAAR.MSolve.XFEM.Output.VTK
{
    class IntersectedMeshOutput: IXfemOutput
    {
        private readonly ICrackDescription crackGeometry;
        private readonly Model2D_old model;
        private readonly string pathNoExtension;
        private readonly Triangulator2D<CartesianPoint> triangulator;

        public IntersectedMeshOutput(Model2D_old model, ICrackDescription crackGeometry, string pathNoExtension)
        {
            this.crackGeometry = crackGeometry;
            this.model = model;
            this.pathNoExtension = pathNoExtension;
            this.triangulator = new Triangulator2D<CartesianPoint>((x, y) => new CartesianPoint(x, y));
        }

        public void WriteOutputData(IDofOrderer dofOrderer, Vector freeDisplacements, Vector constrainedDisplacements, int step)
        {
            // TODO: guess initial capacities from previous steps or from the model
            var allPoints = new List<VtkPoint>();
            var allCells = new List<VtkCell>();
            var displacements = new List<double[]>();
            var strains = new List<Tensor2D>();
            var stresses = new List<Tensor2D>();
            int pointCounter = 0;

            foreach (XContinuumElement2D element in model.Elements)
            {
                Vector standardDisplacements = dofOrderer.ExtractDisplacementVectorOfElementFromGlobal(element,
                    freeDisplacements, constrainedDisplacements);
                Vector enrichedDisplacements =
                    dofOrderer.ExtractEnrichedDisplacementsOfElementFromGlobal(element, freeDisplacements);
                bool mustTriangulate = MustBeTriangulated(element, out ISingleCrack intersectingCrack);

                if (!mustTriangulate)
                {
                    // Mesh
                    var cellPoints = new VtkPoint[element.Nodes.Count];
                    for (int p = 0; p < cellPoints.Length; ++p)
                    {
                        cellPoints[p] = new VtkPoint(pointCounter++, element.Nodes[p]);
                        allPoints.Add(cellPoints[p]);
                    }
                    allCells.Add(new VtkCell(element.CellType, cellPoints));

                    // Displacements
                    for (int p = 0; p < cellPoints.Length; ++p)
                    {
                        displacements.Add(new double[] { standardDisplacements[2 * p], standardDisplacements[2 * p + 1]});
                    }

                    // Strains and stresses at Gauss points of element
                    // WARNING: do not use the quadrature object, since GPs are sorted differently.
                    IReadOnlyList<GaussPoint> gaussPoints = element.GaussPointExtrapolation.Quadrature.IntegrationPoints;
                    var strainsAtGPs = new Tensor2D[gaussPoints.Count];
                    var stressesAtGPs = new Tensor2D[gaussPoints.Count];
                    for (int gp = 0; gp < gaussPoints.Count; ++gp)
                    {
                        EvalInterpolation2D evalInterpol =
                                element.Interpolation.EvaluateAllAt(element.Nodes, gaussPoints[gp]);
                        (Tensor2D strain, Tensor2D stress) = ComputeStrainStress(element, gaussPoints[gp],
                                evalInterpol, standardDisplacements, enrichedDisplacements);
                        strainsAtGPs[gp] = strain;
                        stressesAtGPs[gp] = stress;
                    }

                    // Extrapolate strains and stresses to element nodes. This is exact, since the element is not enriched
                    IReadOnlyList<Tensor2D> strainsAtNodes = element.GaussPointExtrapolation.
                        ExtrapolateTensorFromGaussPointsToNodes(strainsAtGPs, element.Interpolation);
                    IReadOnlyList<Tensor2D> stressesAtNodes = element.GaussPointExtrapolation.
                        ExtrapolateTensorFromGaussPointsToNodes(stressesAtGPs, element.Interpolation);
                    for (int p = 0; p < cellPoints.Length; ++p)
                    {
                        strains.Add(strainsAtNodes[p]);
                        stresses.Add(stressesAtNodes[p]);
                    }
                }
                else
                {
                    // Triangulate and then operate on each triangle
                    SortedSet<CartesianPoint> triangleVertices = intersectingCrack.FindTriangleVertices(element);
                    IReadOnlyList<Triangle2D<CartesianPoint>> triangles = triangulator.CreateMesh(triangleVertices);

                    foreach (Triangle2D<CartesianPoint> triangle in triangles)
                    {
                        // Mesh
                        int numTriangleNodes = 3;
                        var cellPoints = new VtkPoint[numTriangleNodes];
                        for (int p = 0; p < numTriangleNodes; ++p)
                        {
                            CartesianPoint point = triangle.Vertices[p];
                            cellPoints[p] = new VtkPoint(pointCounter++, point.X, point.Y, point.Z);
                            allPoints.Add(cellPoints[p]);
                        }
                        allCells.Add(new VtkCell(CellType.Tri3, cellPoints));

                        // Displacements, strains and stresses are not defined on the crack, thus they must be evaluated at GPs   
                        // and extrapolated to each point of interest. However how should I choose the Gauss points? Here I take 
                        // the Gauss points of the subtriangles.
                        IGaussPointExtrapolation2D extrapolation = ExtrapolationGaussTriangular3Points.UniqueInstance;
                        IIsoparametricInterpolation2D interpolation = InterpolationTri3.UniqueInstance;

                        // Find the Gauss points of the triangle in the natural system of the element
                        IInverseInterpolation2D inverseMapping = element.Interpolation.CreateInverseMappingFor(element.Nodes);
                        var triangleNodesNatural = new NaturalPoint[numTriangleNodes];
                        for (int p = 0; p < numTriangleNodes; ++p)
                        {
                            triangleNodesNatural[p] = inverseMapping.TransformPointCartesianToNatural(cellPoints[p]);
                        }
                        NaturalPoint[] triangleGPsNatural = 
                            FindTriangleGPsNatural(triangleNodesNatural, extrapolation.Quadrature.IntegrationPoints);

                        // Find the field values at the Gauss points of the triangle (their coordinates are in the natural 
                        // system of the element)
                        var displacementsAtGPs = new double[triangleGPsNatural.Length][];
                        var strainsAtGPs = new Tensor2D[triangleGPsNatural.Length];
                        var stressesAtGPs = new Tensor2D[triangleGPsNatural.Length];
                        for (int gp = 0; gp < triangleGPsNatural.Length; ++gp)
                        {
                            EvalInterpolation2D evalInterpol =
                                    element.Interpolation.EvaluateAllAt(element.Nodes, triangleGPsNatural[gp]);
                            displacementsAtGPs[gp] = element.CalculateDisplacementField(triangleGPsNatural[gp],
                                evalInterpol, standardDisplacements, enrichedDisplacements).CopyToArray();
                            (Tensor2D strain, Tensor2D stress) = ComputeStrainStress(element, triangleGPsNatural[gp],
                                evalInterpol, standardDisplacements, enrichedDisplacements);
                            strainsAtGPs[gp] = strain;
                            stressesAtGPs[gp] = stress;
                        }

                        // Extrapolate the field values to the triangle nodes. We need their coordinates in the auxiliary 
                        // system of the triangle. We could use the inverse interpolation of the triangle to map the natural 
                        // (element local) coordinates of the nodes to the auxiliary system of the triangle. Fortunately they 
                        // can be accessed by the extrapolation object directly.
                        IReadOnlyList<double[]> displacementsAtTriangleNodes =
                            extrapolation.ExtrapolateVectorFromGaussPointsToNodes(displacementsAtGPs, interpolation);
                        IReadOnlyList<Tensor2D> strainsAtTriangleNodes =
                            extrapolation.ExtrapolateTensorFromGaussPointsToNodes(strainsAtGPs, interpolation);
                        IReadOnlyList<Tensor2D> stressesAtTriangleNodes =
                            extrapolation.ExtrapolateTensorFromGaussPointsToNodes(stressesAtGPs, interpolation);
                        for (int p = 0; p < numTriangleNodes; ++p)
                        {
                            displacements.Add(displacementsAtTriangleNodes[p]);
                            strains.Add(strainsAtTriangleNodes[p]);
                            stresses.Add(stressesAtTriangleNodes[p]);
                        }
                    }
                }
            }

            using (var writer = new VtkFileWriter($"{pathNoExtension}_{step}.vtk"))
            {
                writer.WriteMesh(allPoints, allCells);
                writer.WriteVector2DField("displacement", displacements);
                writer.WriteTensor2DField("strain", strains);
                writer.WriteTensor2DField("stress", stresses);
            }
        }

        private bool MustBeTriangulated(XContinuumElement2D element, out ISingleCrack intersectingCrack)
        {
            bool isTipEnrichedOrBlending = false;


            var enrichments = new HashSet<IEnrichmentItem2D>();
            foreach (XNode node in element.Nodes)
            {
                foreach (IEnrichmentItem2D enrichment in node.EnrichmentItems.Keys)
                {
                    enrichments.Add(enrichment);
                    if (enrichment is CrackTipEnrichments2D) isTipEnrichedOrBlending = true;
                }
            }

            if (enrichments.Count == 0)
            {
                intersectingCrack = null;
                return false;
            }

            var singleCracks = new List<ISingleCrack>();
            foreach (ISingleCrack crack in crackGeometry.SingleCracks)
            {
                if (enrichments.Contains(crack.CrackBodyEnrichment) || enrichments.Contains(crack.CrackTipEnrichments))
                {
                    singleCracks.Add(crack);
                }
            }

            if (singleCracks.Count == 1)
            {
                if (isTipEnrichedOrBlending)
                {
                    intersectingCrack = singleCracks[0];
                    return true;
                }
                else
                {
                    int positiveNodes = 0;
                    int negativeNodes = 0;
                    foreach (XNode node in element.Nodes)
                    {
                        if (singleCracks[0].SignedDistanceOf(node) > 0.0) ++positiveNodes;
                        else if (singleCracks[0].SignedDistanceOf(node) < 0.0) ++negativeNodes;
                        else
                        {
                            ++positiveNodes;
                            ++negativeNodes;
                        }
                    }

                    if ((positiveNodes > 0) && (negativeNodes > 0))
                    {
                        intersectingCrack = singleCracks[0];
                        return true;
                    }
                    else
                    {
                        intersectingCrack = null; ;
                        return false;
                    }
                }
            }
            else
            {
                throw new Exception($"This element must be intersected by exactly 1 crack, but {singleCracks.Count} were found");
            }
        }

        private NaturalPoint[] FindTriangleGPsNatural(IReadOnlyList<NaturalPoint> triangleNodesNatural, 
            IReadOnlyList<GaussPoint> triangleGPsAuxiliary)
        {
            //TODO: consider removing all this type safety in the coordinate systems or use generics
            // Copy the triangle nodes' natural coordinates to pseudo-cartesian nodes, so that interpolations can be used. 
            var nodesPseudoCartesian = new XNode[3];
            for (int i = 0; i < 3; ++i)
            {
                nodesPseudoCartesian[i] = new XNode(i, triangleNodesNatural[i].Xi, triangleNodesNatural[i].Eta);
            }

            var triangleGPsNatural = new NaturalPoint[triangleGPsAuxiliary.Count];
            for (int i = 0; i < triangleGPsAuxiliary.Count; ++i)
            {
                // Map the triangle's Gauss points from the auxiliary system (natural system of the triangle) to the natural system
                // (pseudo Cartesian system of the triangle)
                CartesianPoint pseudoCartesian = InterpolationTri3.UniqueInstance.TransformNaturalToCartesian(
                    nodesPseudoCartesian, triangleGPsAuxiliary[i]);

                // Copy the pseudo cartesian coordinates to natural system
                triangleGPsNatural[i] = new NaturalPoint(pseudoCartesian.X, pseudoCartesian.Y);
            }

            return triangleGPsNatural;
        }

        private (Tensor2D strain, Tensor2D stress) ComputeStrainStress(XContinuumElement2D element, NaturalPoint gaussPoint,
            EvalInterpolation2D evaluatedInterpolation, Vector standardNodalDisplacements,
            Vector enrichedNodalDisplacements)
        {
            Matrix constitutive =
                element.Material.CalculateConstitutiveMatrixAt(gaussPoint, evaluatedInterpolation);
            Matrix2by2 displacementGradient = element.CalculateDisplacementFieldGradient(gaussPoint, evaluatedInterpolation,
                standardNodalDisplacements, enrichedNodalDisplacements);

            double strainXX = displacementGradient[0, 0];
            double strainYY = displacementGradient[1, 1];
            double strainXYtimes2 = displacementGradient[0, 1] + displacementGradient[1, 0];

            double stressXX = constitutive[0, 0] * strainXX + constitutive[0, 1] * strainYY;
            double stressYY = constitutive[1, 0] * strainXX + constitutive[1, 1] * strainYY;
            double stressXY = constitutive[2, 2] * strainXYtimes2;

            return (new Tensor2D(strainXX, strainYY, 0.5 * strainXYtimes2), new Tensor2D(stressXX, stressYY, stressXY));
        }
    }
}
