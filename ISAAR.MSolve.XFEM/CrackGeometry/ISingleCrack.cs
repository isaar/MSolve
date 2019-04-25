using System.Collections.Generic;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.XFEM.CrackGeometry.HeavisideSingularityResolving;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.CrackGeometry
{
    /// <summary>
    /// A crack that has one tip point and: a) a mouth point (exterior crack), b) another tip point (interior crack) or 
    /// c) a junction point (branch of a crack tree)
    /// </summary>
    public interface ISingleCrack: ICrackDescription
    {
        CrackBodyEnrichment2D CrackBodyEnrichment { get; }
        CrackTipEnrichments2D CrackTipEnrichments { get; }

        IHeavisideSingularityResolver SingularityResolver { get; }

        double SignedDistanceOf(XNode node);
        double SignedDistanceOf(NaturalPoint point, XContinuumElement2D element,
            EvalInterpolation2D interpolation);

        //CartesianPoint GetCrackTip(CrackTipPosition tipPosition);
        /// TODO: an interface is needed for TipSystems. Then the explicit (global, local, polar) systems or the level 
        /// sets could be used for the transformations (points, vectors, derivatives)
        //TipCoordinateSystem GetTipSystem(CrackTipPosition tipPosition);
        //IReadOnlyList<XContinuumElement2D> GetTipElements(CrackTipPosition tipPosition);

        void InitializeGeometry(PolyLine2D initialCrack);
        SortedSet<CartesianPoint> FindTriangleVertices(XContinuumElement2D element);
    }
}
