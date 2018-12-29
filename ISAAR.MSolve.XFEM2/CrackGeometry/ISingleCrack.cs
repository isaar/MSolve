using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.CrackGeometry.HeavisideSingularityResolving;
using ISAAR.MSolve.XFEM.CrackPropagation;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Geometry.Shapes;
using ISAAR.MSolve.XFEM.Geometry.Triangulation;
using ISAAR.MSolve.XFEM.Interpolation;

namespace ISAAR.MSolve.XFEM.CrackGeometry
{
    /// <summary>
    /// A crack that has one tip point and: a) a mouth point (exterior crack), b) another tip point (interior crack) or 
    /// c) a junction point (branch of a crack tree)
    /// </summary>
    interface ISingleCrack: ICrackDescription
    {
        CrackBodyEnrichment2D CrackBodyEnrichment { get; }
        CrackTipEnrichments2D CrackTipEnrichments { get; }

        IHeavisideSingularityResolver SingularityResolver { get; }

        double SignedDistanceOf(XNode2D node);
        double SignedDistanceOf(INaturalPoint2D point, XContinuumElement2D element,
            EvaluatedInterpolation2D interpolation);

        //ICartesianPoint2D GetCrackTip(CrackTipPosition tipPosition);
        /// TODO: an interface is needed for TipSystems. Then the explicit (global, local, polar) systems or the level 
        /// sets could be used for the transformations (points, vectors, derivatives)
        //TipCoordinateSystem GetTipSystem(CrackTipPosition tipPosition);
        //IReadOnlyList<XContinuumElement2D> GetTipElements(CrackTipPosition tipPosition);

        void InitializeGeometry(PolyLine2D initialCrack);
        SortedSet<ICartesianPoint2D> FindTriangleVertices(XContinuumElement2D element);
    }
}
