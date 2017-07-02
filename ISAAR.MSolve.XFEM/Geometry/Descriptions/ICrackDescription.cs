using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Enrichments.Items.CrackTip;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Shapes;
using ISAAR.MSolve.XFEM.Geometry.Triangulation;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Interpolation;


namespace ISAAR.MSolve.XFEM.Geometry.Descriptions
{
    // TODO: this only works for cracks with a single tip
    interface ICrackDescription
    {
        /// TODO: an interface is needed for TipSystems. Then the explicit (global, local, polar) systems or the level 
        /// sets could be used for the transformations (points, vectors, derivatives)
        TipCoordinateSystem TipSystem { get; } 
        ICartesianPoint2D CrackTip { get; }
        List<XContinuumElement2D> TipElements { get; }

        void InitializeGeometry(ICartesianPoint2D crackMouth, ICartesianPoint2D crackTip);
        void UpdateGeometry(double localGrowthAngle, double growthLength); // Perhaps the global angle should be passed in

        double SignedDistanceOf(XNode2D node);
        double SignedDistanceOf(INaturalPoint2D point, IReadOnlyList<XNode2D> elementNodes,
             EvaluatedInterpolation2D interpolation);
        Tuple<double, double> SignedDistanceGradientThrough(INaturalPoint2D point, IReadOnlyList<XNode2D> elementNodes,
             EvaluatedInterpolation2D interpolation); // Not necessary. Should be removed.

        IReadOnlyList<TriangleCartesian2D> TriangulateAreaOf(XContinuumElement2D element);
        void UpdateEnrichments();

        //PolarPoint2D ToPolar(INaturalPoint2D point, IReadOnlyList<XNode2D> elementNodes,
        //     EvaluatedInterpolation2D interpolation);
    }
}
