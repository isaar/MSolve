using System.Collections.Generic;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.XFEM.CrackGeometry
{
    // TODO: this only works for cracks with a single tip
    interface IExteriorCrack: ISingleCrack
    {
        IReadOnlyList<CartesianPoint> CrackPath { get; }

        //TODO: remove it. It is obsolete and should be handled by ICrackGeometry.InitializeGeometry(PolyLine2D initialCrack)
        void InitializeGeometry(CartesianPoint crackMouth, CartesianPoint crackTip);

        //TODO: remove it. It is obsolete and should be handled by ICrackGeometry.Propagate()
        void UpdateGeometry(double localGrowthAngle, double growthLength); // Perhaps the global angle should be passed in
        
        //PolarPoint2D ToPolar(NaturalPoint point, IReadOnlyList<XNode2D> elementNodes,
        //     EvaluatedInterpolation2D interpolation);
    }
}
