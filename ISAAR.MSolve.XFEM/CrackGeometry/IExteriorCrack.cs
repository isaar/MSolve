using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Shapes;
using ISAAR.MSolve.XFEM.Geometry.Triangulation;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Interpolation;


namespace ISAAR.MSolve.XFEM.CrackGeometry
{ 
    // TODO: this only works for cracks with a single tip
    interface IExteriorCrack: ICrackGeometry
    {
        void InitializeGeometry(ICartesianPoint2D crackMouth, ICartesianPoint2D crackTip);
        void UpdateGeometry(double localGrowthAngle, double growthLength); // Perhaps the global angle should be passed in
        
        //PolarPoint2D ToPolar(INaturalPoint2D point, IReadOnlyList<XNode2D> elementNodes,
        //     EvaluatedInterpolation2D interpolation);
    }
}
