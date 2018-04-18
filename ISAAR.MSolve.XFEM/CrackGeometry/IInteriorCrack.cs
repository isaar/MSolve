using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Elements;

namespace ISAAR.MSolve.XFEM.CrackGeometry
{
    interface IInteriorCrack: ICrackGeometry
    {
        void InitializeGeometry(ICartesianPoint2D startTip, ICartesianPoint2D endTip);
        void UpdateGeometry(double localGrowthAngleStart, double growthLengthStart,
            double localGrowthAngleEnd, double growthLengthEnd); // Perhaps the global angle should be passed in
    }
}
