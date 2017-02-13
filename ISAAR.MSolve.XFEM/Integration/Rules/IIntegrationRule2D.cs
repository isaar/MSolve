using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Integration.Points;

namespace ISAAR.MSolve.XFEM.Integration.Rules
{
    interface IIntegrationRule2D
    {
        IReadOnlyList<GaussPoint2D> GenerateIntegrationPoints();
    }
}
