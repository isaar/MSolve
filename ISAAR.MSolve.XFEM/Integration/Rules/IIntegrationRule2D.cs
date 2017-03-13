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
        /// <summary>
        /// The integration points are immutable and can be used across many elements. 
        /// The only need to store them for each element is caching. 
        /// However standard integration rules (e.g Gauss-Legendre) for regular natural coordinate systems have already
        /// precached them and the generation process is inexpensive.
        /// </summary>
        /// <returns></returns>
        IReadOnlyList<GaussPoint2D> GenerateIntegrationPoints();
    }
}
