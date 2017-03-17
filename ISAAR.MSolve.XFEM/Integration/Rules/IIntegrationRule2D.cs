using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Integration.Points;

namespace ISAAR.MSolve.XFEM.Integration.Rules
{
    /// <summary>
    /// Algorithms for complex integration rules for specific finite element types. These need the data from each 
    /// finite element to generate integration points for use only by that finite element. 
    /// They typically make use of the standard quadrature rules.
    /// </summary>
    /// <returns></returns>
    interface IIntegrationRule2D<TElement> //TODO: add type constraints or use a generic FE interface and remove the generics
    {
        IReadOnlyList<GaussPoint2D> GenerateIntegrationPoints(TElement element);
    }
}
