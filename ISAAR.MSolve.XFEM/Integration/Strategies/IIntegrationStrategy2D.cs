using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.XFEM.Integration.Strategies
{
    // Not sure if generics are needed. I could just have different interfaces for different element types. 
    // It would make sense if some of the actual integration strategy implementations would work for many element types.

    /// <summary>
    /// Algorithms for complex integration rules for specific finite element types. These need the data from each 
    /// finite element to generate integration points for use only by that finite element. 
    /// They typically make use of the standard quadrature rules.
    /// </summary>
    /// <typeparam name="TElement"></typeparam>
    interface IIntegrationStrategy2D<TElement>
    {
        IReadOnlyList<GaussPoint2D> GenerateIntegrationPoints(TElement element);
    }
}
