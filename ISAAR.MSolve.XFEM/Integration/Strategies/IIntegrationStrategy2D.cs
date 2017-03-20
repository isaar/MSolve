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
    /// Classes that handle storing and updating the integration rules, integration points and materials of each 
    /// element. The element only  needs to access the integration point-material pairs from this class. The 
    /// integration point generation is delegated to a suitable integration rule usually.
    /// </summary>
    /// </summary>
    /// <typeparam name="TElement"></typeparam>
    interface IIntegrationStrategy2D<TElement>
    {
        // TODO: Determine the semantics of this method: is it expected to do significant work or 
        // is it just an accessor (in which case it should be a property)
        // Notes for implementation: Keep in mind that in XFEM this will be called twice: once for the std part of  
        // the stiffness matrix and once for the enriched parts. Thus it would be beneficial to cache the results if 
        // the calculation is expensive. Actually most non trivial integration strategies need to cache the Gauss
        // points and materials anyway.
        IReadOnlyDictionary<GaussPoint2D, IFiniteElementMaterial2D> GetIntegrationPointsAndMaterials(TElement element);

        // TODO: The update should be done for some element, rule and new condition. 
        // The first 2 could be fields, but the new condition should be passed as an argument.
        // Perhaps this argument should be generic to facilitate different problem types.
        void Update(TElement element); 
    }
}
