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
    // Perhaps I can pass the element instance (actually the element will pass "this") to the GetIntegrationPoints()
    // and Update(). Alternatively I could ensure that the element passed in the constructor is not used unti the GetIntegrationPoints is called.

    /// <summary>
    /// Classes that handle storing and updating the integration rules, integration points and materials of each 
    /// element. The element only  needs to access the integration point-material pairs form this class. The 
    /// integration point generation is delegated to a suitable integration rule usually.
    /// </summary>    
    interface IIntegrationStrategy2D
    {
        /// <summary>
        /// Keep in mind that in XFEM this will be called twice: once for the std part of the stiffness matrix and once 
        /// for the enriched parts. Thus it would be beneficial to cache the results if the calculation is expensive.
        /// </summary>
        /// <returns></returns>
        IEnumerable<Tuple<GaussPoint2D, IFiniteElementMaterial2D>> GetIntegrationPointsAndMaterials();

        // TODO: The update should be done for some element, rule and new condition. 
        // The first 2 could be fields, but the new condition should be passed as an argument.
        // Perhaps this argument should be generic to facilitate different problem types.
        void Update(); 
    }

    /// <summary>
    /// TODO: Find a better design. For now this will simplify element construction. The problem is that the element 
    /// type class stores the integration strategy as a field in its constructor 
    /// and the integration strategy class needs stuff from the element type in its constructor. 
    /// Instead of having mutators for the fields of at least 1 of the 2 classes, I use this instead. 
    /// The whole element construction design needs refactoring though.
    /// </summary>
    interface IIntegrationStrategyFactory2D
    {
        /// <summary>
        /// WARNING: This method will probably try to access the fields of <see cref="ContinuumElement2D"/>. 
        /// Thus they must be readied in <see cref="ContinuumElement2D"/>'s constructor before this method used.
        /// 1) I really hate this dependency over the order of execution of methods in 2 different classes. 
        /// 2) Also there is an insidious bug, where this method results in calling a member of 
        /// <see cref="ContinuumElement2D"/> that is not accessible in the constructor of 
        /// <see cref="ContinuumElement2D"/>, but in the constructor of one of its derived classes. Spaghetti 
        /// inheritance code.
        /// </summary>
        /// <param name="elementType"></param>
        /// <returns></returns>
        IIntegrationStrategy2D CreateStrategy(ContinuumElement2D elementType);
    }
}
