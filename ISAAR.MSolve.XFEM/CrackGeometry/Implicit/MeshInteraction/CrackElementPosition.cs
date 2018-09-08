using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.XFEM.CrackGeometry.Implicit.MeshInteraction
{
    /// <summary>
    /// Represents the type of enrichment that will be applied to all nodes of the element. In LSM with linear 
    /// interpolations, an element enriched with tip functions does not need to be enriched with Heaviside too. 
    /// This is because, even if there are kinks inside the element, the linear interpolation cannot reproduce them.
    /// </summary>
    enum CrackElementPosition
    {
        Irrelevant, Intersected, ContainsTip
    }
}
