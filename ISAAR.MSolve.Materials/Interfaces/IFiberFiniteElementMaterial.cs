using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Materials.Interfaces
{
    public interface IFiberFiniteElementMaterial : IFiniteElementMaterial
    {
        IList<IFiberMaterial> FiberMaterials { get; }
    }
}
