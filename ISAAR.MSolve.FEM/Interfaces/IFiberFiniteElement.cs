using ISAAR.MSolve.Materials.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IFiberFiniteElement : IFiniteElement_v2
    {
        IFiberFiniteElementMaterial Material { get; }
        IList<IFiber> Fibers { get; }
    }
}
