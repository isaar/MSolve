using ISAAR.MSolve.Materials.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IFiberFiniteElement : IFiniteElement
    {
        IFiberFiniteElementMaterial Material { get; }
        IList<IFiber> Fibers { get; }
    }
}
