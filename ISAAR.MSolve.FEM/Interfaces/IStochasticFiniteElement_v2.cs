using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IStochasticFiniteElement_v2 : IFiniteElement_v2
    {
        IStochasticCoefficientsProvider CoefficientsProvider { get; set; }
    }
}
