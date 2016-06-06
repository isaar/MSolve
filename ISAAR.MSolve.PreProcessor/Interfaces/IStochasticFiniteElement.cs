using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.PreProcessor.Interfaces
{
    public interface IStochasticFiniteElement : IFiniteElement
    {
        IStochasticCoefficientsProvider CoefficientsProvider { get; set; }
    }
}
